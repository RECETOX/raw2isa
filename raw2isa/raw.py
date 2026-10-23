"""Thermo RAW file metadata extractor.

This module contains the ``ThermoRawFile`` class, which uses the ``epyraw``
library (a Python wrapper around the Thermo Fisher ``RawFileReader`` .NET
library) to extract metadata from a Thermo ``.raw`` file, and returns it
as a dictionary following the same conventions as the ``mzml2isa`` package:

- a CV term: a dictionary with the ``name``, ``accession`` and ``ref`` keys,
  and optionally ``value`` and ``unit`` keys;
- a multi-value field: a dictionary with a single ``entry_list`` key;
- a plain value: a dictionary with a single ``value`` key.

All PSI-MS accessions are stored URL-ized (e.g. ``MS:1000529`` is stored as
``http://purl.obolibrary.org/obo/MS_1000529``).

License:
    GNU General Public License version 3.0 (GPLv3)
"""

import collections
import os
import re

from .cv import (
    ANALYZER_CV,
    DATA_FILE_CONTENT_CV,
    FIXED_CV,
    ION_SOURCE_CV,
    MS_PURL,
    ORBITRAP_FAMILY_PATTERNS,
    POLARITY_CV,
    REPRESENTATION_CV,
    SAMPLE_TYPE_CV,
    UO_PURL,
    PSI_MSOntology,
)

# --- Regexes for the Thermo scan filter strings --------------------------------
# Example filters:
#   "FTMS + p ESI Full ms [350.0000-1800.0000]"
#   "ITMS + c ESI d Full ms2 245.11@cid35.00 [110.000-1425.000]"
#   "TQMS - c NSI SRM MS2 245.11@cid35.00 [110.000-1425.000]"

FILTER_ANALYZER_REGEX = re.compile(
    r"^(?P<analyzer>[A-Z]{2,6}MS)\s", re.IGNORECASE
)
FILTER_POLARITY_REGEX = re.compile(r"\s(?P<polarity>[+-])\s")
FILTER_PROFILE_REGEX = re.compile(r"[+-]\s(?P<representation>[cp])\s", re.IGNORECASE)
FILTER_MSTYPE_REGEX = re.compile(r"\bms\s*(?P<level>\d*)\b", re.IGNORECASE)
FILTER_IONSOURCE_REGEX = re.compile(
    r"\s(?P<source>ESI|NSI|APCI|APPI|CI|EI|FAB|FD|GD|MALDI)\b"
)


def _urlize(term):
    """Return a copy of a CV term dict with its accession URL-ized."""
    term = dict(term)
    accession = term.get("accession")
    if accession and ":" in accession:
        prefix, local = accession.split(":", 1)
        if prefix == "MS":
            term["accession"] = MS_PURL.format(local.replace(":", "_"))
        elif prefix == "UO":
            term["accession"] = UO_PURL.format(local)
    return term


def _most_frequent(counter):
    """Return the most common value of a counter, or None if empty."""
    return counter.most_common(1)[0][0] if counter else None


class ThermoRawFile:
    """Extract metadata from a Thermo RAW file using ``epyraw``.

    Arguments:
        path (str): the local filesystem path to the RAW file.
    """

    def __init__(self, path):
        try:
            from epyraw import Device, RawFile
        except ImportError as err:
            raise ImportError(
                "the `epyraw` package is required to parse Thermo RAW files; "
                "install it with `pip install epyraw` (https://github.com/hechth/epyraw), "
                "and make sure a .NET or Mono runtime is available"
            ) from err
        self.path = path
        self._epyraw = RawFile
        self._device = Device

    @property
    def metadata(self):
        """Extract and return the metadata dictionary of the RAW file.

        Returns:
            dict: a dictionary containing the extracted metadata, with the
            same structure as the one produced by ``mzml2isa``.
        """
        RawFile = self._epyraw
        with RawFile(self.path) as raw_file:
            raw_file.select_instrument(self._device.MS, 1)
            file_header = raw_file.file_header
            run_header = raw_file.run_header
            sample_information = raw_file.sample_information
            instrument = raw_file.get_instrument_data()
            filters = self._collect_filters(raw_file)

            meta = {}
            self._extract_instrument(meta, instrument)
            self._extract_sample(meta, sample_information)
            self._extract_run(meta, run_header, file_header)
            self._extract_filters(meta, filters)
            self._extract_files(meta)

            # URL-ize all the accessions of the metadata dictionary
            for key, value in meta.items():
                if "accession" in value:
                    meta[key] = _urlize(value)
                elif "entry_list" in value:
                    meta[key] = {
                        "entry_list": [
                            _urlize(entry) if "accession" in entry else entry
                            for entry in value["entry_list"]
                        ]
                    }
                if "unit" in value:
                    meta[key]["unit"] = _urlize(value["unit"])
            return meta

    # --- Instrument ------------------------------------------------------------

    def _extract_instrument(self, meta, instrument):
        ontology = PSI_MSOntology.instance()
        model_term = ontology.find_instrument_model(instrument.model)
        if model_term is not None:
            meta["Instrument"] = model_term
        elif instrument.model:
            meta["Instrument"] = {"name": instrument.model}
        elif instrument.name:
            meta["Instrument"] = {"name": instrument.name}
        if "Instrument" in meta:
            meta["Instrument manufacturer"] = ontology.find_manufacturer(model_term)
        if instrument.serial_number:
            meta["Instrument serial number"] = {"value": instrument.serial_number}
        if instrument.software_version:
            meta["Instrument software version"] = {"value": instrument.software_version}
        if instrument.hardware_version:
            meta["Instrument hardware version"] = {"value": instrument.hardware_version}
        self._instrument_model = instrument.model or ""

    # --- Sample ------------------------------------------------------------------

    def _extract_sample(self, meta, sample_information):
        filename = os.path.basename(self.path)
        stem = os.path.splitext(filename)[0]
        name = sample_information.sample_name or stem
        meta["Sample Name"] = {"value": name}
        meta["MS Assay Name"] = {"value": name}
        sample_type = SAMPLE_TYPE_CV.get(
            int(sample_information.sample_type or 0), "Unknown"
        )
        if sample_type != "Unknown":
            meta["Sample Type"] = {"value": sample_type}
        if sample_information.injection_volume:
            meta["Injection volume"] = {
                "value": sample_information.injection_volume,
                "unit": dict(FIXED_CV["microlitre unit"]),
            }
        if sample_information.comment:
            meta["Sample comment"] = {"value": sample_information.comment}
        if sample_information.sample_id:
            meta["Sample ID"] = {"value": sample_information.sample_id}

    # --- Run ----------------------------------------------------------------------

    def _extract_run(self, meta, run_header, file_header):
        if run_header.low_mass or run_header.high_mass:
            meta["Scan m/z range"] = {
                "value": f"{run_header.low_mass}-{run_header.high_mass}",
                "unit": dict(FIXED_CV["m/z unit"]),
            }
        if run_header.start_time or run_header.end_time:
            meta["Time range"] = {
                "value": f"{run_header.start_time:.3f}-{run_header.end_time:.3f}",
                "unit": dict(FIXED_CV["minute unit"]),
            }
        meta["Number of scans"] = {
            "value": int(run_header.spectra_count or 0)
        }
        if file_header.creation_date:
            meta["Acquisition date"] = {"value": file_header.creation_date}
        if file_header.who_created_id:
            meta["Acquisition operator"] = {"value": file_header.who_created_id}

    # --- Scan filters --------------------------------------------------------------

    def _collect_filters(self, raw_file):
        """Collect the unique scan filter strings of a RAW file.

        Arguments:
            raw_file (epyraw.RawFile): the opened RAW file.

        Returns:
            list: the list of unique filter strings of the file.
        """
        filters = []
        seen = set()
        first = raw_file.run_header.first_spectrum or 1
        last = raw_file.run_header.last_spectrum or first
        for scan in range(first, last + 1):
            try:
                scan_filter = raw_file.get_filter_for_scan_number(scan)
            except Exception:  # pragma: no cover
                continue
            string = scan_filter.filter_string
            if string and string not in seen:
                seen.add(string)
                filters.append(string)
        return filters

    def _extract_filters(self, meta, filters):
        """Extract acquisition metadata from the Thermo scan filter strings.

        Arguments:
            meta (dict): the metadata dictionary to update.
            filters (list): the unique scan filter strings of the file.
        """
        polarities = collections.Counter()
        analyzers = collections.Counter()
        representations = collections.Counter()
        ion_sources = collections.Counter()
        ms_levels = set()

        for filter_string in filters:
            match = FILTER_ANALYZER_REGEX.match(filter_string)
            if match is not None:
                analyzer = match.group("analyzer").upper()
                analyzers[analyzer] += 1
            match = FILTER_POLARITY_REGEX.search(filter_string)
            if match is not None:
                polarities["positive" if match.group("polarity") == "+" else "negative"] += 1
            match = FILTER_PROFILE_REGEX.search(filter_string)
            if match is not None:
                representations[
                    "centroid" if match.group("representation").lower() == "c" else "profile"
                ] += 1
            match = FILTER_IONSOURCE_REGEX.search(filter_string)
            if match is not None:
                ion_sources[match.group("source").upper()] += 1
            for match in FILTER_MSTYPE_REGEX.finditer(filter_string):
                level = match.group("level") or "1"
                ms_levels.add(int(level) if level.isdigit() else 1)

        # polarity
        if polarities:
            if len(polarities) == 1:
                polarity = next(iter(polarities))
            else:
                polarity = "alternating"
            meta["Scan polarity"] = dict(POLARITY_CV[polarity])

        # mass analyzer, resolving ambiguous FTMS with the instrument model
        analyzer = _most_frequent(analyzers)
        if analyzer is not None:
            if analyzer == "FTMS":
                model = getattr(self, "_instrument_model", "")
                if any(pattern in model.lower() for pattern in ORBITRAP_FAMILY_PATTERNS):
                    analyzer = "FTMS-ORBITRAP"
                elif "FT" in model.upper():
                    analyzer = "FTMS-FTICR"
            term = ANALYZER_CV.get(analyzer)
            if term is not None:
                meta["Mass analyzer"] = dict(term)

        # spectrum representation
        representation = _most_frequent(representations)
        if representation is not None:
            meta["Spectrum representation"] = dict(
                REPRESENTATION_CV[representation]
            )

        # ion source
        source = _most_frequent(ion_sources)
        if source is not None:
            term = ION_SOURCE_CV.get(source)
            if term is not None:
                meta["Ion source"] = dict(term)
            else:
                meta["Ion source"] = {"name": source}

        # data file content, based on the MS levels found in the filters
        entries = []
        if ms_levels:
            if 1 in ms_levels:
                entries.append(dict(DATA_FILE_CONTENT_CV["MS1"]))
            if any(level >= 2 for level in ms_levels):
                entries.append(dict(DATA_FILE_CONTENT_CV["MSn"]))
        else:
            entries.append(dict(DATA_FILE_CONTENT_CV["MS1"]))
        meta["Data file content"] = {"entry_list": entries}

    # --- Files ---------------------------------------------------------------------

    def _extract_files(self, meta):
        filename = os.path.basename(self.path)
        meta["Raw Spectral Data File"] = {"entry_list": [{"value": filename}]}
        meta["Raw data file format"] = dict(FIXED_CV["raw data file format"])
        meta["Native spectrum identifier format"] = dict(
            FIXED_CV["native spectrum identifier format"]
        )
