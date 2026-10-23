"""Controlled vocabulary helpers for raw2isa.

This module provides a minimal parser for OBO files and a set of static
mapping tables translating Thermo RAW metadata (scan filter strings, sample
types, ...) into PSI-MS controlled vocabulary terms.

The PSI-MS ontology file shipped with raw2isa in the ``ontologies`` package
is only used to resolve *instrument models* by name, since the model string
found in a RAW file can take many forms; all other terms are resolved from
the static tables below, which were manually curated from the PSI-MS CV.

License:
    GNU General Public License version 3.0 (GPLv3)
"""

import collections

from . import ontologies
from ._impl import resource_files

#: Base PURL used to URL-ize PSI-MS accessions in the metadata dictionary.
MS_PURL = "http://purl.obolibrary.org/obo/MS_{}"

#: PURL used to URL-ize Units of Measurement ontology accessions.
UO_PURL = "http://purl.obolibrary.org/obo/UO_{}"

# --- Static PSI-MS term tables -------------------------------------------------
# All terms were verified against the PSI-MS controlled vocabulary.

#: CV terms for the ion polarity extracted from RAW scan filter strings.
POLARITY_CV = {
    "positive": {"name": "positive scan", "accession": "MS:1000130", "ref": "MS"},
    "negative": {"name": "negative scan", "accession": "MS:1000129", "ref": "MS"},
    "alternating": {"name": "alternating scan"},
}

#: CV terms for the mass analyzer tokens found in RAW scan filter strings.
#: ``FTMS`` is ambiguous (Orbitrap vs FT-ICR) and is resolved at runtime
#: using the instrument model, see ``raw2isa.raw``.
ANALYZER_CV = {
    "ITMS": {"name": "quadrupole ion trap", "accession": "MS:1000082", "ref": "MS"},
    "SQMS": {"name": "quadrupole", "accession": "MS:1000081", "ref": "MS"},
    "TQMS": {"name": "quadrupole", "accession": "MS:1000081", "ref": "MS"},
    "TOFMS": {"name": "time-of-flight", "accession": "MS:1000084", "ref": "MS"},
    "FTMS-ORBITRAP": {"name": "orbitrap", "accession": "MS:1000484", "ref": "MS"},
    "FTMS-FTICR": {
        "name": "fourier transform ion cyclotron resonance",
        "accession": "MS:1000079",
        "ref": "MS",
    },
    "FTMS": {"name": "fourier transform mass spectrometer"},
}

#: CV terms for the ion source tokens found in RAW scan filter strings.
ION_SOURCE_CV = {
    "ESI": {
        "name": "electrospray ionization",
        "accession": "MS:1000073",
        "ref": "MS",
    },
    "NSI": {"name": "nanoelectrospray", "accession": "MS:1000398", "ref": "MS"},
    "APCI": {
        "name": "atmospheric pressure chemical ionization",
        "accession": "MS:1000070",
        "ref": "MS",
    },
    "APPI": {
        "name": "atmospheric pressure photoionization",
        "accession": "MS:1000382",
        "ref": "MS",
    },
}

#: CV terms for the spectrum representation extracted from scan filters.
REPRESENTATION_CV = {
    "centroid": {"name": "centroid spectrum", "accession": "MS:1000127", "ref": "MS"},
    "profile": {"name": "profile spectrum", "accession": "MS:1000128", "ref": "MS"},
}

#: CV terms for the content of the data file (MS level based).
DATA_FILE_CONTENT_CV = {
    "MS1": {"name": "MS1 spectrum", "accession": "MS:1000579", "ref": "MS"},
    "MSn": {"name": "MSn spectrum", "accession": "MS:1000580", "ref": "MS"},
}

#: Names of the RAW file sample types (Thermo sample table).
SAMPLE_TYPE_CV = {
    0: "Unknown",
    1: "Blank",
    2: "Standard",
    3: "Quality Control",
    4: "Sample",
    5: "Calibration Standard",
    6: "System Suitability",
}

#: Fixed CV terms for file-level metadata.
FIXED_CV = {
    "raw data file format": {
        "name": "Thermo RAW format",
        "accession": "MS:1000563",
        "ref": "MS",
    },
    "native spectrum identifier format": {
        "name": "scan number only nativeID format",
        "accession": "MS:1000776",
        "ref": "MS",
    },
    "m/z unit": {"name": "m/z", "accession": "MS:1000040", "ref": "MS"},
    "minute unit": {"name": "minute", "accession": "UO:0000031", "ref": "UO"},
    "microlitre unit": {"name": "microlitre", "accession": "UO:0000101", "ref": "UO"},
}

#: Thermo instrument manufacturer terms, used as ancestors when resolving
#: an instrument model in the ontology.
MANUFACTURER_CV = {
    "MS:1000125": {"name": "Thermo Finnigan instrument model"},
    "MS:1000483": {"name": "Thermo Fisher Scientific instrument model"},
    "MS:1000492": {"name": "Thermo Electron instrument model"},
    "MS:1000494": {"name": "Thermo Scientific instrument model"},
}

#: Terms identifying a Thermo instrument manufacturer; default used when
#: the model is not found in the ontology.
DEFAULT_MANUFACTURER = {
    "name": "Thermo Fisher Scientific instrument model",
    "accession": "MS:1000483",
    "ref": "MS",
}

#: Regex patterns matching Orbitrap-family instrument models, used to
#: disambiguate the ``FTMS`` scan filter analyzer token.
ORBITRAP_FAMILY_PATTERNS = (
    "orbitrap",
    "exactive",
    "astral",
    "exploris",
    "fusion",
    "ascend",
    "stella",
)


#: A controlled vocabulary term parsed from an OBO file.
CVTerm = collections.namedtuple("CVTerm", ["accession", "name", "parents", "synonyms", "obsolete"])


def parse_obo(path):
    """Parse an OBO file and yield its terms.

    Arguments:
        path (str): the path to the OBO file to parse.

    Yields:
        CVTerm: the next term in the ontology.
    """
    term = None
    with open(path, encoding="utf-8") as file:
        for line in file:
            line = line.rstrip("\n").rstrip("\r")
            if line == "[Term]":
                if term is not None:
                    yield CVTerm(**term)
                term = {
                    "accession": None,
                    "name": None,
                    "parents": [],
                    "synonyms": [],
                    "obsolete": False,
                }
            elif term is not None:
                if line.startswith("id: "):
                    term["accession"] = line[4:]
                elif line.startswith("name: "):
                    term["name"] = line[6:]
                elif line.startswith("is_a: "):
                    term["parents"].append(line[6:].split("!")[0].strip())
                elif line.startswith("synonym: "):
                    term["synonyms"].append(line[len("synonym: "):].split('"')[1])
                elif line.startswith("is_obsolete: true"):
                    term["obsolete"] = True
            elif line == "":
                term = None
    if term is not None:
        yield CVTerm(**term)


class PSI_MSOntology:
    """A lazily-loaded index of the PSI-MS ontology.

    Only the subset of the ontology needed to resolve Thermo instrument
    models is kept in memory: terms which are descendants of one of the
    Thermo manufacturer terms.
    """

    _instance = None

    def __init__(self):
        path = resource_files(ontologies.__name__).joinpath("psi-ms.obo")
        with open(path, encoding="utf-8"):
            self._terms = {
                term.accession: term
                for term in parse_obo(path)
                if not term.obsolete and term.name
            }
        self._by_name = {term.name.lower(): term for term in self._terms.values() if not term.obsolete}

    @classmethod
    def instance(cls):
        """Return the singleton instance of the ontology, loading it if needed."""
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def _is_descendant_of(self, accession, ancestors, _seen=None):
        seen = _seen if _seen is not None else set()
        if accession in seen:
            return False
        seen.add(accession)
        term = self._terms.get(accession)
        if term is None:
            return False
        for parent in term.parents:
            if parent in ancestors:
                return True
            if self._is_descendant_of(parent, ancestors, seen):
                return True
        return False

    def find_instrument_model(self, model):
        """Resolve a Thermo instrument model string to a CV term.

        The model string found in RAW files (e.g. ``Q Exactive HF-X Orbitrap``)
        is matched against the names and synonyms of the PSI-MS terms that are
        descendants of the Thermo manufacturer terms; the most specific (i.e.
        longest) matching term name wins.

        Arguments:
            model (str): the instrument model string.

        Returns:
            dict or None: a dictionary with the ``name``, ``accession`` and
            ``ref`` of the matching term, or None if no match was found.
        """
        if not model:
            return None

        lowered = model.lower()

        # collect matching Thermo instrument model terms, with a priority:
        # tier 0 for exact matches (name or synonym), tier 1 when the term
        # name is contained in the model string (e.g. the RAW model
        # "Q Exactive HF-X Orbitrap" contains the term "Q Exactive HF-X"),
        # and tier 2 when the model string is contained in the term name.
        matches = []

        for candidate in self._terms.values():
            if candidate.obsolete or not candidate.name:
                continue
            name = candidate.name.lower()
            tier = None
            if name == lowered or lowered in {s.lower() for s in candidate.synonyms}:
                tier = 0
            elif name in lowered:
                tier = 1
            elif lowered in name:
                tier = 2
            if tier is None:
                continue
            if not self._is_descendant_of(
                candidate.accession, set(MANUFACTURER_CV)
            ):
                continue
            matches.append((tier, candidate))

        if not matches:
            return None

        # pick the best match: lowest tier first, then the most specific
        # (i.e. longest) term name within tier 0 and 1, and the shortest
        # within tier 2
        tier, term = min(
            matches,
            key=lambda match: (
                match[0],
                -len(match[1].name) if match[0] < 2 else len(match[1].name),
            ),
        )
        return {"name": term.name, "accession": term.accession, "ref": "MS"}

    def find_manufacturer(self, model_term):
        """Resolve the manufacturer term for a given instrument model term.

        Arguments:
            model_term (dict): the instrument model term, as returned by
                ``find_instrument_model``.

        Returns:
            dict: a dictionary with the ``name``, ``accession`` and ``ref``
            of the manufacturer term (a default Thermo term if no ancestor
            matched).
        """
        if model_term is None:
            return dict(DEFAULT_MANUFACTURER)
        accession = model_term["accession"]
        seen = set()
        stack = list(self._terms.get(accession, CVTerm(None, None, [], [], False)).parents)
        while stack:
            ancestor = stack.pop()
            if ancestor in seen:
                continue
            seen.add(ancestor)
            if ancestor in MANUFACTURER_CV:
                return {
                    "name": MANUFACTURER_CV[ancestor]["name"],
                    "accession": ancestor,
                    "ref": "MS",
                }
            stack.extend(self._terms.get(ancestor, CVTerm(None, None, [], [], False)).parents)
        return dict(DEFAULT_MANUFACTURER)
