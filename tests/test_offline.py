"""Offline tests of raw2isa: CV mappings, filter regexes and the ISA writer."""

import csv
import json
import os

import pytest

from raw2isa import __version__
from raw2isa.cv import ANALYZER_CV, POLARITY_CV, PSI_MSOntology
from raw2isa.isa import ISA_Tab
from raw2isa.raw import (
    FILTER_ANALYZER_REGEX,
    FILTER_IONSOURCE_REGEX,
    FILTER_MSTYPE_REGEX,
    FILTER_POLARITY_REGEX,
    FILTER_PROFILE_REGEX,
)
from raw2isa.usermeta import UserMetaLoader


def fake_meta(sample="sample1", polarity="positive scan", analyzer="FTMS-ORBITRAP"):
    """Build a metadata dictionary mimicking a parsed Thermo RAW file."""

    def urlize(term):
        term = dict(term)
        if "accession" in term and ":" in term["accession"]:
            prefix, local = term["accession"].split(":", 1)
            if prefix == "MS":
                term["accession"] = "http://purl.obolibrary.org/obo/MS_{}".format(
                    local.replace(":", "_")
                )
            elif prefix == "UO":
                term["accession"] = f"http://purl.obolibrary.org/obo/UO_{local}"
        return term

    return {
        "Sample Name": {"value": sample},
        "MS Assay Name": {"value": sample},
        "Sample Type": {"value": "Sample"},
        "Raw Spectral Data File": {"entry_list": [{"value": sample + ".raw"}]},
        "Raw data file format": urlize(
            {
                "name": "Thermo RAW format",
                "accession": "MS:1000563",
                "ref": "MS",
            }
        ),
        "Native spectrum identifier format": urlize(
            {
                "name": "scan number only nativeID format",
                "accession": "MS:1000776",
                "ref": "MS",
            }
        ),
        "Instrument": urlize(
            {
                "name": "Q Exactive Plus",
                "accession": "MS:1002857",
                "ref": "MS",
            }
        ),
        "Instrument manufacturer": urlize(
            {
                "name": "Thermo Scientific instrument model",
                "accession": "MS:1000494",
                "ref": "MS",
            }
        ),
        "Instrument serial number": {"value": "SN12345"},
        "Instrument software version": {"value": "Xcalibur 4.2"},
        "Scan polarity": dict(POLARITY_CV["positive" if polarity == "positive scan" else "negative"]),
        "Scan m/z range": {
            "value": "70.0-1050.0",
            "unit": {
                "name": "m/z",
                "accession": "http://purl.obolibrary.org/obo/MS_1000040",
                "ref": "MS",
            },
        },
        "Time range": {
            "value": "0.005-10.014",
            "unit": {
                "name": "minute",
                "accession": "http://purl.obolibrary.org/obo/UO_0000031",
                "ref": "UO",
            },
        },
        "Number of scans": {"value": 24000},
        "Mass analyzer": urlize(dict(ANALYZER_CV[analyzer])),
        "Ion source": urlize(
            {
                "name": "electrospray ionization",
                "accession": "MS:1000073",
                "ref": "MS",
            }
        ),
        "Spectrum representation": {
            "name": "profile spectrum",
            "accession": "http://purl.obolibrary.org/obo/MS_1000128",
            "ref": "MS",
        },
        "Data file content": {
            "entry_list": [
                {
                    "name": "MS1 spectrum",
                    "accession": "http://purl.obolibrary.org/obo/MS_1000579",
                    "ref": "MS",
                },
                {
                    "name": "MSn spectrum",
                    "accession": "http://purl.obolibrary.org/obo/MS_1000580",
                    "ref": "MS",
                },
            ]
        },
        "Acquisition date": {"value": "10/22/2026 11:56:14 AM"},
    }


# --- controlled vocabulary ------------------------------------------------------


@pytest.fixture(scope="session")
def ontology():
    return PSI_MSOntology.instance()


class TestOntology:
    def test_find_known_models(self, ontology):
        for model in [
            "Q Exactive",
            "Q Exactive Plus",
            "LTQ Orbitrap Classic",
            "Orbitrap Fusion",
            "TSQ Quantum",
            "LTQ",
        ]:
            term = ontology.find_instrument_model(model)
            assert term is not None, model
            assert term["accession"].startswith("MS:")
            assert term["ref"] == "MS"

    def test_find_unknown_model(self, ontology):
        assert ontology.find_instrument_model("Not An Instrument") is None
        assert ontology.find_instrument_model("") is None

    def test_find_manufacturer(self, ontology):
        term = ontology.find_instrument_model("Q Exactive Plus")
        manufacturer = ontology.find_manufacturer(term)
        assert manufacturer is not None
        assert "Thermo" in manufacturer["name"]
        assert manufacturer["accession"].startswith("MS:")

    def test_find_manufacturer_without_model(self, ontology):
        manufacturer = ontology.find_manufacturer(None)
        assert "Thermo" in manufacturer["name"]


# --- scan filter regexes ---------------------------------------------------------


class TestFilterRegexes:
    FILTERS = [
        "FTMS + p ESI Full ms [350.0000-1800.0000]",
        "ITMS + c ESI d Full ms2 245.11@cid35.00 [110.000-1425.000]",
        "TQMS - c NSI SRM MS2 245.11@cid35.00 [110.000-1425.000]",
        "- c SID ESI Full ms [50.0000-2000.0000]",
    ]

    def test_analyzer(self):
        assert FILTER_ANALYZER_REGEX.match(self.FILTERS[0]).group("analyzer") == "FTMS"
        assert FILTER_ANALYZER_REGEX.match(self.FILTERS[1]).group("analyzer") == "ITMS"
        assert FILTER_ANALYZER_REGEX.match(self.FILTERS[2]).group("analyzer") == "TQMS"
        assert FILTER_ANALYZER_REGEX.match(self.FILTERS[3]) is None

    def test_polarity(self):
        assert FILTER_POLARITY_REGEX.search(self.FILTERS[0]).group("polarity") == "+"
        assert FILTER_POLARITY_REGEX.search(self.FILTERS[2]).group("polarity") == "-"

    def test_representation(self):
        assert FILTER_PROFILE_REGEX.search(self.FILTERS[0]).group("representation") == "p"
        assert FILTER_PROFILE_REGEX.search(self.FILTERS[1]).group("representation") == "c"

    def test_mstype(self):
        assert FILTER_MSTYPE_REGEX.search(self.FILTERS[0]).group("level") == ""
        assert FILTER_MSTYPE_REGEX.search(self.FILTERS[1]).group("level") == "2"

    def test_ionsource(self):
        assert FILTER_IONSOURCE_REGEX.search(self.FILTERS[0]).group("source") == "ESI"
        assert FILTER_IONSOURCE_REGEX.search(self.FILTERS[2]).group("source") == "NSI"


# --- ISA writer -------------------------------------------------------------------


class TestISAWriter:
    @pytest.fixture
    def metalist(self):
        positive = fake_meta("sample1")
        negative = fake_meta("sample2", polarity="negative scan")
        return [positive, negative]

    @pytest.fixture
    def outdir(self, tmp_path):
        return tmp_path / "out"

    def test_write_all_files(self, metalist, outdir):
        isa_tab = ISA_Tab(str(outdir), "MTBLS0001")
        isa_tab.write(metalist, "raw")
        files = sorted(os.listdir(outdir))
        assert files == [
            "a_MTBLS0001_NEG_metabolite_profiling_mass_spectrometry.txt",
            "a_MTBLS0001_POS_metabolite_profiling_mass_spectrometry.txt",
            "i_Investigation.txt",
            "s_MTBLS0001.txt",
        ]

    def test_assay_rows(self, metalist, outdir):
        isa_tab = ISA_Tab(str(outdir), "MTBLS0001")
        isa_tab.write(metalist, "raw")
        for polarity, sample in [("POS", "sample1"), ("NEG", "sample2")]:
            path = outdir / f"a_MTBLS0001_{polarity}_metabolite_profiling_mass_spectrometry.txt"
            with open(path, newline="") as f:
                rows = list(csv.reader(f, delimiter="\t"))
            assert len(rows) == 2  # header + data
            assert len(rows[0]) == len(rows[1])
            # entry_list of the first meta has 2 entries -> 2 triplets of columns
            assert rows[0].count("Parameter Value[Data file content]") == 2
            columns = list(zip(rows[0], rows[1], strict=False))

            def cells(header, value=None, *, columns=columns):
                return [
                    data for head, data in columns if head == header and (value is None or data == value)
                ]

            def cell_after(header, offset, *, rows=rows):
                index = rows[0].index(header)
                return rows[1][index + offset]

            assert cells("Sample Name") == [sample]
            assert cells("MS Assay Name") == [sample]
            assert cells("Raw Spectral Data File") == [sample + ".raw"]
            assert cells("Parameter Value[Instrument]") == ["Q Exactive Plus"]
            assert (
                cell_after("Parameter Value[Instrument]", 2)
                == "http://purl.obolibrary.org/obo/MS_1002857"
            )
            assert cells("Parameter Value[Number of scans]") == ["24000"]
            assert cells("Parameter Value[Scan polarity]")[0].startswith(
                "positive" if polarity == "POS" else "negative"
            )
            assert cells("Parameter Value[Mass analyzer]") == ["orbitrap"]

    def test_study_file(self, metalist, outdir):
        isa_tab = ISA_Tab(str(outdir), "MTBLS0001")
        isa_tab.write(metalist, "raw")
        with open(outdir / "s_MTBLS0001.txt") as f:
            lines = f.read().splitlines()
        assert len(lines) == 3  # header + 2 samples
        assert lines[1].split("\t")[0] == '"sample1"'
        assert lines[2].split("\t")[0] == '"sample2"'

    def test_investigation_file(self, metalist, outdir):
        isa_tab = ISA_Tab(str(outdir), "MTBLS0001")
        isa_tab.write(metalist, "raw")
        content = (outdir / "i_Investigation.txt").read_text()
        assert 'Study File Name\t"s_MTBLS0001.txt"' in content
        assert f"raw2isa {__version__}" in content
        # both assay files are referenced
        assert "a_MTBLS0001_POS_metabolite_profiling_mass_spectrometry.txt" in content
        assert "a_MTBLS0001_NEG_metabolite_profiling_mass_spectrometry.txt" in content
        # protocol parameters have a consistent number of separators
        param_lines = {}
        for line in content.splitlines():
            if line.startswith("Study Protocol Parameters Name"):
                key, *fields = line.split("\t")
                param_lines[key] = [field.strip('"') for field in fields]
                assert len(fields) == 5
        assert len(param_lines) == 3
        for key, fields in param_lines.items():
            for field in fields:
                params, separators = field.split(";"), field.count(";")
                assert len(params) == separators + 1, (key, field)
        # the separator-only lines have as many separators as parameters
        for key, fields in param_lines.items():
            if key.endswith(("Term Accession Number", "Term Source REF")):
                reference = param_lines["Study Protocol Parameters Name"]
                for field, ref in zip(fields, reference, strict=False):
                    assert field.count(";") == ref.count(";"), (key, field, ref)

    def test_no_split(self, metalist, outdir):
        isa_tab = ISA_Tab(str(outdir), "MTBLS0001")
        isa_tab.write(metalist, "raw", split=False)
        files = sorted(os.listdir(outdir))
        assert files == [
            "a_MTBLS0001_metabolite_profiling_mass_spectrometry.txt",
            "i_Investigation.txt",
            "s_MTBLS0001.txt",
        ]
        with open(outdir / "a_MTBLS0001_metabolite_profiling_mass_spectrometry.txt", newline="") as f:
            rows = list(csv.reader(f, delimiter="\t"))
        assert len(rows) == 3  # header + both samples

    def test_usermeta(self, metalist, outdir):
        usermeta = {
            "characteristics": {
                "organism": {
                    "name": "Mus musculus",
                    "ref": "NCBITAXON",
                    "accession": "http://purl.bioontology.org/ontology/NCBITAXON/10090",
                }
            }
        }
        isa_tab = ISA_Tab(str(outdir), "MTBLS0001", usermeta=usermeta)
        isa_tab.write(metalist, "raw")
        with open(outdir / "s_MTBLS0001.txt") as f:
            content = f.read()
        assert "Mus musculus" in content


# --- usermeta loader ----------------------------------------------------------------


class TestUserMetaLoader:
    def test_json_string(self):
        loader = UserMetaLoader('{"study": {"title": "My study"}}')
        assert loader.usermeta == {"study": {"title": "My study"}}

    def test_json_file(self, tmp_path):
        path = tmp_path / "meta.json"
        path.write_text(json.dumps({"study": {"title": "My study"}}))
        loader = UserMetaLoader(str(path))
        assert loader.usermeta == {"study": {"title": "My study"}}

    def test_none(self):
        assert UserMetaLoader(None).usermeta is None

    def test_example_file(self):
        here = os.path.join(os.path.dirname(__file__), "..")
        path = os.path.join(here, "examples", "usermeta.json")
        if os.path.exists(path):
            loader = UserMetaLoader(path)
            assert loader.usermeta["study"]["title"] == "Example metabolomics study"


# --- CLI ------------------------------------------------------------------------------


class TestCLI:
    def test_main(self, tmp_path, monkeypatch):
        from raw2isa import parsing

        # create a fake study folder with fake RAW files
        study = tmp_path / "study"
        study.mkdir()
        for name in ["a.raw", "b.RAW", "notaraw.txt"]:
            (study / name).touch()

        # mock the RAW parser so no .NET runtime is needed
        class FakeRawFile:
            def __init__(self, path):
                self.path = path

            @property
            def metadata(self):
                return fake_meta(os.path.basename(self.path).split(".")[0])

        monkeypatch.setattr(parsing, "ThermoRawFile", FakeRawFile)

        out = tmp_path / "out"
        parsing.main(["-i", str(study), "-o", str(out), "-s", "MTBLS0001"])
        files = sorted(os.listdir(out))
        assert "i_Investigation.txt" in files
        assert "s_MTBLS0001.txt" in files
        assert len(files) == 3  # single polarity -> single assay file

    def test_main_no_files(self, tmp_path, capsys):
        from raw2isa import parsing

        study = tmp_path / "empty"
        study.mkdir()
        with pytest.warns(UserWarning):
            parsing.main(["-i", str(study), "-o", str(tmp_path / "out"), "-s", "X"])
