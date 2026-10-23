"""Online tests of raw2isa: parsing real Thermo RAW files with epyraw.

These tests require the ``epyraw`` package (and a .NET or Mono runtime,
plus the RawFileReader assemblies, which epyraw downloads on first use)
and a sample RAW file. Point the ``RAW2ISA_TEST_RAW`` environment variable
to a RAW file to use a custom one; otherwise the test data bundled with a
local checkout of epyraw (``~/dev/epyraw/tests/testdata/sample.raw``) is
used if present.
"""

import os

import pytest

pytest.importorskip("epyraw")

from raw2isa.parsing import convert  # noqa: E402
from raw2isa.raw import ThermoRawFile  # noqa: E402


def find_sample_raw():
    candidates = []
    env = os.environ.get("RAW2ISA_TEST_RAW")
    if env:
        candidates.append(env)
    home = os.path.expanduser("~")
    candidates.append(
        os.path.join(home, "dev", "epyraw", "tests", "testdata", "sample.raw")
    )
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    return None


SAMPLE_RAW = find_sample_raw()

pytestmark = pytest.mark.skipif(
    SAMPLE_RAW is None, reason="no Thermo RAW test file available"
)


@pytest.fixture(scope="module")
def metadata():
    return ThermoRawFile(SAMPLE_RAW).metadata


class TestMetadata:
    def test_basic_keys(self, metadata):
        for key in [
            "Sample Name",
            "MS Assay Name",
            "Raw Spectral Data File",
            "Number of scans",
            "Data file content",
        ]:
            assert key in metadata, key

    def test_instrument(self, metadata):
        instrument = metadata["Instrument"]
        assert instrument["name"]
        assert instrument.get("accession", "").startswith(
            "http://purl.obolibrary.org/obo/MS_"
        )
        manufacturer = metadata["Instrument manufacturer"]
        assert "Thermo" in manufacturer["name"]

    def test_files(self, metadata):
        entry = metadata["Raw Spectral Data File"]["entry_list"][0]
        assert entry["value"].lower().endswith(".raw")
        assert metadata["Raw data file format"]["name"] == "Thermo RAW format"

    def test_urlized_accessions(self, metadata):
        for term in [metadata["Raw data file format"], metadata["Instrument"]]:
            if "accession" in term:
                assert term["accession"].startswith("http"), term

    def test_data_file_content(self, metadata):
        entries = metadata["Data file content"]["entry_list"]
        assert 1 <= len(entries) <= 2
        assert all(entry["name"].startswith("MS") for entry in entries)

    def test_polarity(self, metadata):
        if "Scan polarity" in metadata:
            assert metadata["Scan polarity"]["name"] in (
                "positive scan",
                "negative scan",
                "alternating scan",
            )


class TestConvert:
    def test_convert(self, tmp_path):
        study = tmp_path / "study"
        study.mkdir()
        # copy the sample into the study folder (the .NET memory-mapped
        # loader used by RawFileReader does not follow symlinks)
        import shutil

        link = study / os.path.basename(SAMPLE_RAW)
        shutil.copy(SAMPLE_RAW, link)
        out = tmp_path / "out"
        convert(str(study), str(out), "MTBLS0001", verbose=False)
        files = sorted(os.listdir(out))
        assert "i_Investigation.txt" in files
        assert "s_MTBLS0001.txt" in files
        assert any(name.startswith("a_MTBLS0001") for name in files)
