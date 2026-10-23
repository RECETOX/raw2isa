# raw2isa

`raw2isa` is a parser/converter that translates Thermo Fisher Scientific
`.raw` mass spectrometry files into the [ISA-Tab](http://isa-tools.org/format.html)
format (Investigation–Study–Assay with TAB separated values).

It is the Thermo RAW counterpart of
[mzml2isa](https://github.com/ISA-tools/mzml2isa): given a folder of `.raw`
files, it extracts the acquisition metadata of every file and pre-populates
the ISA investigation, study and assay files with it, so that only the
metadata that cannot be found in the RAW files (organism, contacts,
publications, ...) has to be filled in manually.

Metadata extraction is powered by [epyraw](https://pypi.org/project/epyraw/),
a Python wrapper around Thermo Fisher's
[RawFileReader](https://planetorbitrap.com/rawfilereader) .NET library, which
requires a .NET 8 runtime to be available on the system.

## Features

For every `.raw` file in the input folder, `raw2isa` extracts:

- **Sample metadata**: sample name, sample type (blank/standard/QC),
  injection volume, comments;
- **Instrument metadata**: instrument model (resolved to a PSI-MS
  controlled vocabulary term, along with its manufacturer), serial
  number, acquisition software version;
- **Acquisition metadata**: scan polarity, scan m/z range, time range,
  number of scans, mass analyzer, ion source, spectrum representation
  (centroid/profile), data file content (MS1/MSn), acquisition date and
  operator — mostly parsed from the Thermo scan filter strings and
  resolved to PSI-MS terms.

and writes:

- `i_Investigation.txt` — the ISA investigation file;
- `s_<study_id>.txt` — the ISA study file, one row per sample;
- `a_<study_id>[_POL]_metabolite_profiling_mass_spectrometry.txt` — the ISA
  assay file(s), split by scan polarity by default.

The generated ISA-Tab validates with zero errors with the
[ISA-API](https://github.com/ISA-tools/ISA-api) (`isatools.isatab.validate`);
remaining warnings concern metadata that cannot be extracted from RAW files
and has to be filled in manually (see below).

## Installation

`raw2isa` and its dependency `epyraw` are available on PyPI:

```bash
pip install raw2isa
```

The underlying Thermo Fisher RawFileReader assemblies target **.NET 8**, so
a .NET 8 runtime is required (recent Mono versions may also work):

```bash
# Ubuntu/Debian: install the .NET 8 runtime (or use your package manager)
curl -sSL https://dot.net/v1/dotnet-install.sh | bash -s -- --channel 8.0 --runtime dotnet --install-dir ~/.dotnet

epyraw download  # fetch the RawFileReader .NET assemblies
```

## Usage

Command line:

```bash
raw2isa -i IN_PATH -o OUT_PATH -s STUDY_ID [options]
```

| Option | Description |
| --- | --- |
| `-i` | input folder containing Thermo RAW files (required) |
| `-o` | output folder in which the ISA-Tab files will be created (required) |
| `-s` | study identifier, e.g. `MTBLSxxx` (required) |
| `-m` | additional user-provided metadata (JSON file, XLSX file or JSON string) |
| `-j` | number of parallel jobs used for parsing the RAW files |
| `-n` | do NOT split assay files based on scan polarity |
| `-t` | directory containing custom template files |
| `-v` | show more output |

Programmatic use:

```python
from raw2isa.parsing import convert

convert("raw_files/", "out_folder/", "MTBLS1234")
```

## Development with uv

The project is managed with [uv](https://docs.astral.sh/uv/). Clone the
repository and set up the environment with:

```bash
uv sync                             # create .venv, install raw2isa (editable) + dev tools
uv sync --group test --group pb    # also install isatools and the tqdm progress bar
```

Common tasks:

```bash
uv run raw2isa --help               # run the CLI
uv run pytest                       # run the test suite
uv run --group test --group pb pytest  # ... including the isatools validation tests
uv run ruff check raw2isa tests     # lint
uv lock                             # update the lockfile
```

The test suite needs a Thermo RAW file to run the online tests; it uses the
`sample.raw` bundled with a local checkout of epyraw
(`~/dev/epyraw/tests/testdata/sample.raw`), or the file pointed to by the
`RAW2ISA_TEST_RAW` environment variable.

## User-provided metadata

RAW files only contain acquisition-related metadata; the biological and
administrative metadata (organism, contacts, publications, protocol
descriptions, chromatography instrument, ...) can be provided in a separate
file with the `-m` flag, either as JSON:

```json
{
    "characteristics": {
        "organism": {
            "name": "Mus musculus",
            "ref": "NCBITAXON",
            "accession": "http://purl.bioontology.org/ontology/NCBITAXON/10090"
        }
    },
    "study": {
        "title": "My study",
        "description": "..."
    },
    "study_contacts": [
        {
            "last_name": "Doe",
            "first_name": "Jane",
            "email": "jane.doe@example.org",
            "affiliation": "University of Example",
            "roles": {
                "name": "principal investigator",
                "ref": "PRO",
                "accession": "http://purl.obolibrary.org/obo/PRO_000000001"
            }
        }
    ]
}
```

or as an XLSX file; a template can be generated with
`python -m raw2isa.usermeta` (writes `usermeta.xlsx` to the current
directory).

## License

`raw2isa` is distributed under the GNU General Public License version 3.0
(GPLv3); see the `COPYING` file. It is derived from
[mzml2isa](https://github.com/ISA-tools/mzml2isa), by Tom Lawson and
Martin Larralde, also GPLv3.
