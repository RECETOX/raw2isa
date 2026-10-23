#!/bin/bash
# End-to-end CLI run: convert the epyraw sample.raw into an ISATab.
set -e
cd "$(dirname "$0")/.."
export PATH="$HOME/.local/bin:$PATH"
rm -rf examples/out
mkdir -p examples/study
cp -f ~/dev/epyraw/tests/testdata/sample.raw examples/study/
uv run --group pb raw2isa -i examples/study -o examples/out -s MTBLS1234 -m examples/usermeta.json -v
echo "=== files ==="
ls examples/out
echo "=== study ==="
cat examples/out/s_MTBLS1234.txt
