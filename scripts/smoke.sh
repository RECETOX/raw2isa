#!/bin/bash
# Final smoke test: lockfile, CLI, and build of the sdist/wheel.
set -e
cd "$(dirname "$0")/.."
export PATH="$HOME/.local/bin:$PATH"
echo "--- uv.lock check ---"
uv lock --check
echo "--- CLI ---"
uv run raw2isa --version
echo "--- build ---"
rm -rf dist
uv build
echo "--- package data in wheel ---"
uv run python - <<'EOF'
import zipfile
import glob

wheel = glob.glob("dist/*.whl")[0]
names = zipfile.ZipFile(wheel).namelist()
data = [n for n in names if n.endswith((".obo", ".txt")) and ("templates" in n or "ontologies" in n)]
print(wheel)
for n in sorted(data):
    print(" ", n)
assert any(n.endswith("psi-ms.obo") for n in data), "ontology missing from wheel"
assert any(n.endswith("a_raw.txt") for n in data), "templates missing from wheel"
print("package data OK")
EOF
