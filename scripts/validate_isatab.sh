#!/bin/bash
# Validate the example output with the ISA-API (isatools).
set -e
cd "$(dirname "$0")/.."
export PATH="$HOME/.local/bin:$PATH"
uv run --group test python - <<'EOF'
from isatools import isatab

with open("examples/out/i_Investigation.txt") as fp:
    report = isatab.validate(fp)
errors = report.get("errors", [])
warnings = report.get("warnings", [])
print("validation errors:", len(errors))
for e in errors[:10]:
    print("  -", e.get("message"))
print("validation warnings:", len(warnings))
for w in warnings[:10]:
    print("  -", w.get("message"))
EOF
