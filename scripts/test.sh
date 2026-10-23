#!/bin/bash
# Run the raw2isa test suite with uv.
set -e
cd "$(dirname "$0")/.."
export PATH="$HOME/.local/bin:$PATH"
uv run --group test --group pb python scripts/check_templates.py
uv run --group test --group pb pytest tests/ "$@"
