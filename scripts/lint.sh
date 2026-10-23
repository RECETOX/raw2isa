#!/bin/bash
# Lint the project with ruff.
set -e
cd "$(dirname "$0")/.."
export PATH="$HOME/.local/bin:$PATH"
uv run ruff check raw2isa tests "$@"
