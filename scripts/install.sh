#!/bin/bash
# Set up the uv-managed development environment.
set -e
cd "$(dirname "$0")/.."
export PATH="$HOME/.local/bin:$PATH"
if ! command -v uv > /dev/null; then
  bash scripts/install_uv.sh
fi
# start from a clean, uv-managed virtual environment
rm -rf .venv
uv sync
uv sync --group test --group pb
uv pip list
