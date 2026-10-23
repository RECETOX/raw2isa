#!/bin/bash
# Install uv (if missing) using the official standalone installer.
set -e
if [ -x "$HOME/.local/bin/uv" ]; then
  "$HOME/.local/bin/uv" --version
  exit 0
fi
curl -LsSf https://astral.sh/uv/install.sh | sh
"$HOME/.local/bin/uv" --version
