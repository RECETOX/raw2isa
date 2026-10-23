#!/bin/bash
# Install the .NET 8 runtime into ~/.dotnet (needed by epyraw RawFileReader).
set -e
if [ -x "$HOME/.dotnet/dotnet" ]; then
  echo ".NET already installed:"
  "$HOME/.dotnet/dotnet" --list-runtimes
  exit 0
fi
cd /tmp
curl -sSL https://dot.net/v1/dotnet-install.sh -o dotnet-install.sh
chmod +x dotnet-install.sh
./dotnet-install.sh --channel 8.0 --runtime dotnet --install-dir "$HOME/.dotnet"
"$HOME/.dotnet/dotnet" --list-runtimes
