#!/usr/bin/env bash
set -euo pipefail
# Ubuntu 24.04 / WSL2. Run explicitly; no Windows driver changes.
repo="${1:-$HOME/vesuvius/villa}"
pin=2dcfaf6a08c3bc796fde726c4fa32050c8fc90e7
if [[ -e "$repo" ]]; then
  echo "Destination already exists: $repo. Use a new path; nothing overwritten." >&2
  exit 1
fi
sudo apt-get update
sudo apt-get install -y git build-essential g++-13 cmake ninja-build python3-venv rclone
python3 -m venv "$HOME/vesuvius-tools"
"$HOME/vesuvius-tools/bin/python" -m pip install uv
mkdir -p "$(dirname "$repo")"
git clone https://github.com/ScrollPrize/villa.git "$repo"
git -C "$repo" checkout --detach "$pin"
cd "$repo/spiral-fitting"
CC=gcc-13 CXX=g++-13 "$HOME/vesuvius-tools/bin/uv" sync --frozen
.venv/bin/python -c 'import torch; assert torch.cuda.is_available(); print(torch.cuda.get_device_name()); x=torch.empty(1, device="cuda"); torch.cuda.synchronize()'
echo 'SETUP DONE'
