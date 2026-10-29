#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
make build/jet3d
extra_data=()
if [[ -d results/physics_v04_baseline ]]; then
  extra_data+=(--add-data results/physics_v04_baseline:examples/baseline)
fi
.venv/bin/python -m PyInstaller --noconfirm --clean --windowed --name ElSpin-Linux \
  --add-binary build/jet3d:build --add-data configs:configs --add-data config:config \
  --add-data src:src --add-data docs:docs --add-data templates:templates \
  --collect-all numpy \
  --hidden-import scipy._cyutility \
  --exclude-module tkinter --exclude-module IPython --exclude-module pandas \
  "${extra_data[@]}" app.py
MPLCONFIGDIR=/tmp/elspin-mpl dist/ElSpin-Linux/ElSpin-Linux --self-test
