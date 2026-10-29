#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ ! -x .venv/bin/python ]]; then python3 -m venv .venv; fi
.venv/bin/python -m pip install --disable-pip-version-check -r requirements.lock.txt
make build/jet3d
make build/jet3d_debug
printf '\nRun: .venv/bin/python app.py\nCheck without display: .venv/bin/python app.py --self-test\n'
