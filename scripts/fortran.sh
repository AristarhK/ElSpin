#!/usr/bin/env bash
set -euo pipefail
if command -v gfortran >/dev/null 2>&1; then
    exec gfortran "$@"
elif [[ -x /workspace/toolchains/gfortran/usr/bin/gfortran-14 ]]; then
    exec /workspace/toolchains/gfortran/usr/bin/gfortran-14 -B/usr/lib/gcc/x86_64-linux-gnu/14/ "$@"
else
    echo 'GNU Fortran missing. Install gfortran using your OS package manager, then retry.' >&2
    exit 127
fi
