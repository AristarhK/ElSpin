#!/usr/bin/env bash
# Rootless Debian 13 / GCC 14 bootstrap, for the managed cloud machine only.
# apt authenticates signed repository metadata and verifies package hashes.
set -euo pipefail
if command -v gfortran >/dev/null 2>&1; then exit 0; fi
if [[ -x /workspace/toolchains/gfortran/usr/bin/gfortran-14 ]]; then
    /workspace/toolchains/gfortran/usr/bin/gfortran-14 --version
    exit 0
fi
[[ -d /workspace ]] || { echo 'Install GNU Fortran using your OS package manager.' >&2; exit 1; }
[[ $(gcc -dumpversion) == 14 ]] || { echo 'This bootstrap requires Debian GCC 14; install matching gfortran with your OS package manager.' >&2; exit 1; }
compiler_root=/workspace/toolchains
mkdir -p "$compiler_root/apt/lists/partial" "$compiler_root/apt/cache/archives/partial" "$compiler_root/deb" "$compiler_root/gfortran"
cat > "$compiler_root/apt/sources.list" <<'SOURCES'
deb https://deb.debian.org/debian trixie main
SOURCES
apt_flags=(-o "Dir::Etc::sourcelist=$compiler_root/apt/sources.list" -o Dir::Etc::sourceparts=-
    -o "Dir::State::Lists=$compiler_root/apt/lists" -o "Dir::Cache=$compiler_root/apt/cache")
apt-get "${apt_flags[@]}" update
cd "$compiler_root/deb"
apt-get "${apt_flags[@]}" download gfortran-14=14.2.0-19 gfortran-14-x86-64-linux-gnu=14.2.0-19 libgfortran-14-dev=14.2.0-19
for package in ./*.deb; do dpkg-deb -x "$package" "$compiler_root/gfortran"; done
"$compiler_root/gfortran/usr/bin/gfortran-14" --version
