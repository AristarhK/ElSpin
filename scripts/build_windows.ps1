$ErrorActionPreference = 'Stop'
Set-Location (Join-Path $PSScriptRoot '..')
# Use native Windows CPython 3.12 x64 and a MinGW-w64 gfortran on PATH.
New-Item -ItemType Directory -Force build | Out-Null
gfortran -O2 -std=f2008 -ffree-line-length-none -static -Jbuild src/jet3d.f90 -o build/jet3d.exe
if ($LASTEXITCODE -ne 0) { throw 'Fortran compilation failed' }
py -3.12 -m venv .venv-win
& .venv-win/Scripts/python.exe -m pip install -r requirements-desktop.txt 'pyinstaller==6.12.0'
if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed' }
& .venv-win/Scripts/python.exe -m PyInstaller --noconfirm --clean --windowed --name ElSpin `
  --add-binary 'build/jet3d.exe;build' --add-data 'configs;configs' --add-data 'config;config' `
  --add-data 'src;src' --add-data 'docs;docs' --add-data 'templates;templates' --collect-all numpy --hidden-import scipy._cyutility app.py
if ($LASTEXITCODE -ne 0) { throw 'Packaging failed' }
& dist/ElSpin/ElSpin.exe --self-test
if ($LASTEXITCODE -ne 0) { throw 'Packaged self-test failed' }
