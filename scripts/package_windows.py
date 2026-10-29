"""Assemble a real Windows portable app from Ubuntu using MinGW and signed-index packages.

Usage: python scripts/package_windows.py --python-zip /path/python-3.12.10-embed-amd64.zip
  --wheels /path/windows-wheels --cc x86_64-w64-mingw32-gcc --fc x86_64-w64-mingw32-gfortran
This is a directory distribution: keep ElSpin.exe and _runtime together.
"""
from __future__ import annotations
import argparse,base64,csv,hashlib,io,json,os,shutil,subprocess
from pathlib import Path,PurePosixPath
from zipfile import ZipFile,ZIP_DEFLATED

ROOT=Path(__file__).resolve().parents[1]


def unpack_verified_wheel(wheel: Path,target: Path):
    with ZipFile(wheel) as z:
        record=[n for n in z.namelist() if n.endswith('.dist-info/RECORD')]
        if len(record)!=1:raise ValueError('Wheel RECORD missing: '+str(wheel))
        checks={row[0]:row[1] for row in csv.reader(io.StringIO(z.read(record[0]).decode('utf-8')))}
        for name in z.namelist():
            parts=PurePosixPath(name).parts
            if '..' in parts or name.startswith('/') or '\\' in name:raise ValueError('Unsafe wheel path')
            if name.endswith('/'):continue
            data=z.read(name);check=checks.get(name,'')
            if check:
                method,value=check.split('=',1)
                actual=base64.urlsafe_b64encode(hashlib.new(method,data).digest()).rstrip(b'=').decode()
                if actual!=value:raise ValueError('Wheel integrity failure: '+name)
            if parts[0].endswith('.data'):
                if parts[1] not in ('purelib','platlib'):continue
                parts=parts[2:]
            dest=target.joinpath(*parts);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)


def main():
    p=argparse.ArgumentParser();p.add_argument('--python-zip',required=True,type=Path);p.add_argument('--wheels',required=True,type=Path)
    p.add_argument('--cc',default='x86_64-w64-mingw32-gcc');p.add_argument('--fc',default='x86_64-w64-mingw32-gfortran')
    p.add_argument('--out',default=str(ROOT/'dist'/'ElSpin-Windows-x64'),type=Path);args=p.parse_args()
    blob=args.python_zip.read_bytes()
    if hashlib.md5(blob).hexdigest()!='fe8ef205f2e9c3ba44d0cf9954e1abd3':raise ValueError('CPython release checksum mismatch')
    out=args.out.resolve()
    if out.exists():raise FileExistsError('Choose a new output directory: '+str(out))
    runtime=out/'_runtime';runtime.mkdir(parents=True);libs=runtime/'Lib'/'site-packages';libs.mkdir(parents=True)
    with ZipFile(io.BytesIO(blob)) as z:
        for n in z.namelist():
            if '..' in PurePosixPath(n).parts or n.startswith('/'):raise ValueError('Unsafe Python archive')
        z.extractall(runtime)
    (runtime/'python312._pth').write_text('python312.zip\n.\nLib/site-packages\n..\nimport site\n',encoding='ascii')
    wheels=sorted(args.wheels.glob('*.whl'))
    if not wheels:raise ValueError('No Windows wheels')
    for wheel in wheels:unpack_verified_wheel(wheel,libs)
    for name in ['electrospin','configs','config','docs','src','templates']:
        shutil.copytree(ROOT/name,out/name,ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(ROOT/'app.py',out/'app.py');shutil.copy2(ROOT/'README.md',out/'README.md')
    baseline=ROOT/'results/physics_v04_baseline'
    if not baseline.is_dir():baseline=ROOT/'results/release_baseline'
    if baseline.is_dir():shutil.copytree(baseline,out/'examples/baseline')
    if (ROOT/'results/release_report').is_dir():
        shutil.copytree(ROOT/'results/release_report',out/'results/release_report')
    if (ROOT/'results/ui_v03').is_dir():
        shutil.copytree(ROOT/'results/ui_v03',out/'results/ui_v03')
    if (ROOT/'results/physics_report_v04').is_dir():
        shutil.copytree(ROOT/'results/physics_report_v04',out/'results/physics_report_v04')
    if (ROOT/'results/experiment_v05').is_dir():
        shutil.copytree(ROOT/'results/experiment_v05',out/'results/experiment_v05',ignore=shutil.ignore_patterns('__pycache__'))
    (out/'build').mkdir()
    subprocess.run([args.fc,'-O2','-std=f2008','-ffree-line-length-none','-static','-J'+str(out/'build'),str(ROOT/'src/jet3d.f90'),'-o',str(out/'build/jet3d.exe')],check=True)
    subprocess.run([args.cc,'-O2','-static','-municode','-mwindows',str(ROOT/'src/windows_launcher.c'),'-o',str(out/'ElSpin.exe')],check=True)
    (out/'CHECK.cmd').write_text('@echo off\r\ncd /d "%~dp0"\r\n_runtime\\python.exe -X utf8 app.py --self-test\r\npause\r\n',encoding='ascii')
    (out/'START_HERE.txt').write_text('ElSpin 0.5 / Windows 10-11 x64\n\nРаспакуйте весь архив в папку. Запустите ElSpin.exe.\nНе переносите .exe отдельно от _runtime и остальных файлов.\nPython, Qt и ядро включены; сервер и доступ к Интернету не нужны.\nКнопка «Эксперимент и прогноз…»: свои измерения, подбор зависимости, проверка ошибки.\nИнструкция: docs/EXPERIMENT_RU.md. Пустые CSV-шаблоны: templates.\nПроверка установки: CHECK.cmd.\nЭта сборка подготовлена в Ubuntu. Запуск на настоящей Windows в этой среде не проверялся.\nФизическая модель исследовательская: docs/PHYSICS_AUDIT_V04_RU.pdf.\n',encoding='utf-8-sig')
    manifest={'python_url':'https://www.python.org/ftp/python/3.12.10/python-3.12.10-embed-amd64.zip',
              'python_sha256':hashlib.sha256(blob).hexdigest(),'wheels':{w.name:hashlib.sha256(w.read_bytes()).hexdigest() for w in wheels},
              'native_windows_tested':False,'build_host':'Ubuntu/Linux cross compilation'}
    (out/'BUILD_PROVENANCE.json').write_text(json.dumps(manifest,indent=2)+'\n')
    checksums={str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest() for f in out.rglob('*') if f.is_file()}
    (out/'SHA256.json').write_text(json.dumps(checksums,indent=2)+'\n')
    archive=out.parent/(out.name+'.zip')
    with ZipFile(archive,'w',ZIP_DEFLATED,compresslevel=6) as z:
        for f in out.rglob('*'):
            if f.is_file():z.write(f,str(Path(out.name)/f.relative_to(out)))
    print(archive,archive.stat().st_size)

if __name__=='__main__':main()
