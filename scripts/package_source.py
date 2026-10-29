"""Package current source and compact reproducible examples, excluding environments."""
from pathlib import Path
import hashlib
import json
import tomllib
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[1]


def main():
    version = tomllib.loads((ROOT/'pyproject.toml').read_text())['project']['version']
    release = '.'.join(version.split('.')[:2])
    out = ROOT/f'dist/ElSpin-source-v{release}.zip'
    if out.exists():
        raise FileExistsError(out)
    files = set()
    for name in ['README.md', 'Makefile', '.gitignore', 'app.py', 'pyproject.toml',
                 'requirements.txt', 'requirements.lock.txt', 'requirements-desktop.txt',
                 'electrospin', 'src', 'config', 'configs', 'docs', 'templates', 'tests', 'scripts', 'original']:
        path = ROOT/name
        for f in ([path] if path.is_file() else path.rglob('*')):
            rel = f.relative_to(ROOT)
            if not f.is_file() or '__pycache__' in rel.parts or f.suffix in ('.pyc', '.mod', '.o', '.aux', '.log', '.out', '.toc'):
                continue
            if rel.parts[:2] == ('original', 'validation') or f.name in ('electrospinning', 'electrospinning_debug'):
                continue
            files.add(f)
    for name in ['physics_v04_baseline', 'physics_report_v04', 'experiment_v05']:
        files.update(f for f in (ROOT/'results'/name).rglob('*') if f.is_file())
    for name in ['physics_audit_v04', 'physics_material_audit_v04']:
        files.update(f for f in (ROOT/'results'/name).rglob('*') if f.is_file() and
                     (f.parent == ROOT/'results'/name or f.name in
                      ('config.json', 'analysis.json', 'manifest.json', 'drying_history.csv', 'branching_screen.csv')))
    for pattern in ['tests-*v05.log', 'packaging-*v05.log']:
        files.update((ROOT/'results').glob(pattern))
    out.parent.mkdir(exist_ok=True)
    checksums = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)}
    with ZipFile(out, 'x', ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(files):
            z.write(p, str(Path('ElSpin')/p.relative_to(ROOT)))
        z.writestr('ElSpin/SHA256.json', json.dumps(checksums, indent=2)+'\n')
    with ZipFile(out) as z:
        bad = z.testzip()
        if bad:
            raise ValueError('Corrupt ZIP entry: '+bad)
    print(out, len(files), out.stat().st_size)


if __name__ == '__main__':
    main()
