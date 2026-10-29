"""Create a portable source + worked-example archive, without environments/caches."""
from pathlib import Path
import hashlib
import zipfile

root=Path(__file__).resolve().parents[1]
out=root/'dist/electrospinning_3d.zip'
out.parent.mkdir(exist_ok=True)
files=[]
for name in ['README.md','Makefile','.gitignore','pyproject.toml','requirements.txt','requirements.lock.txt',
             'src','config','configs','docs','electrospin','tests','scripts','validation','original']:
    p=root/name
    candidates=[p] if p.is_file() else list(p.rglob('*'))
    for f in candidates:
        rel=f.relative_to(root)
        if not f.is_file() or '__pycache__' in rel.parts or f.suffix in ('.pyc','.mod','.o'): continue
        if rel.parts[:2]==('original','validation'):continue
        if f.name in ('electrospinning','electrospinning_debug'):continue
        files.append(f)
for name in ['demo','whipping_example']:
    files.extend(p for p in (root/'results'/name).rglob('*') if p.is_file())
for name in ['sweep','lhs_smoke','convergence_full']:
    files.extend(p for p in (root/'results'/name).glob('*') if p.is_file())
files=sorted(set(files))
sums=''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(root)}\n' for p in files)
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files:z.write(p,Path('electrospinning_3d')/p.relative_to(root))
    z.writestr('electrospinning_3d/SHA256SUMS',sums)
with zipfile.ZipFile(out) as z:
    bad=z.testzip()
    if bad:raise ValueError(f'Corrupt archive entry: {bad}')
print(f'{out}: {len(files)} files, {out.stat().st_size/1e6:.2f} MB; ZIP integrity checked')
