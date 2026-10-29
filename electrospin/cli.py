"""Terminal workflow, reproducible runs, and Cartesian/Latin-hypercube studies."""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from .config import ROOT, CORE_DEFAULTS, FIELD_DEFAULTS, load_config, validate, write_config, write_namelist


def build(debug: bool = False) -> Path:
    suffix = '.exe' if sys.platform == 'win32' else ''
    packaged = ROOT / 'build' / (('jet3d_debug' if debug else 'jet3d') + suffix)
    if getattr(sys, 'frozen', False) or sys.platform == 'win32':
        if not packaged.is_file(): raise RuntimeError('Missing packaged Fortran solver: ' + str(packaged))
        return packaged
    target = 'build/jet3d_debug' if debug else 'build/jet3d'
    subprocess.run(['make', target], cwd=ROOT, check=True)
    return ROOT / target


def provenance(config: dict) -> dict:
    hashes = {}
    for file in sorted((ROOT / 'src').glob('*.f90')) + sorted((ROOT / 'electrospin').glob('*.py')):
        hashes[str(file.relative_to(ROOT))] = hashlib.sha256(file.read_bytes()).hexdigest()
    import numpy, scipy
    return {'status':'RUNNING', 'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
            'python':sys.version, 'platform':platform.platform(),
            'versions':{'numpy':numpy.__version__, 'scipy':scipy.__version__},
            'source_sha256':hashes, 'model':f"3D reduced {config['constitutive_model']} filament; axisymmetric static needle field",
            'experimental_validation':False,
            'notes':[('Prescribed total inlet current with reduced longitudinal ohmic charge redistribution.' if config['charge_transport'] else 'Prescribed convective inlet current; no longitudinal charge redistribution.'),
                     'Relative permittivity is a relaxation diagnostic; dielectric polarization is not resolved.',
                     'Electrode is a straight hollow 22G-like needle with conducting hemispherical meniscus; bevel omitted.']}


def run_case(config: dict, out: Path, *, executable: Path | None = None, report: bool = True,
             progress=None, cancelled=None) -> dict:
    config = validate(config)
    out = out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(f'Output directory is not empty: {out}. Choose a new --out to preserve earlier results.')
    out.mkdir(parents=True, exist_ok=True)
    exe = (executable or build()).resolve()
    write_config(config,out/'config.json')
    manifest = provenance(config)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    try:
        if cancelled and cancelled(): raise InterruptedError('Расчёт отменён')
        effective = dict(config)
        if config['field_mode'] == 'needle':
            from .field import solve_field
            solve_field(config,out/'field')
            effective['field_file'] = 'field/potential_rz.csv'
        elif config['field_mode'] == 'uniform': effective['field_file'] = ''
        else:
            src = Path(config['field_file']).expanduser().resolve()
            if not src.is_file(): raise ValueError(f'External field file missing: {src}')
            (out/'field').mkdir()
            shutil.copy2(src,out/'field'/'potential_rz.csv')
            effective['field_file'] = 'field/potential_rz.csv'
        write_namelist(effective,out/'jet3d.nml')
        if cancelled and cancelled(): raise InterruptedError('Расчёт отменён')
        started = time.monotonic()
        with (out/'console.log').open('w') as log:
            # Relative ASCII paths keep Fortran I/O independent of Windows ANSI
            # code pages when the selected results folder has Cyrillic names.
            proc = subprocess.Popen([str(exe),'jet3d.nml'],cwd=out,stdout=log,stderr=subprocess.STDOUT,
                                    creationflags=subprocess.CREATE_NO_WINDOW if sys.platform=='win32' else 0)
            try:
                while proc.poll() is None:
                    if cancelled and cancelled(): raise InterruptedError('Расчёт отменён; частичные данные сохранены')
                    if progress:
                        p=out/'diagnostics.csv'
                        if p.exists():
                            with p.open('rb') as f:
                                f.seek(max(0,p.stat().st_size-4096))
                                lines=f.read().decode('ascii',errors='ignore').splitlines()
                            if len(lines)>1:
                                try: progress(min(1.,float(lines[-1].split(',')[0])/config['duration']))
                                except ValueError: pass
                    time.sleep(.1)
            finally:
                if proc.poll() is None: proc.terminate();proc.wait()
        manifest['wall_time_solver_s'] = time.monotonic()-started
        manifest['executable_sha256'] = hashlib.sha256(exe.read_bytes()).hexdigest()
        status = (out/'run_status.txt').read_text() if (out/'run_status.txt').exists() else 'NO_STATUS'
        if proc.returncode or not status.startswith('COMPLETED'):
            tail = (out/'console.log').read_text(errors='replace')[-3500:]
            raise RuntimeError(f'3D solver failed (exit {proc.returncode}) in {out}:\n{tail}')
        from .analysis import summarize_run
        summary = summarize_run(out)
        from .instability import analyze_instability
        summary.update(analyze_instability(out))
        from .drying import write_drying,summarize_drying
        from .branching import run_scenario
        from .analysis import read_table
        data=read_table(out/'trajectory.csv')
        write_drying(out,data,config)
        summary.update(summarize_drying(data,config))
        from .stability import write_screening
        summary.update(write_screening(out,data,config))
        summary.update(run_scenario(out,data,config))
        (out/'analysis.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        if report:
            from .visualize import render_report
            render_report(out)
        manifest['status'] = 'COMPLETED'
        return summary
    except Exception as exc:
        manifest.update(status='CANCELLED' if isinstance(exc,InterruptedError) else 'FAILED',error=str(exc))
        (out/'run_status.txt').write_text(manifest['status']+': '+str(exc)+'\n',encoding='utf-8')
        raise
    finally:
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')


# Text input uses common laboratory units; JSON/Fortran always use SI.
PROMPTS = [
    ('voltage','Напряжение [кВ]',1e3,0), ('gap','Расстояние от мениска до коллектора [см]',.01,0),
    ('emitter_radius','Наружный диаметр иглы [мм]',.0005,0),
    ('needle_inner_radius','Внутренний диаметр иглы [мм]',.0005,0),
    ('needle_length','Длина моделируемой выступающей иглы [мм]',.001,0),
    ('field_radius','Радиус заземлённой камеры [см]',.01,0),
    ('flow_rate','Расход раствора [мл/ч]',1e-6/3600,0),
    ('jet_radius','Радиус струи после мениска [мкм]',1e-6,0),
    ('polymer_fraction','Массовая доля полимера [%]',.01,0),
    ('density_polymer','Плотность полимера [кг/м³]',1,0),
    ('density_solvent','Плотность растворителя [кг/м³]',1,0),
    ('viscosity','Динамическая вязкость раствора [Па·с]',1,0),
    ('relaxation_time','Время вязкоупругой релаксации [мс]',.001,0),
    ('surface_tension','Поверхностное натяжение [мН/м]',.001,0),
    ('conductivity','Электропроводность [мкСм/м] (диагностика зарядовой релаксации)',1e-6,0),
    ('relative_permittivity','Относительная диэлектрическая проницаемость',1,0),
    ('current','Конвективный ток струи [нА], НЕ полный ток установки',1e-9,0),
    ('temperature','Температура раствора [°C]',1,273.15),
    ('ambient_temperature','Температура воздуха [°C]',1,273.15),
    ('humidity','Относительная насыщенность воздуха выбранным растворителем [%]',.01,0),
    ('molar_mass','Молярная масса растворителя [г/моль]',.001,0),
    ('latent_heat','Удельная теплота испарения [кДж/кг]',1000,0),
    ('psat_ref','Давление насыщенного пара при опорной температуре [Па]',1,0),
    ('psat_temperature','Опорная температура давления пара [°C]',1,273.15),
    ('mass_transfer','Коэффициент массоотдачи [м/с]',1,0),
    ('heat_transfer','Коэффициент теплоотдачи [Вт/(м²·К)]',1,0),
    ('duration','Длительность расчёта [мс]',.001,0),
    ('injection_interval','Интервал дискретизации поступающего раствора [мс]',.001,0),
]


def wizard(args) -> int:
    c = load_config(args.config)
    print('Игла 22G по умолчанию: OD 0.70 мм, ID 0.40 мм. Enter сохраняет значение.')
    print('Пример раствора иллюстративный. Вводите измеренные свойства своего раствора.')
    for key,label,scale,offset in PROMPTS:
        while True:
            raw = input(f'{label} [{(c[key]-offset)/scale:.7g}]: ').strip()
            if not raw: break
            try:
                value = float(raw.replace(',','.'))*scale+offset
                if not __import__('math').isfinite(value): raise ValueError()
                c[key] = value
                break
            except ValueError: print('Нужно конечное число.')
    c['max_beads'] = max(c['max_beads'],int(c['duration']/c['injection_interval'])+4)
    c = validate(c)
    path = Path(args.save)
    if path.exists(): raise ValueError(f'Configuration already exists: {path}; use another --save')
    path.parent.mkdir(parents=True,exist_ok=True)
    write_config(c,path)
    print(f'Параметры сохранены: {path}')
    if args.run:
        summary = run_case(c,Path(args.out),report=not args.no_report)
        print(json.dumps(summary,ensure_ascii=False,indent=2))
    return 0


def sweep(args) -> int:
    config = load_config(args.config,args.set)
    out = Path(args.out).resolve()
    if out.exists() and any(out.iterdir()): raise ValueError(f'Output directory not empty: {out}')
    out.mkdir(parents=True,exist_ok=True)
    axes = {}
    for param in args.param or []:
        key,sep,text = param.partition('=')
        if not sep or key not in config or isinstance(config[key],(list,str,dict,bool)):
            raise ValueError('--param expects numeric key=value1,value2,...')
        values = [float(v) for v in text.split(',')]
        if not values: raise ValueError('Empty parameter list')
        if key in axes: raise ValueError(f'Duplicate parameter axis: {key}')
        axes[key] = values
    cases = [dict(zip(axes,values)) for values in itertools.product(*axes.values())]
    if args.lhs:
        if axes: raise ValueError('Use --param or --lhs/--range, not both')
        ranges = {}
        for text in args.ranges or []:
            key,sep,val = text.partition('=')
            if not sep or key not in config or isinstance(config[key],(list,str,dict,bool)):
                raise ValueError('--range expects numeric key=lower:upper')
            lo,hi = map(float,val.split(':'))
            if not lo < hi: raise ValueError('Each range requires lower < upper')
            ranges[key] = (lo,hi)
        if not ranges: raise ValueError('--lhs requires at least one --range')
        from scipy.stats import qmc
        points = qmc.LatinHypercube(d=len(ranges),seed=args.seed).random(args.lhs)
        points = qmc.scale(points,[v[0] for v in ranges.values()],[v[1] for v in ranges.values()])
        cases = [dict(zip(ranges,row.tolist())) for row in points]
    elif args.ranges: raise ValueError('--range requires --lhs')
    if not axes and not args.lhs: raise ValueError('Specify --param or --lhs with --range')
    if len(cases)>args.max_runs: raise ValueError(f'{len(cases)} runs exceeds --max-runs={args.max_runs}')
    # Validate every case before launching any expensive computation.
    for case in cases: validate({**config,**case})
    exe = build(args.debug)
    def one(item):
        i,changes = item
        dest = out/f'case_{i:04d}'
        row = {'case':i,**changes,'status':'FAILED','directory':str(dest)}
        try:
            result = run_case({**config,**changes},dest,executable=exe,report=args.reports)
            row.update(status='COMPLETED')
            for key,value in result.items():
                if value is None or isinstance(value,(float,int,str,bool)): row[key] = value
        except Exception as exc: row['error'] = str(exc)
        print(f'[{i+1}/{len(cases)}] {row["status"]} {changes}',flush=True)
        return row
    if args.workers < 1 or args.workers > 8: raise ValueError('--workers must be between 1 and 8')
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(one,enumerate(cases)))
    names = list(dict.fromkeys(k for row in rows for k in row))
    with (out/'study.csv').open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=names);writer.writeheader();writer.writerows(rows)
    write_config(config,out/'base_config.json')
    (out/'study_design.json').write_text(json.dumps({'seed':args.seed,'cases':cases},indent=2)+'\n')
    from .study import render_study
    render_study(out,rows,list(axes) or list(ranges))
    failed = sum(row['status']!='COMPLETED' for row in rows)
    print(f'Study: {len(rows)-failed} completed, {failed} failed. {out}/study.csv')
    return 1 if failed else 0


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description='3D electrospinning: Fortran dynamics, syringe needle, SI configuration')
    commands = p.add_subparsers(dest='command',required=True)
    run = commands.add_parser('run',help='Run a configured simulation')
    run.add_argument('--config',type=Path);run.add_argument('--out',type=Path,required=True)
    run.add_argument('--set',action='append',default=[],metavar='KEY=JSON')
    run.add_argument('--no-report',action='store_true');run.add_argument('--debug',action='store_true')
    wiz = commands.add_parser('wizard',help='Interactive Russian terminal input of solution and setup')
    wiz.add_argument('--config',type=Path);wiz.add_argument('--save',default='configs/my_solution.json')
    wiz.add_argument('--run',action='store_true');wiz.add_argument('--out',default='results/my_solution')
    wiz.add_argument('--no-report',action='store_true')
    field = commands.add_parser('field',help='Solve electrostatic field of a hollow syringe needle')
    field.add_argument('--config',type=Path);field.add_argument('--out',type=Path,required=True)
    field.add_argument('--set',action='append',default=[],metavar='KEY=JSON')
    report = commands.add_parser('report',help='Regenerate figures from existing CSV data')
    report.add_argument('directory',type=Path)
    study = commands.add_parser('sweep',help='Cartesian or Latin hypercube parameter study')
    study.add_argument('--config',type=Path);study.add_argument('--out',required=True)
    study.add_argument('--set',action='append',default=[],metavar='KEY=JSON')
    study.add_argument('--param',action='append',metavar='KEY=VALUE1,VALUE2')
    study.add_argument('--lhs',type=int,default=0,metavar='N')
    study.add_argument('--range',dest='ranges',action='append',metavar='KEY=LOW:HIGH')
    study.add_argument('--seed',type=int,default=42);study.add_argument('--workers',type=int,default=1)
    study.add_argument('--max-runs',type=int,default=100)
    study.add_argument('--reports',action='store_true');study.add_argument('--debug',action='store_true')
    commands.add_parser('build',help='Build both original and 3D Fortran solvers')
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    os.environ.setdefault('MPLCONFIGDIR','/tmp/electrospin-mpl')
    try:
        if args.command == 'build': subprocess.run(['make','all'],cwd=ROOT,check=True);return 0
        if args.command == 'wizard': return wizard(args)
        if args.command == 'sweep': return sweep(args)
        if args.command == 'report':
            from .visualize import render_report
            for file in render_report(args.directory.resolve()): print(file)
            return 0
        c = load_config(args.config,args.set)
        if args.command == 'field':
            from .field import solve_field
            if args.out.exists() and any(args.out.iterdir()): raise ValueError('Field output directory is not empty')
            print(solve_field(c,args.out.resolve()));return 0
        summary = run_case(c,args.out,executable=build(args.debug),report=not args.no_report)
        print(json.dumps(summary,ensure_ascii=False,indent=2))
        print(f'Results: {args.out.resolve()}')
        return 0
    except (ValueError,RuntimeError,OSError,subprocess.CalledProcessError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        return 1
    except (KeyboardInterrupt,EOFError):
        print('Interrupted.',file=sys.stderr);return 130
