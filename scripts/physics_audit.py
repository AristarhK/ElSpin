"""Matched controls, bounded hypotheses, independent drying benchmark and refinement."""
from __future__ import annotations
import argparse,csv,itertools,json,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from electrospin.config import validate
from electrospin.cli import build,run_case
from electrospin.materials import pva_water
from electrospin.drying_reference import integrate_filament


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='results/physics_audit_v04');parser.add_argument('--material-only',action='store_true');args=parser.parse_args()
    root=Path(args.out).resolve();root.mkdir(parents=True,exist_ok=True);base=pva_water();exe=build()
    jobs=[]
    for r,d in itertools.product([.0002,.001,.005],[1e-10,1e-9,1e-8]):
        jobs.append((f'film_R{r:g}_D{d:g}',{**base,'gas_boundary_radius':r,'liquid_diffusivity':d},'uncertainty'))
    jobs.extend([
        ('constant_legacy_transfer',{**base,'transport_model':'constant','activity_model':'volume_fraction','liquid_diffusion':False},'transfer'),
        ('film_no_liquid_resistance',{**base,'liquid_diffusion':False},'transfer'),
        ('film_athermal_activity',{**base,'flory_huggins_chi':0.},'activity'),
        ('crossflow_domain_control',{**base,'transport_model':'crossflow'},'transfer')])
    for dt in [4e-6,2e-6,1e-6]:jobs.append((f'time_{dt:g}',{**base,'duration':.012,'dt_max':dt,'output_interval':.0002},'time'))
    for di in [.0002,.0001,.00005]:jobs.append((f'space_{di:g}',{**base,'duration':.012,'injection_interval':di,'output_interval':.0002},'space'))
    material_jobs=[]
    for key,values in [('viscosity',[.05,.8]),('relaxation_time',[.001,.1]),('conductivity',[2e-5,2e-4])]:
        for value in values:material_jobs.append((f'{key}_{value:g}',{**base,key:value,'duration':.012,'output_interval':.0002},'material'))
    for radius in [.0002,.001,.005]:
        material_jobs.append((f'slow_diffusion_R{radius:g}',{**base,'gas_boundary_radius':radius,'liquid_diffusivity':1e-12},'uncertainty'))
    jobs=material_jobs if args.material_only else jobs+material_jobs
    def one(job):
        name,c,group=job;c=validate(c);out=root/name
        row={'name':name,'group':group,'transport_model':c['transport_model'],'gas_boundary_radius':c['gas_boundary_radius'],
             'liquid_diffusivity':c['liquid_diffusivity'],'dt_max':c['dt_max'],'injection_interval':c['injection_interval'],'status':'FAILED'}
        try:
            if (out/'analysis.json').exists():
                if json.loads((out/'config.json').read_text(encoding='utf-8'))!=c:raise ValueError('Config mismatch on resume')
                s=json.loads((out/'analysis.json').read_text(encoding='utf-8'))
            else:s=run_case(c,out,executable=exe,report=False)
            row.update({k:v for k,v in s.items() if v is None or isinstance(v,(str,float,int,bool))});row['status']='COMPLETED'
        except Exception as exc:row['error']=str(exc)
        print(name,row['status'],row.get('drying_remaining_solvent_fraction_at_first_impact'),flush=True)
        return row
    rows=[]
    with ThreadPoolExecutor(max_workers=2) as pool:
        for row in pool.map(one,jobs):
            rows.append(row);(root/'audit.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (root/'audit.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    reference=[]
    for radius in ([] if args.material_only else [3e-6,1e-6,.2e-6]):
        for rtol in [1e-5,1e-7,1e-9]:
            a=integrate_filament(base,radius,5.,np.linspace(0,.1,201),rtol=rtol,atol=rtol*.01)
            serial={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in a.items()}
            (root/f'reference_a{radius:g}_tol{rtol:g}.json').write_text(json.dumps(serial,indent=2)+'\n')
            reference.append({'radius0_m':radius,'rtol':rtol,'remaining_solvent_at_100mm':float(a['remaining_solvent_fraction'][-1]),'inlet_flow_ml_h':a['inlet_flow_ml_h']})
    (root/'reference_summary.json').write_text(json.dumps(reference,indent=2)+'\n')
    return int(any(r['status']!='COMPLETED' for r in rows))

if __name__=='__main__':raise SystemExit(main())
