"""Reproducible model controls and a 3x3 voltage/flow study, max two solvers."""
from __future__ import annotations
import argparse,csv,itertools,json,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from electrospin.cli import build,run_case
from electrospin.materials import pva_water


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',default='results/research');parser.add_argument('--workers',type=int,default=2);args=parser.parse_args()
    root=Path(args.out).resolve();root.mkdir(parents=True,exist_ok=True);exe=build();base=pva_water()
    jobs=[('baseline',base,'baseline')]
    for i,(voltage,flow) in enumerate(itertools.product([8000.,10000.,12000.],[.3,.5,.7])):
        jobs.append((f'grid_{i:02d}',{**base,'voltage':voltage,'flow_rate':flow*1e-6/3600},'grid'))
    for seed,on in itertools.product([7,42,123],[True,False]):
        jobs.append((f'seed_{seed}_coulomb_{int(on)}_window12',{**base,'random_seed':seed,'coulomb':on,'duration':.012,'output_interval':.0002},'control'))
    # Compare identical noise path in time; the injector samples at fixed events.
    for dt in [4e-6,2e-6,1e-6]:
        jobs.append((f'time_{dt:g}',{**base,'duration':.012,'dt_max':dt,'output_interval':.0002},'time'))
    for di in [.0002,.0001,.00005]:
        jobs.append((f'space_{di:g}',{**base,'duration':.012,'injection_interval':di,'output_interval':.0002,'max_beads':500},'space'))
    for soft in [1.5e-6,3e-6,6e-6]:
        jobs.append((f'soften_{soft:g}',{**base,'duration':.012,'coulomb_softening':soft,'output_interval':.0002},'softening'))
    def one(job):
        name,c,group=job;out=root/name
        row={'name':name,'group':group,'voltage':c['voltage'],'flow_rate':c['flow_rate'],
             'random_seed':c['random_seed'],'coulomb':c['coulomb'],'dt_max':c['dt_max'],
             'injection_interval':c['injection_interval'],'coulomb_softening':c['coulomb_softening'],'status':'FAILED','directory':str(out)}
        try:
            if out.exists():
                old=json.loads((out/'config.json').read_text(encoding='utf-8'))
                if old!=c:raise ValueError('Resume config mismatch: '+name)
                if json.loads((out/'manifest.json').read_text())['status']!='COMPLETED':raise ValueError('Prior incomplete run retained: '+name)
                result=json.loads((out/'analysis.json').read_text(encoding='utf-8'))
            else:result=run_case(c,out,executable=exe,report=False)
            row.update({k:v for k,v in result.items() if v is None or isinstance(v,(int,float,str,bool))});row['status']='COMPLETED'
        except Exception as e:row['error']=str(e)
        print(name,row['status'],row.get('deposited_parcels',''),flush=True);return row
    rows=[]
    with ThreadPoolExecutor(max_workers=max(1,min(args.workers,2))) as pool:
        for row in pool.map(one,jobs):
            rows.append(row);(root/'research.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (root/'research.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    return int(any(r['status']!='COMPLETED' for r in rows))

if __name__=='__main__':raise SystemExit(main())
