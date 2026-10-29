"""Refine the plasma patch while keeping liquid field and physical initial data fixed."""
from pathlib import Path
import argparse,csv,subprocess,time


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--exe',default='./electrospinning')
    parser.add_argument('--out',default='validation/plasma/mesh')
    args=parser.parse_args()
    project=Path(__file__).resolve().parent
    exe=Path(args.exe).resolve()
    out=Path(args.out).resolve(); out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for nr,nz in [(32,64),(64,128),(128,256)]:
        work=out/f'{nr}x{nz}'; work.mkdir(exist_ok=True)
        cfg=work/'plasma.nml'
        cfg.write_text('&air_plasma\n'+f" p_nr={nr}, p_nz={nz}, p_dt_max=2.5e-13,\n"
                       +f" p_transport_file='{project/'air_transport.csv'}',\n/\n")
        start=time.perf_counter()
        proc=subprocess.run([str(exe),'--config',str(project/'config_spinning.nml'),'--preset','fine',
                             '--plasma-config',str(cfg)],cwd=work,text=True,capture_output=True)
        (work/'console.log').write_text(proc.stdout+proc.stderr)
        if proc.returncode:
            raise RuntimeError(f'{nr}x{nz}: inspect {work}/console.log')
        with (work/'plasma_history.csv').open() as stream:
            h=list(csv.DictReader(stream))[-1]
        row={'Nr':nr,'Nz':nz,'dt_max_s':2.5e-13,'end_time_s':float(h['time_s']),
             'electrons':float(h['electrons']),'peak_ne_m3':float(h['ne_peak_m3']),
             'peak_Te_eV':float(h['Te_peak_eV']),
             'charge_error':float(h['charge_balance_relative']),
             'energy_error':float(h['energy_balance_relative']),
             'seconds':time.perf_counter()-start}
        rows.append(row)
        print(row,flush=True)
    with (out/'plasma_mesh_convergence.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)


if __name__=='__main__':
    main()
