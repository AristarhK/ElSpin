"""Independent conservation, time-refinement, empty-seed and input checks for air plasma."""
from pathlib import Path
import argparse
import csv
import math
import subprocess
import tempfile


def read_csv(path):
    with Path(path).open() as stream:
        return [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', default='./electrospinning_debug')
    args = parser.parse_args()
    exe = Path(args.exe).resolve()
    project = Path(__file__).resolve().parents[1]
    print(subprocess.run([str(exe), '--plasma-test'], cwd=project,
                         text=True, capture_output=True, check=True).stdout.strip())
    table = project / 'air_transport.csv'
    flow = project / 'config_spinning.nml'
    results = []
    with tempfile.TemporaryDirectory(prefix='air-plasma-tests-') as tmp:
        root = Path(tmp)

        def run(name, extra='', good=True):
            work = root/name
            work.mkdir()
            cfg = work/'plasma.nml'
            cfg.write_text('&air_plasma\n p_nr=32, p_nz=64, p_end_time=1e-10, p_dt_max=5e-13,\n'
                           f" p_transport_file='{table}',\n" + extra + '\n/\n')
            proc = subprocess.run([str(exe), '--config', str(flow), '--preset', 'quick',
                                   '--plasma-config', str(cfg)], cwd=work,
                                  text=True, capture_output=True, timeout=180)
            if not good:
                assert proc.returncode != 0, name
                assert 'FAILED' in (work/'run_status.txt').read_text()
                assert not (work/'plasma_final.csv').exists()
                return
            if proc.returncode:
                raise AssertionError(f'{name} failed:\n{proc.stdout}\n{proc.stderr}')
            assert 'PLASMA_TRANSIENT_COMPLETED' in (work/'run_status.txt').read_text()
            initial = read_csv(work/'plasma_initial.csv')
            final = read_csv(work/'plasma_final.csv')
            history = read_csv(work/'plasma_history.csv')
            assert len(initial) == len(final) == 32*64
            assert abs(history[-1]['time_s']/1e-10-1) < 1e-10
            assert abs(final[-1]['time_s']/history[-1]['time_s']-1) < 1e-12
            for r in final:
                assert all(math.isfinite(x) for x in r.values())
                assert min(r['ne_m3'],r['n_positive_m3'],r['n_negative_m3'],r['energy_eV_m3']) >= 0
                if r['gas_cell'] == 0:
                    assert r['ne_m3'] == r['n_positive_m3'] == r['n_negative_m3'] == 0
                if r['Te_defined']:
                    assert abs(r['Te_eV']-2*r['energy_eV_m3']/r['ne_m3']/3) < 1e-12
                else:
                    assert r['Te_eV'] == 0
            assert max(r['charge_balance_relative'] for r in history) < 1e-9
            assert max(r['energy_balance_relative'] for r in history) < 1e-9
            return initial, final, history

        for dt in (1e-12,5e-13,2.5e-13):
            _, final, history = run(f'dt{dt}',f'p_dt_max={dt},')
            results.append((dt,history[-1]['electrons'],max(r['Te_eV'] for r in final)))
        coarse = abs(results[0][1]-results[1][1])
        fine = abs(results[1][1]-results[2][1])
        assert fine < .75*coarse, results
        assert fine/results[-1][1] < .02, results
        print('TEST time refinement (dt_s, integral_electrons, peak_Te_eV):',results)

        initial, final, history = run('empty','p_background_ne=0, p_seed_ne=0,')
        assert max(r['ne_m3'] for r in final) == 0
        assert max(r['Te_eV'] for r in final) == 0
        print('TEST zero electrons stay zero; no fabricated density or temperature: passed')
        run('bad_radius','p_radius=3e-4,',good=False)
        run('missing_table',"p_transport_file='not_a_table.csv',",good=False)
        print('TEST geometry/transport input rejection and failed-run status: passed')
    print('ALL PLASMA INTEGRATION TESTS PASSED')


if __name__ == '__main__':
    main()
