"""Integration limits and conservation checks of the actual Fortran solver.

Run: python -m unittest discover -s tests -p test_core.py -v
Uses debug executable when available; make build/jet3d_debug first.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / 'build' / 'jet3d_debug'
if not EXE.exists():
    EXE = ROOT / 'build' / 'jet3d'


def nml_value(value):
    if isinstance(value, bool):
        return '.true.' if value else '.false.'
    if isinstance(value, list):
        return ', '.join(nml_value(v) for v in value)
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    return str(value)


def rows(path):
    with path.open() as stream:
        data = list(csv.DictReader(stream))
    for row in data:
        if None in row or any(v is None or v.strip() == '' for v in row.values()):
            raise AssertionError(f'Malformed CSV row in {path}: {row}')
        for value in row.values():
            if not math.isfinite(float(value)):
                raise AssertionError(f'Nonfinite CSV field in {path}')
    return [{k: float(v) for k, v in r.items()} for r in data]


@unittest.skipUnless(EXE.exists(), 'build Fortran executable before running tests')
class CoreIntegration(unittest.TestCase):
    def run_case(self, overrides, expect_ok=True):
        work = Path(tempfile.mkdtemp(prefix='electrospin-core-'))
        self.addCleanup(lambda: __import__('shutil').rmtree(work, ignore_errors=True))
        config = json.loads((ROOT / 'config' / 'core_defaults.json').read_text())
        config.update(overrides)
        cfg = work / 'case.nml'
        cfg.write_text('&jet3d\n' + '\n'.join(f'{k}={nml_value(v)}' for k, v in config.items()) + '\n/\n')
        proc = subprocess.run([str(EXE), str(cfg)], cwd=work, capture_output=True, text=True, timeout=90)
        if not expect_ok:
            self.assertNotEqual(proc.returncode, 0, proc.stdout)
            self.assertTrue((work / 'run_status.txt').read_text().startswith('FAILED'))
            return work
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertEqual((work / 'run_status.txt').read_text().strip(), 'COMPLETED')
        trajectory = rows(work / 'trajectory.csv')
        diagnostic = rows(work / 'diagnostics.csv')
        deposited = rows(work / 'deposition.csv')
        for field in ('mass_relative_error', 'polymer_relative_error', 'charge_relative_error'):
            self.assertLess(max(r[field] for r in diagnostic), 1e-10, field)
        return trajectory, diagnostic, deposited

    @staticmethod
    def neutral():
        return dict(voltage=0., current=0., surface_tension=0., gravity=[0., 0., 0.],
                    drag_coefficient=0., heat_transfer=0., evaporation=False, coulomb=False,
                    image_force=False, perturbation_amplitude=0., duration=.001,
                    injection_interval=.0001, output_interval=.0001, dt_max=1e-6)

    def test_actual_forcefree_inlet_has_no_artificial_strain(self):
        cfg = self.neutral()
        trajectory, diagnostic, deposited = self.run_case(cfg)
        final = [r for r in trajectory if abs(r['time_s'] - cfg['duration']) < 1e-12]
        self.assertGreater(len(final), 5)
        u0 = 1.388888888888889e-10 / (math.pi * (3e-6)**2)
        for r in final:
            self.assertAlmostEqual(r['vz_m_s'] / u0, 1., delta=1e-11)
            self.assertAlmostEqual(r['radius_m'] / 3e-6, 1., delta=1e-10)
            self.assertLess(abs(r['stress_Pa']), 1e-6)
            self.assertEqual(r['x_m'], 0.)
            self.assertEqual(r['y_m'], 0.)
        self.assertFalse(deposited)
        self.assertAlmostEqual(diagnostic[-1]['injected_mass_kg'], diagnostic[-1]['retained_mass_kg'], delta=1e-22)

    def test_evaporation_cools_and_preserves_polymer(self):
        cfg = self.neutral()
        cfg.update(evaporation=True, mass_transfer=2., duration=.002)
        trajectory, diagnostic, _ = self.run_case(cfg)
        self.assertGreater(diagnostic[-1]['evaporated_mass_kg'], 0.)
        last = [r for r in trajectory if abs(r['time_s'] - cfg['duration']) < 1e-12]
        first = min(last, key=lambda r: r['bead_id'])
        self.assertLess(first['temperature_K'], 298.15)
        self.assertLess(first['radius_m'], 3e-6)
        self.assertGreater(first['polymer_mass_kg'], 0.)

    def test_polarity_symmetry_and_coulomb_dynamics(self):
        cfg = dict(duration=.004, injection_interval=.0002, output_interval=.001, dt_max=3e-6,
                   evaporation=False, image_force=False)
        positive, _, _ = self.run_case(cfg)
        negative, _, _ = self.run_case(dict(cfg, voltage=-10000., current=-2e-8))
        self.assertEqual(len(positive), len(negative))
        for a, b in zip(positive, negative):
            for key in ('x_m', 'y_m', 'z_m', 'vx_m_s', 'vy_m_s', 'vz_m_s', 'radius_m'):
                self.assertAlmostEqual(a[key], b[key], delta=1e-12 * max(1., abs(a[key])))
            self.assertAlmostEqual(a['charge_C'], -b['charge_C'], delta=1e-24)
        no_pairs, _, _ = self.run_case(dict(cfg, coulomb=False))
        # Actual nonzero Coulomb force changes the integrated state, beyond seed geometry.
        pa = max(r['z_m'] for r in positive)
        pn = max(r['z_m'] for r in no_pairs)
        self.assertGreater(abs(pa - pn), 1e-5)

    def test_ballistic_crossings_and_moving_collector_frame(self):
        cfg = self.neutral()
        cfg.update(gap=.005, duration=.002, viscosity=1e-12,
                   collector_velocity=[.01, -.02], output_interval=.0002)
        trajectory, diagnostic, deposited = self.run_case(cfg)
        self.assertGreater(len(deposited), 2)
        u0 = 1.388888888888889e-10 / (math.pi * (3e-6)**2)
        first = min(deposited, key=lambda r: r['time_s'])
        self.assertAlmostEqual(first['time_s'], .005/u0, delta=1e-10)
        for r in deposited:
            self.assertAlmostEqual(r['collector_x_m'], r['x_m'] - .01*r['time_s'], delta=1e-14)
            self.assertAlmostEqual(r['collector_y_m'], r['y_m'] + .02*r['time_s'], delta=1e-14)
            self.assertGreater(r['dry_radius_m'], 0.)
            self.assertLess(r['dry_radius_m'], r['radius_m'])
        self.assertGreater(diagnostic[-1]['deposited'], 0)

    def test_images_attract_actual_parcel(self):
        cfg = self.neutral()
        cfg.update(current=1e-8, image_force=False, duration=.00015, gap=.005,
                   output_interval=.00015)
        off, _, _ = self.run_case(cfg)
        on, _, _ = self.run_case(dict(cfg, image_force=True))
        self.assertGreater(on[-1]['z_m'], off[-1]['z_m'] + 1e-12)

    def test_small_output_interval_before_first_injection(self):
        cfg = self.neutral()
        cfg.update(output_interval=.00002, duration=.00035)
        _, diagnostic, _ = self.run_case(cfg)
        self.assertGreater(len(diagnostic), 10)
        times = [r['time_s'] for r in diagnostic]
        self.assertEqual(times, sorted(set(times)))

    def test_timestep_refinement_of_actual_charged_trajectory(self):
        positions = []
        for dt in (4e-6, 2e-6, 1e-6):
            trajectory, _, _ = self.run_case(dict(duration=.002, output_interval=.002,
                                                 dt_max=dt, evaporation=False))
            first = [r for r in trajectory if abs(r['time_s']-.002) < 1e-12
                     and r['bead_id'] == 1][0]
            positions.append(first['z_m'])
        coarse = abs(positions[0] - positions[1])
        fine = abs(positions[1] - positions[2])
        self.assertGreater(coarse, 1e-9)
        self.assertLess(fine, .7*coarse)
        self.assertLess(fine/positions[-1], 1e-3)

    def test_invalid_parameters_fail(self):
        self.run_case(dict(self.neutral(), humidity=1.1), expect_ok=False)
        self.run_case(dict(self.neutral(), viscosity=0.), expect_ok=False)
        self.run_case(dict(self.neutral(), max_beads=2), expect_ok=False)


if __name__ == '__main__':
    unittest.main()
