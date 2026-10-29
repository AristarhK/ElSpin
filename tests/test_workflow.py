"""Cross-layer test: real needle field -> Fortran integration -> numeric exports -> offline report."""
from pathlib import Path
import json
import tempfile
import unittest
from electrospin.cli import run_case
from electrospin.config import DEFAULTS, ROOT


@unittest.skipUnless((ROOT/'build/jet3d').exists(),'make all first')
class WorkflowTest(unittest.TestCase):
    def test_actual_needle_workflow_and_offline_report(self):
        with tempfile.TemporaryDirectory() as work:
            out=Path(work)/'case'
            config={**DEFAULTS,'duration':.006,'output_interval':.001,'field_nr':48,'field_nz':80}
            summary=run_case(config,out,executable=ROOT/'build/jet3d',report=True)
            self.assertGreater(summary['deposited_parcels'],0)
            self.assertLess(abs(summary['coverage_mass_integral_error_kg']),1e-20)
            self.assertLess(summary['maximum_mass_relative_error'],1e-10)
            self.assertEqual(json.loads((out/'manifest.json').read_text())['status'],'COMPLETED')
            field=json.loads((out/'field/field_summary.json').read_text())
            self.assertEqual(field['geometry']['outer_radius_m'],.00035)
            self.assertEqual(field['geometry']['inner_radius_m'],.00020)
            self.assertLess(field['relative_boundary_charge_imbalance'],1e-8)
            self.assertGreater((out/'report.png').stat().st_size,10000)
            html=(out/'report.html').read_text()
            self.assertNotIn('<script src="https://cdn',html)
            self.assertTrue((out/'trajectory.vtp').is_file())
            with self.assertRaises(ValueError):run_case(config,out)
