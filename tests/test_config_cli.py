from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest
from electrospin.config import CORE_DEFAULTS, DEFAULTS, ROOT, load_config, validate, write_namelist


class ConfigurationTests(unittest.TestCase):
    def test_defaults_match_fortran_manifest(self):
        self.assertEqual(CORE_DEFAULTS,json.loads((ROOT/'config/core_defaults.json').read_text()))

    def test_invalid_physics_or_typos_rejected(self):
        for changes in ({'polymer_fraction':12}, {'humidity':40}, {'needle_inner_radius':.001},
                        {'volatge':3000}, {'duration':float('nan')}, {'evaporation':'false'},
                        {'field_nr':10}, {'max_beads':3}, {'gravity':[1,2]}, {'dt_max':-1}):
            with self.subTest(changes=changes),self.assertRaises(ValueError): validate(changes)

    def test_overrides_are_json_and_si(self):
        config = load_config(overrides=['voltage=12000','collector_velocity=[0.1,0]','evaporation=false'])
        self.assertEqual(config['voltage'],12000)
        self.assertEqual(config['collector_velocity'],[.1,0])
        self.assertFalse(config['evaporation'])
        with self.assertRaises(ValueError): load_config(overrides=['gap=10cm'])

    def test_namelist_excludes_field_python_keys_and_quotes_paths(self):
        c = validate({'field_file':"/tmp/a'b.csv"})
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'case.nml';write_namelist(c,path)
            text=path.read_text()
            self.assertNotIn('field_mode',text)
            self.assertIn("field_file = '/tmp/a''b.csv'",text)
            self.assertIn('evaporation = .true.',text)

    def test_existing_results_preserved(self):
        from electrospin.cli import run_case
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);p=out/'important.txt';p.write_text('keep')
            with self.assertRaisesRegex(ValueError,'not empty'):
                run_case(DEFAULTS,out)
            self.assertEqual(p.read_text(),'keep')

    def test_terminal_wizard_unit_conversions(self):
        from electrospin.cli import PROMPTS
        fields = {key:(scale,offset) for key,_,scale,offset in PROMPTS}
        scale,offset=fields['emitter_radius'];self.assertAlmostEqual(.70*scale+offset,.00035)
        scale,offset=fields['flow_rate'];self.assertAlmostEqual(.5*scale+offset,1.388888888888889e-10)
        scale,offset=fields['temperature'];self.assertAlmostEqual(25*scale+offset,298.15)


if __name__=='__main__':unittest.main()
