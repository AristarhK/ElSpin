import unittest
import numpy as np
from electrospin.stability import stokes_azimuthal_rate,EPS0
from electrospin.materials import pva_water
from electrospin.branching import run_scenario
from pathlib import Path
import tempfile


class StabilityLimits(unittest.TestCase):
    def test_uncharged_stokes_cylinder_relaxes(self):
        a=2e-6;eta=.2;gamma=.055
        for m in [2,3,4]:
            rate,bond=stokes_azimuthal_rate(a,1000,eta,gamma,0,m)
            self.assertEqual(bond,0)
            self.assertAlmostEqual(rate,-gamma*m/(2*eta*a))

    def test_neutral_electric_bond_is_mode_plus_one(self):
        a=2e-6;gamma=.055
        for m in [2,3,4]:
            linear=np.sqrt((m+1)*4*np.pi**2*EPS0*gamma*a)
            rate,bond=stokes_azimuthal_rate(a,1000,.2,gamma,linear,m)
            self.assertAlmostEqual(bond,m+1)
            self.assertAlmostEqual(rate,0,delta=1e-8)
            self.assertGreater(stokes_azimuthal_rate(a,1000,.2,gamma,linear*1.1,m)[0],0)
            self.assertLess(stokes_azimuthal_rate(a,1000,.2,gamma,linear*.9,m)[0],0)

    def test_viscosity_changes_rate_not_neutral_threshold(self):
        r1,b1=stokes_azimuthal_rate(2e-6,1000,.2,.055,1e-8)
        r2,b2=stokes_azimuthal_rate(2e-6,1000,.4,.055,1e-8)
        self.assertAlmostEqual(r1,2*r2);self.assertEqual(b1,b2)

    def test_disabled_scenario_does_not_report_a_fake_physical_count(self):
        with tempfile.TemporaryDirectory() as d:
            s=run_scenario(Path(d),np.empty(0),pva_water())
            self.assertIsNone(s['branching_expected_count_at_last_point'])
            self.assertEqual(s['branching_scenario_status'],'DISABLED')
            self.assertFalse((Path(d)/'branching_scenario.csv').exists())
