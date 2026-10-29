import unittest
import numpy as np
from electrospin.pump import Pump
from electrospin.branching import population_balance
from electrospin.materials import pva_water


class PumpTests(unittest.TestCase):
    def test_volume_per_step_and_frequency_reproduce_requested_flow(self):
        p=Pump(syringe_inner_diameter_mm=20,screw_lead_mm_rev=2,motor_steps_rev=200,microsteps=32,flow_ml_h=.5)
        r=p.calculate()
        self.assertAlmostEqual(r['step_frequency_Hz']*r['volume_per_step_ml']*3600,.5)
        self.assertAlmostEqual(r['plunger_speed_mm_s']*np.pi*100*3600/1000,.5)
        self.assertGreater(r['pressure_drop_Pa'],0)

    def test_hardware_microstepping_limit_is_enforced(self):
        with self.assertRaises(ValueError):Pump(driver='DRV8825',microsteps=64).calculate()
        with self.assertRaises(ValueError):Pump(driver='A4988',microsteps=32).calculate()

    def test_syringe_diameter_changes_feed_calibration(self):
        a=Pump(syringe_inner_diameter_mm=10).calculate();b=Pump(syringe_inner_diameter_mm=20).calculate()
        self.assertAlmostEqual(a['step_frequency_Hz']/b['step_frequency_Hz'],4)
        self.assertAlmostEqual(b['plunger_hydraulic_force_N']/a['plunger_hydraulic_force_N'],4)


class PopulationTests(unittest.TestCase):
    def test_no_branching_and_no_evaporation_preserve_flux(self):
        c=pva_water();c.update(evaporation=False,heat_transfer=0,branching_enabled=False)
        a,status=population_balance(c,np.linspace(.001,.1,30),np.full(30,5.))
        self.assertEqual(status,'COMPLETED_SCENARIO');np.testing.assert_allclose(a[:,1],1)
        np.testing.assert_allclose(a[:,3],a[0,3]);np.testing.assert_allclose(a[:,9],c['temperature'])

    def test_equal_branching_conserves_polymer_solvent_and_charge_fluxes(self):
        c=pva_water();c.update(evaporation=False,heat_transfer=0,branching_enabled=True,current=2e-7,branching_rate_s_inv=20.)
        a,status=population_balance(c,np.linspace(.001,.1,50),np.full(50,5.))
        self.assertGreater(a[-1,1],1.01)
        for total,per in [(2,5),(3,6),(4,7)]:np.testing.assert_allclose(a[:,1]*a[:,per],a[:,total],rtol=1e-12)
        np.testing.assert_allclose(a[:,1]*a[:,8]**2,a[0,8]**2,rtol=1e-6)

    def test_population_limit_terminates_without_clipping_count(self):
        c=pva_water();c.update(branching_enabled=True,current=2e-7,branching_rate_s_inv=1000.,branching_max_expected=1.1)
        a,status=population_balance(c,np.linspace(.001,.1,100),np.full(100,5.))
        self.assertEqual(status,'LIMIT_REACHED');self.assertLessEqual(a[-1,1],1.1)
