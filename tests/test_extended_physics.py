"""Exercise the compiled solver: new closures, conservation, seeds and failure paths."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import test_core
from electrospin.config import validate,DEFAULTS
from electrospin.materials import pva_water,save_material,load_material,LEGACY_MATERIAL_KEYS


class ExtendedPhysics(unittest.TestCase):
    run_case=test_core.CoreIntegration.run_case

    def test_fene_force_free_has_no_artificial_stress(self):
        c=test_core.CoreIntegration.neutral();c.update(constitutive_model='fene_p')
        a,_,_=self.run_case(c)
        self.assertLess(max(abs(r['stress_Pa']) for r in a),1e-5)
        self.assertLess(max(abs(r['radius_m']/3e-6-1) for r in a),1e-8)

    def test_fene_finite_extensibility_changes_actual_solution(self):
        c=dict(duration=.004,output_interval=.001,evaporation=False,constitutive_model='fene_p')
        finite,_,_=self.run_case({**c,'finite_extensibility':10.})
        loose,_,_=self.run_case({**c,'finite_extensibility':1000.})
        self.assertGreater(abs(finite[-1]['stress_Pa']-loose[-1]['stress_Pa']),.01)

    def test_colored_inlet_is_reproducible_and_seed_dependent(self):
        c=dict(duration=.004,output_interval=.0002,perturbation_mode='colored_noise',random_seed=42)
        a,_,_=self.run_case(c);b,_,_=self.run_case(c);d,_,_=self.run_case({**c,'random_seed':43})
        self.assertEqual(a,b)
        self.assertGreater(max(abs(x['x_m']-y['x_m']) for x,y in zip(a,d)),1e-8)
        newest=[r for r in a if abs(r['z_m']-.000982437920320341)<1e-10]
        norms=np.array([np.hypot(r['x_m'],r['y_m']) for r in newest])
        self.assertGreater(np.std(norms),1e-8) # not a fixed-radius rotating inlet

    def test_conduction_changes_charges_and_preserves_polarity_symmetry(self):
        c=dict(duration=.004,output_interval=.001,charge_transport=True,conductivity=2e-4)
        a,_,_=self.run_case(c);b,_,_=self.run_case({**c,'conductivity':0.})
        self.assertGreater(max(abs(x['charge_C']-y['charge_C']) for x,y in zip(a,b)),1e-15)
        negative,_,_=self.run_case({**c,'voltage':-10000.,'current':-2e-8})
        self.assertEqual(len(a),len(negative))
        for x,y in zip(a,negative):
            self.assertAlmostEqual(x['z_m'],y['z_m'],delta=1e-11)
            self.assertAlmostEqual(x['charge_C'],-y['charge_C'],delta=1e-20)

    def test_cylinder_drag_responds_to_air_viscosity(self):
        c=dict(duration=.003,output_interval=.001,drag_model='cylinder')
        a,_,_=self.run_case(c);b,_,_=self.run_case({**c,'air_viscosity':3.7e-5})
        self.assertGreater(abs(a[-1]['vz_m_s']-b[-1]['vz_m_s']),1e-4)

    def test_all_closures_run_together_with_grounded_collector(self):
        c=dict(duration=.01,gap=.03,output_interval=.001,constitutive_model='fene_p',
               perturbation_mode='colored_noise',drag_model='cylinder',charge_transport=True,conductivity=2e-4,surface_contact=True)
        a,d,hit=self.run_case(c)
        self.assertGreater(len(hit),3)
        self.assertGreater(d[-1]['deposited'],3)
        self.assertGreater(max(r['stress_Pa'] for r in a),0)

    def test_surface_contact_precedes_centerline_crossing(self):
        c=test_core.CoreIntegration.neutral();c.update(gap=.005,duration=.002,viscosity=1e-12,surface_contact=True)
        _,_,hits=self.run_case(c)
        speed=1.388888888888889e-10/(np.pi*(3e-6)**2)
        self.assertAlmostEqual(min(h['time_s'] for h in hits),(.005-3e-6)/speed,delta=1e-10)

    def test_bad_new_parameters_are_rejected_in_core(self):
        self.run_case(dict(test_core.CoreIntegration.neutral(),finite_extensibility=2),expect_ok=False)
        self.run_case(dict(test_core.CoreIntegration.neutral(),perturbation_mode='fake_spiral'),expect_ok=False)


class MaterialValidation(unittest.TestCase):
    def test_old_material_retains_its_original_constant_transfer_law(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'old.json';c=pva_water()
            p.write_text(json.dumps({'schema_version':1,'name':'old','properties':{k:c[k] for k in LEGACY_MATERIAL_KEYS}}))
            restored=load_material(p,c)
            self.assertEqual(restored['transport_model'],'constant')
            self.assertFalse(restored['liquid_diffusion'])

    def test_material_round_trip_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'BN.json';c=pva_water();c['viscosity']=.7
            save_material(c,'ПВС с BN','Измерение пользователя','BN 0.5%; свойства раствора',p)
            b=load_material(p,pva_water());self.assertEqual(b['viscosity'],.7)
            self.assertEqual(b['_metadata']['material_name'],'ПВС с BN')
            with self.assertRaises(FileExistsError):save_material(c,'x','','',p)

    def test_validate_rejects_invalid_extended_parameters(self):
        for patch in [dict(finite_extensibility=3),dict(noise_correlation_time=0),dict(random_seed=0),dict(coulomb_softening=0)]:
            with self.subTest(patch=patch),self.assertRaises(ValueError):validate({**DEFAULTS,**patch})
