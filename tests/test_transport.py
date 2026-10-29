import math
import unittest
import numpy as np
from electrospin.config import DEFAULTS,validate
from electrospin.transport import coefficients,evaporation,solvent_activity
from electrospin.drying_reference import integrate_filament
from electrospin.materials import pva_water
import test_core


class TransportAnalytical(unittest.TestCase):
    def setUp(self):self.c=validate({**DEFAULTS,'transport_model':'gas_film'})

    def test_radial_annulus_analytical_conductance(self):
        # Total conductance per length is exactly 2 pi D / log(R/a).
        a=3e-6;km,h,_,valid=coefficients(self.c,a,0)
        self.assertAlmostEqual(2*math.pi*a*km,2*math.pi*self.c['vapor_diffusivity']/math.log(.001/a))
        self.assertAlmostEqual(2*math.pi*a*h,2*math.pi*self.c['air_thermal_conductivity']/math.log(.001/a))
        self.assertTrue(valid)

    def test_radius_and_gas_domain_have_physical_effects(self):
        km,_,_,_=coefficients(self.c,3e-6,0);fine=coefficients(self.c,1e-6,0)[0]
        far=coefficients({**self.c,'gas_boundary_radius':.005},3e-6,0)[0]
        self.assertGreater(fine,km);self.assertLess(far,km)

    def test_saturation_and_absence_of_solvent_give_zero_loss(self):
        self.assertAlmostEqual(evaporation({**self.c,'humidity':1},3e-6,298.15,1,0)[0],0)
        self.assertEqual(evaporation(self.c,3e-6,298.15,0,0)[0],0)

    def test_liquid_diffusion_limits_surface_flux(self):
        c={**self.c,'liquid_diffusion':True,'liquid_diffusivity':1e-12,'diffusion_concentration_exponent':0}
        flux,_,_,_,_,surface=evaporation(c,3e-6,298.15,.9,0)
        self.assertLess(surface,.9);self.assertLessEqual(flux,4*c['density_solvent']*1e-12/(3e-6)*.9)
        self.assertLess(flux,evaporation({**c,'liquid_diffusion':False},3e-6,298.15,.9,0)[0])

    def test_flory_huggins_pure_solvent_and_dry_limits(self):
        c={**self.c,'activity_model':'flory_huggins'}
        self.assertEqual(solvent_activity(c,0),0);self.assertEqual(solvent_activity(c,1),1)
        values=[solvent_activity(c,x) for x in np.linspace(.01,.99,100)]
        self.assertTrue(np.all(np.diff(values)>0));self.assertLess(max(values),1)

    def test_crossflow_applicability_is_not_silently_assumed(self):
        c={**self.c,'transport_model':'crossflow'}
        self.assertFalse(coefficients(c,3e-6,0)[3]);self.assertTrue(coefficients(c,3e-6,10)[3])


class CompiledTransport(unittest.TestCase):
    run_case=test_core.CoreIntegration.run_case

    def test_actual_compiled_flux_matches_independent_annulus_model(self):
        for resistance in [False,True]:
            c={**test_core.CoreIntegration.neutral(),'transport_model':'gas_film',
               'activity_model':'flory_huggins','liquid_diffusion':resistance,'evaporation':True}
            a,_,_=self.run_case(c)
            first=a[0];full=validate(c)
            phi=(first['solvent_mass_kg']/full['density_solvent'])/(first['solvent_mass_kg']/full['density_solvent']+first['polymer_mass_kg']/full['density_polymer'])
            reference=evaporation(full,first['radius_m'],first['temperature_K'],phi,0)
            self.assertAlmostEqual(first['mass_transfer_m_s'],reference[2],delta=1e-10)
            self.assertAlmostEqual(first['heat_transfer_W_m2_K'],reference[1],delta=1e-7)
            self.assertAlmostEqual(first['evaporation_flux_kg_m2_s'],reference[0],delta=1e-11)
            self.assertAlmostEqual(first['surface_solvent_volume_fraction'],reference[5],delta=1e-10)

    def test_adiabatic_latent_cooling_matches_closed_form(self):
        c={**test_core.CoreIntegration.neutral(),'duration':.001,'evaporation':True,
           'mass_transfer':2.,'dt_max':1e-7,'output_interval':.0001}
        a,_,_=self.run_case(c)
        history=[r for r in a if r['bead_id']==1 and r['time_s']>=.0002-1e-12]
        first,last=history[0],history[-1];mp=first['polymer_mass_kg']
        exact=first['temperature_K']+2.44e6/4180*math.log((mp*1700+last['solvent_mass_kg']*4180)/(mp*1700+first['solvent_mass_kg']*4180))
        self.assertAlmostEqual(last['temperature_K'],exact,delta=.002)

    def test_inlet_noise_path_is_shared_under_injection_refinement(self):
        paths=[]
        for interval in [.0002,.0001]:
            c={**test_core.CoreIntegration.neutral(),'perturbation_mode':'colored_noise',
               'perturbation_amplitude':1e-7,'injection_interval':interval,'output_interval':.0002}
            a,_,_=self.run_case(c);path={}
            for r in a:
                if abs(r['z_m']-1.388888888888889e-10/(math.pi*(3e-6)**2)*interval)<1e-11:
                    path[round(r['time_s'],8)]=(r['x_m'],r['y_m'])
            paths.append(path)
        common=set(paths[0])&set(paths[1]);self.assertGreater(len(common),2)
        for t in common:np.testing.assert_allclose(paths[0][t],paths[1][t],atol=1e-16,rtol=0)

    def test_compiled_mass_and_temperature_follow_independent_segment_ode(self):
        c={**test_core.CoreIntegration.neutral(),'duration':.001,'dt_max':1e-7,
           'output_interval':.0001,'transport_model':'gas_film','activity_model':'flory_huggins',
           'liquid_diffusion':True,'evaporation':True}
        data,_,_=self.run_case(c);full=validate(c)
        a=[r for r in data if r['bead_id']==1 and r['time_s']>=.0002-1e-12]
        initial,last=a[0],a[-1];speed=initial['vz_m_s']
        reference=integrate_filament(full,initial['radius_m'],speed,[0.,(last['time_s']-initial['time_s'])*speed],
            initial_mass_fraction=initial['polymer_mass_kg']/initial['mass_kg'],initial_temperature=initial['temperature_K'])
        self.assertAlmostEqual(last['solvent_mass_kg']/initial['solvent_mass_kg'],reference['remaining_solvent_fraction'][-1],delta=2e-5)
        self.assertAlmostEqual(last['temperature_K'],reference['temperature_K'][-1],delta=.002)


class IndependentDryingReference(unittest.TestCase):
    def test_tolerance_refinement_and_polymer_conservation(self):
        c=pva_water();z=np.linspace(0,.1,41)
        coarse=integrate_filament(c,1e-6,5.,z,rtol=1e-5)
        fine=integrate_filament(c,1e-6,5.,z,rtol=1e-8)
        self.assertLess(abs(coarse['remaining_solvent_fraction'][-1]-fine['remaining_solvent_fraction'][-1]),1e-4)
        np.testing.assert_allclose(fine['polymer_mass_per_length_kg_m'],fine['polymer_mass_per_length_kg_m'][0])
        self.assertTrue(np.all(np.diff(fine['remaining_solvent_fraction'])<=1e-9))

    def test_thinner_cylinder_drains_faster_at_same_prescribed_speed(self):
        c=pva_water();z=np.linspace(0,.1,41)
        thick=integrate_filament(c,3e-6,5.,z);thin=integrate_filament(c,.2e-6,5.,z)
        self.assertLess(thin['remaining_solvent_fraction'][-1],thick['remaining_solvent_fraction'][-1])
