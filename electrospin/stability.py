"""Local cross-section stability SCREEN, never a daughter-count prediction.

Conducting cylinder, azimuthal mode m>=2, k_z=0, negligible exterior stress.
In the Stokes/Newtonian limit, stress balance gives
 s = gamma*m/(2*eta*a*(m+1)) * (B_e-m-1).
This must not be interpreted as a polymer-jet branching rate outside its
Ohnesorge, Deborah and rapid-charge-relaxation limits. Neither sign nor
m identifies how many daughter jets will actually form.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from .drying import parcel_history
EPS0=8.8541878128e-12


def stokes_azimuthal_rate(radius,density,viscosity,surface_tension,linear_charge,mode=2):
    if radius<=0 or density<=0 or viscosity<=0 or surface_tension<=0 or mode<2:
        raise ValueError('Positive physical properties and azimuthal mode m>=2 required')
    electric_bond=linear_charge**2/(4*np.pi**2*EPS0*surface_tension*radius)
    rate=surface_tension*mode/(2*viscosity*radius*(mode+1))*(electric_bond-mode-1)
    return float(rate),float(electric_bond)


def screening_rows(data,c):
    a,s=parcel_history(data)
    if not len(a) or c['surface_tension']<=0:return np.empty((0,11))
    result=[]
    for row,travel in zip(a,s):
        if row['deposited']>=.5:break
        radius=row['radius_m'];w=row['polymer_mass_kg']/row['mass_kg']
        eta=c['viscosity']*(w/c['polymer_fraction'])**c['viscosity_exponent']*np.exp(c['activation_energy']/8.314462618*(1/row['temperature_K']-1/c['temperature']))
        rho=1/(w/c['density_polymer']+(1-w)/c['density_solvent'])
        linear=row['charge_C']/row['segment_length_m']
        rate,bond=stokes_azimuthal_rate(radius,rho,eta,c['surface_tension'],linear)
        oh=eta/np.sqrt(rho*c['surface_tension']*radius)
        de=abs(rate)*c['relaxation_time']
        charge=abs(rate)*EPS0*c['relative_permittivity']/c['conductivity'] if c['conductivity']>0 else 1e300
        prestress=abs(row['stress_Pa'])*radius/c['surface_tension']
        valid=oh>=10 and de<=.1 and charge<=.1 and prestress<=.1
        result.append([row['time_s'],row['z_m'],travel,bond/2,rate,oh,de,charge,prestress,int(valid),linear])
    return np.asarray(result)


def write_screening(out:Path,data,c):
    values=screening_rows(data,c)
    with (out/'branching_screen.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['time_s','z_m','path_length_m','electric_capillary_pressure_ratio',
             'm2_stokes_conductor_limit_rate_s_inv','Ohnesorge','Deborah_of_limit_rate','charge_relaxation_parameter',
             'axial_stress_capillary_ratio','limit_model_in_domain','linear_charge_C_m']);w.writerows(values)
    return {'physical_branch_count_prediction_available':False,
            'physical_branch_count_status':'NOT_IDENTIFIED: topology-changing kinetics or surface-resolved model and validation data are missing.',
            'branching_screen_in_domain_fraction':float(values[:,9].mean()) if len(values) else None,
            'branching_screen_pressure_ratio_max':float(values[:,3].max()) if len(values) else None,
            'branching_screen_note':'Local conducting-cylinder m=2 neutral ratio=1.5 is a limiting stability criterion, not a daughter count or a calibrated branching law. Rate is only a Newtonian Stokes conductor limit; check all exported applicability flags.'}
