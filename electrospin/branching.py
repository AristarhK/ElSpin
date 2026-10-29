"""Separate, uncalibrated population balance; NOT daughter trajectories or CFD.

Parent speed is prescribed from the first parcel of the single-filament run.
Equal splitting preserves polymer and charge fluxes; solvent loss and heat
are reintegrated because total interface grows with the expected count N.
The kinetic coefficient must be calibrated; the stress ratio alone gives no
experimentally identified branching rate. N is an expectation, not an integer.
"""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from .drying import parcel_history
from .transport import evaporation


def axial_profile(data,branching):
    """Map scenario samples to the parent's actual z, never replace z by travel."""
    if not len(branching):return np.empty(0),np.empty(0)
    if 'z_m' in branching.dtype.names:
        return branching['z_m'],branching['expected_filaments']
    history,travel=parcel_history(data)
    if not len(travel):return np.empty(0),np.empty(0)
    unique=np.r_[True,np.diff(travel)>1e-12]
    return np.interp(branching['path_length_m'],travel[unique],history['z_m'][unique]),branching['expected_filaments']


def population_balance(config,s,speed):
    c=config;s=np.asarray(s,float);speed=np.asarray(speed,float)
    if len(s)<2 or np.any(np.diff(s)<=0) or np.any(speed<=0):raise ValueError('Need increasing travel and positive speed')
    w=c['polymer_fraction'];rho=1/(w/c['density_polymer']+(1-w)/c['density_solvent'])
    mp=rho*c['flow_rate']*w;ms0=rho*c['flow_rate']*(1-w);current=c['current']
    ambient=c['ambient_temperature'];R=8.314462618;eps0=8.8541878128e-12
    psat_air=c['psat_ref']*np.exp(c['latent_heat']*c['molar_mass']/R*(1/c['psat_temperature']-1/ambient))
    def properties(x,y):
        N=np.exp(y[0]);ms=ms0*np.exp(y[1]);T=y[2];u=np.interp(x,s,speed)
        volume=mp/c['density_polymer']+ms/c['density_solvent'];a=np.sqrt(volume/(N*np.pi*u))
        sigma=current/(N*u*2*np.pi*a)
        ratio=sigma*sigma*a/(2*eps0*max(c['surface_tension'],1e-30))
        return N,ms,T,u,volume,a,ratio
    def rhs(x,y):
        N,ms,T,u,volume,a,ratio=properties(x,y)
        rate=(c['branching_rate_s_inv']*max(0.,ratio-c['branching_threshold'])/u if c['branching_enabled'] else 0.)
        activity=(ms/c['density_solvent'])/volume
        flux,h,_,_,_,_=evaporation(c,a,T,activity,0.)
        loss=(flux*N*2*np.pi*a if c['evaporation'] else 0.)
        heat=h*N*2*np.pi*a*(ambient-T)-c['latent_heat']*loss
        return [rate,-loss/max(ms,1e-300),heat/(mp*c['cp_polymer']+ms*c['cp_solvent'])]
    def limit(x,y):return np.log(c['branching_max_expected'])-y[0]
    limit.terminal=True;limit.direction=-1
    solution=solve_ivp(rhs,(s[0],s[-1]),[0.,0.,c['temperature']],t_eval=s,rtol=1e-6,atol=[1e-8,1e-8,1e-5],events=limit)
    if not solution.success:raise RuntimeError('Population scenario: '+solution.message)
    rows=[]
    for x,y in zip(solution.t,solution.y.T):
        N,ms,T,u,volume,a,ratio=properties(x,y)
        rows.append([x,N,mp,ms,current,mp/N,ms/N,current/N,a,T,ratio])
    return np.asarray(rows),('LIMIT_REACHED' if solution.status==1 else 'COMPLETED_SCENARIO')


def run_scenario(out:Path,data,c):
    if not c['branching_enabled']:
        return {'branching_scenario_status':'DISABLED','branching_scenario_enabled':False,
                'branching_expected_count_at_last_point':None,
                'branching_scenario_warning':'No physical daughter-count prediction. Disabled scenario is not a calculated N(z)=1 branching result.'}
    if c['transport_model']=='crossflow':
        return {'branching_scenario_status':'UNSUPPORTED_CROSSFLOW_PROFILE','branching_scenario_enabled':True,
                'branching_scenario_warning':'Parent normal-velocity profile is not transferred to the population scenario; no count is reported.'}
    a,s=parcel_history(data);mask=np.r_[True,np.diff(s)>1e-12] if len(s) else []
    a=a[mask];s=s[mask]
    if len(s)<2:return {'branching_scenario_status':'INSUFFICIENT_HISTORY'}
    speed=np.sqrt(a['vx_m_s']**2+a['vy_m_s']**2+a['vz_m_s']**2)
    keep=(a['deposited']<.5)&(speed>0);a=a[keep];s=s[keep];speed=speed[keep]
    if len(s)<2:return {'branching_scenario_status':'INSUFFICIENT_HISTORY'}
    values,status=population_balance(c,s,speed)
    with (out/'branching_scenario.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['path_length_m','expected_filaments','polymer_flux_kg_s','solvent_flux_kg_s','charge_flux_A','polymer_flux_per_filament_kg_s','solvent_flux_per_filament_kg_s','charge_flux_per_filament_A','radius_per_filament_m','temperature_K','electric_capillary_ratio','z_m'])
        w.writerows([*row,float(np.interp(row[0],s,a['z_m']))] for row in values)
    return {'branching_scenario_status':status,'branching_scenario_enabled':c['branching_enabled'],
            'branching_expected_count_at_last_point':float(values[-1,1]),
            'branching_scenario_warning':'Uncalibrated equal-split population balance with prescribed parent velocity; no branch trajectories, no daughter mutual forces and no feedback to the single-filament 3D run. N is an expectation. A nonzero kinetic coefficient requires experimental calibration.'}
