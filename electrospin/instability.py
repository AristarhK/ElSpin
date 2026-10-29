"""Diagnostics of simulated transverse motion, without manufacturing a trajectory."""
from __future__ import annotations
import csv
import json
from pathlib import Path
import numpy as np
from .analysis import read_table


def analyze_instability(out: Path) -> dict:
    data = read_table(out/'trajectory.csv')
    cfg = json.loads((out/'config.json').read_text(encoding='utf-8'))
    times = np.unique(data['time_s'])
    history = []
    for t in times:
        live = data[(data['time_s']==t) & (data['deposited']<.5)]
        if not len(live): continue
        mass = live['mass_kg']; r2 = live['x_m']**2+live['y_m']**2
        # Mass-weighted RMS; a mean-free measure distinguishes translation.
        cx = np.average(live['x_m'],weights=mass); cy=np.average(live['y_m'],weights=mass)
        rms = np.sqrt(np.average(r2,weights=mass))
        centered = np.sqrt(np.average((live['x_m']-cx)**2+(live['y_m']-cy)**2,weights=mass))
        kinetic = .5*np.sum(mass*(live['vx_m_s']**2+live['vy_m_s']**2))
        history.append([t,rms,centered,np.sqrt(r2.max()),kinetic])
    with (out/'instability.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['time_s','rms_transverse_m','centered_rms_m','max_transverse_m','transverse_kinetic_J']);w.writerows(history)
    # Growth vs z in equal bins, with counts to expose sparse sampling.
    bins=np.linspace(0,cfg['gap'],31); rows=[]
    live=data[data['deposited']<.5]
    for lo,hi in zip(bins[:-1],bins[1:]):
        part=live[(live['z_m']>=lo)&(live['z_m']<hi)]
        if len(part): rows.append([.5*(lo+hi),np.sqrt(np.mean(part['x_m']**2+part['y_m']**2)),len(part)])
    with (out/'bending_profile.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['z_m','rms_transverse_m','samples']);w.writerows(rows)
    result={'inlet_seed_mode':cfg.get('perturbation_mode','periodic'), 'inlet_random_seed':cfg.get('random_seed',42),
            'instability_interpretation':'Transverse amplification diagnostic. A Coulomb on/off control, several seeds, and discretization refinement are required before claiming a bending instability.',
            'rms_transverse_gain':None, 'apparent_growth_rate_s_inv':None}
    if history:
        h=np.asarray(history); seed=cfg['perturbation_amplitude']
        if seed>0: result['rms_transverse_gain']=float(h[-1,1]/seed)
        # Not called a linear eigenvalue: the population changes during startup.
        select=(h[:,0]>.2*cfg['duration'])&(h[:,0]<.7*cfg['duration'])&(h[:,1]>0)
        if np.count_nonzero(select)>=5:
            result['apparent_growth_rate_s_inv']=float(np.polyfit(h[select,0],np.log(h[select,1]),1)[0])
    if len(live):
        a=live['radius_m'];ell=live['segment_length_m']
        result['minimum_active_segment_slenderness']=float(np.min(ell/a))
        if cfg['surface_tension']>0:
            sigma=live['charge_C']/(2*np.pi*a*ell)
            ratio=sigma**2*a/(2*8.8541878128e-12*cfg['surface_tension'])
            result['maximum_electrocapillary_pressure_ratio']=float(np.max(ratio))
            if np.max(ratio)>1:
                result['radial_stress_warning']='Surface-charge pressure exceeds capillary pressure in some cells. Radial stress, fission and breakup are not resolved; a completed run is not proof of a stable continuous jet.'
    return result
