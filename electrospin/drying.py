"""Drying along axial distance and actual parcel travel, distinct from jet arclength."""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np


def parcel_history(data,bead_id=None):
    if not len(data):return data,np.empty(0)
    bead_id=np.min(data['bead_id']) if bead_id is None else bead_id
    a=np.sort(data[data['bead_id']==bead_id],order='time_s')
    # New core integrates travel at each mechanical step. Old archives have
    # only sparse snapshots: use chord-sum reconstruction, explicitly labelled.
    if 'path_length_m' in a.dtype.names:s=a['path_length_m']
    else:
        p=np.column_stack([a[k] for k in ['x_m','y_m','z_m']]);s=np.r_[np.linalg.norm(p[0]),np.linalg.norm(np.diff(p,axis=0),axis=1)].cumsum()
    return a,s


def drying_values(a,c):
    ratio=a['solvent_mass_kg']/a['polymer_mass_kg']*c['polymer_fraction']/(1-c['polymer_fraction'])
    concentration=a['polymer_mass_kg']/(a['polymer_mass_kg']+a['solvent_mass_kg'])
    return ratio,concentration


def write_drying(out: Path,data,c):
    a,s=parcel_history(data)
    if not len(a):return
    retention,concentration=drying_values(a,c)
    with (out/'drying_history.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['time_s','bead_id','z_m','path_length_m','remaining_solvent_fraction','polymer_mass_fraction','temperature_K','radius_m'])
        w.writerows(zip(a['time_s'],a['bead_id'],a['z_m'],s,retention,concentration,a['temperature_K'],a['radius_m']))


def summarize_drying(data,c):
    a,s=parcel_history(data)
    if not len(a):return {}
    retention,_=drying_values(a,c)
    hit=np.flatnonzero(a['deposited']>=.5);dry=np.flatnonzero(retention<=.05)
    result={'drying_transport_model':c['transport_model'],
            'drying_remaining_solvent_fraction_last':float(retention[-1]),
            'drying_first_parcel_reached_collector':bool(len(hit)),
            'drying_remaining_solvent_fraction_at_first_impact':float(retention[hit[0]]) if len(hit) else None,
            'drying_95_percent_loss_z_m':float(a['z_m'][dry[0]]) if len(dry) else None,
            'drying_95_percent_loss_path_m':float(s[dry[0]]) if len(dry) else None,
            'drying_experimentally_validated':False}
    live=data[data['deposited']<.5]
    if len(live) and 'mass_transfer_m_s' in live.dtype.names:
        result.update(drying_mass_transfer_min_m_s=float(live['mass_transfer_m_s'].min()),
                      drying_mass_transfer_max_m_s=float(live['mass_transfer_m_s'].max()),
                      drying_heat_transfer_min_W_m2_K=float(live['heat_transfer_W_m2_K'].min()),
                      drying_heat_transfer_max_W_m2_K=float(live['heat_transfer_W_m2_K'].max()),
                      drying_crossflow_out_of_domain_sample_fraction=float(np.mean(live['transfer_correlation_in_domain']<.5)))
    return result
