"""Independent steady translation benchmark: no electric or mechanical stretching.

A unit material length travels at prescribed speed. Polymer mass/length is
constant; solvent loss shrinks its radius. This isolates the transport model
from the unresolved spatial convergence of the 3D bending discretization.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from .transport import evaporation


def integrate_filament(c,radius0,speed,z,*,initial_mass_fraction=None,initial_temperature=None,
                       rtol=1e-8,atol=1e-9):
    z=np.asarray(z,float)
    if c['transport_model']=='crossflow':raise ValueError('Crossflow reference requires a specified normal gas velocity profile')
    if radius0<=0 or speed<=0 or len(z)<2 or np.any(np.diff(z)<=0):raise ValueError('Positive radius/speed and increasing z required')
    w=c['polymer_fraction'] if initial_mass_fraction is None else initial_mass_fraction
    if not 0<w<1:raise ValueError('Initial mass fraction must be in (0,1)')
    total=np.pi*radius0**2/(w/c['density_polymer']+(1-w)/c['density_solvent'])
    mp=w*total;ms0=(1-w)*total
    temperature=c['temperature'] if initial_temperature is None else initial_temperature
    def rhs(x,y):
        ms=ms0*np.exp(y[0]);volume=mp/c['density_polymer']+ms/c['density_solvent'];a=np.sqrt(volume/np.pi)
        phi=(ms/c['density_solvent'])/volume
        flux,h,_,_,_,_=evaporation(c,a,y[1],phi,0.)
        loss=2*np.pi*a*flux if c['evaporation'] else 0.
        heat=2*np.pi*a*h*(c['ambient_temperature']-y[1])-c['latent_heat']*loss
        return [-loss/(ms*speed),heat/((mp*c['cp_polymer']+ms*c['cp_solvent'])*speed)]
    sol=solve_ivp(rhs,(z[0],z[-1]),[0.,temperature],t_eval=z,rtol=rtol,atol=atol,method='Radau')
    if not sol.success:raise RuntimeError(sol.message)
    ms=ms0*np.exp(sol.y[0]);radius=np.sqrt((mp/c['density_polymer']+ms/c['density_solvent'])/np.pi)
    return {'z_m':sol.t,'time_s':(sol.t-z[0])/speed,'remaining_solvent_fraction':np.exp(sol.y[0]),
            'temperature_K':sol.y[1],'radius_m':radius,'polymer_mass_per_length_kg_m':np.full(len(z),mp),
            'solvent_mass_per_length_kg_m':ms,'polymer_mass_fraction':mp/(mp+ms),
            'inlet_flow_ml_h':np.pi*radius0**2*speed*3600e6}
