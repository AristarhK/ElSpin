"""SI heat/mass-transfer reference, independent of the mechanical integrator.

gas_film: exact steady radial diffusion/conduction between a cylinder of
radius a and a prescribed well-mixed boundary at R. The boundary radius
is an explicit uncertain physical assumption, not a hidden fitted k_m.
crossflow: Churchill--Bernstein (1977), used only with a domain diagnostic.
Liquid resistance is an optional quasi-steady parabolic radial closure.
"""
from __future__ import annotations
import math

R_GAS=8.314462618


def churchill_bernstein(re,pr):
    if re<0 or pr<=0:raise ValueError('Re>=0 and Pr/Sc>0 required')
    return .3+.62*math.sqrt(re)*pr**(1/3)/(1+(.4/pr)**(2/3))**.25*(1+(re/282000.)**.625)**.8


def coefficients(c,radius,normal_speed):
    if radius<=0 or normal_speed<0:raise ValueError('Invalid radius or speed')
    re=2*radius*c['air_density']*normal_speed/c['air_viscosity']
    if c['transport_model']=='constant':return c['mass_transfer'],c['heat_transfer'],re,True
    if c['transport_model']=='gas_film':
        if c['gas_boundary_radius']<=radius:raise ValueError('Gas boundary must be outside the jet')
        resistance=radius*math.log(c['gas_boundary_radius']/radius)
        return c['vapor_diffusivity']/resistance,c['air_thermal_conductivity']/resistance,re,True
    sc=c['air_viscosity']/(c['air_density']*c['vapor_diffusivity'])
    pr=c['air_viscosity']*c['air_heat_capacity']/c['air_thermal_conductivity']
    return (churchill_bernstein(re,sc)*c['vapor_diffusivity']/(2*radius),
            churchill_bernstein(re,pr)*c['air_thermal_conductivity']/(2*radius),
            re,re*min(pr,sc)>=.2)


def solvent_activity(c,phi):
    if phi<=0:return 0.
    if phi>=1:return 1.
    if c['activity_model']=='volume_fraction':return phi
    # High-polymerization Flory-Huggins limit; chi is a hypothesis to calibrate.
    return math.exp(math.log(phi)+(1-phi)+c['flory_huggins_chi']*(1-phi)**2)


def vapor_mass_fraction(c,pressure):
    if not 0<=pressure<c['ambient_pressure']:raise ValueError('Vapor pressure outside dilute-gas domain')
    x=pressure/c['ambient_pressure']
    return c['molar_mass']*x/(c['molar_mass']*x+c['air_molar_mass']*(1-x))


def evaporation(c,radius,temperature,phi,normal_speed):
    km,h,re,valid=coefficients(c,radius,normal_speed)
    ps=c['psat_ref']*math.exp(c['latent_heat']*c['molar_mass']/R_GAS*(1/c['psat_temperature']-1/temperature))
    pa=c['psat_ref']*math.exp(c['latent_heat']*c['molar_mass']/R_GAS*(1/c['psat_temperature']-1/c['ambient_temperature']))
    activity=solvent_activity(c,phi)
    if c['transport_model']=='constant':
        flux=max(0.,km*c['molar_mass']/R_GAS*(activity*ps/temperature-c['humidity']*pa/c['ambient_temperature']))
        return flux,h,km,re,valid,phi
    ya=vapor_mass_fraction(c,c['humidity']*pa)
    def gas_flux(surface_phi):
        ys=vapor_mass_fraction(c,solvent_activity(c,surface_phi)*ps)
        return max(0.,c['air_density']*km*math.log((1-ya)/(1-ys)))
    surface_phi=phi;flux=gas_flux(phi)
    if c['liquid_diffusion'] and flux>0 and phi>0:
        phi0=(1-c['polymer_fraction'])/c['density_solvent']/(
            c['polymer_fraction']/c['density_polymer']+(1-c['polymer_fraction'])/c['density_solvent'])
        diffusivity=c['liquid_diffusivity']*(phi/phi0)**c['diffusion_concentration_exponent']
        conductance=4*c['density_solvent']*diffusivity/radius
        lo,hi=0.,phi
        for _ in range(45):
            mid=(lo+hi)/2
            if gas_flux(mid)>conductance*(phi-mid):hi=mid
            else:lo=mid
        surface_phi=(lo+hi)/2;flux=gas_flux(surface_phi)
    return flux,h,km,re,valid,surface_phi
