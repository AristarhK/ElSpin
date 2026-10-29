"""Validated flat JSON configuration and strict Fortran namelist serialization."""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path

ROOT = Path(getattr(sys, '_MEIPASS', Path(__file__).resolve().parents[1]))
# Synchronized with the standalone Fortran solver; JSON examples may override these.
CORE_DEFAULTS = {
    'voltage': 10000.0, 'gap': 0.10, 'flow_rate': 1.388888888888889e-10,
    'jet_radius': 3e-6, 'density_polymer': 1200.0, 'density_solvent': 997.0,
    'polymer_fraction': 0.12, 'viscosity': 0.2, 'relaxation_time': 0.01,
    'surface_tension': 0.055, 'current': 2e-8, 'conductivity': 2e-6,
    'relative_permittivity': 40.0, 'temperature': 298.15, 'ambient_temperature': 298.15,
    'humidity': 0.4, 'mass_transfer': 0.01, 'heat_transfer': 100.0,
    'latent_heat': 2.44e6, 'molar_mass': 0.01801528, 'psat_ref': 3169.0,
    'psat_temperature': 298.15, 'cp_polymer': 1700.0, 'cp_solvent': 4180.0,
    'viscosity_exponent': 3.0, 'activation_energy': 20000.0,
    'air_density': 1.184, 'air_viscosity': 1.85e-5, 'drag_coefficient': 0.01,
    'air_velocity': [0.0, 0.0, 0.0], 'gravity': [0.0, 0.0, 9.81],
    'duration': 0.03, 'dt_max': 2e-6, 'dt_min': 1e-12,
    'injection_interval': 0.0002, 'output_interval': 0.0002,
    'max_beads': 500, 'max_steps': 2000000,
    'perturbation_amplitude': 1e-6, 'perturbation_frequency': 250.0,
    'collector_velocity': [0.0, 0.0], 'coulomb_softening': 3e-6,
    'evaporation': True, 'coulomb': True, 'image_force': False, 'field_file': '',
    'perturbation_mode': 'periodic', 'random_seed': 42, 'noise_correlation_time': 4e-4,
    'constitutive_model': 'maxwell', 'finite_extensibility': 100.,
    'drag_model': 'effective', 'charge_transport': False,
    'surface_contact': False,
    'transport_model':'constant','activity_model':'volume_fraction',
    'gas_boundary_radius':.001,'vapor_diffusivity':2.5e-5,
    'air_thermal_conductivity':.0262,'air_heat_capacity':1005.,'air_molar_mass':.02897,
    'ambient_pressure':101325.,'flory_huggins_chi':.5,
    'liquid_diffusion':False,'liquid_diffusivity':1e-9,'diffusion_concentration_exponent':3.,
    'noise_sample_interval':1e-5,
}
FIELD_DEFAULTS = {
    'field_mode': 'needle', 'emitter_radius': 0.00035, 'needle_inner_radius': 0.00020,
    'needle_length': 0.01, 'field_radius': 0.15, 'field_nr': 96, 'field_nz': 160,
}
SCENARIO_DEFAULTS = {'branching_enabled':False,'branching_rate_s_inv':100.,'branching_threshold':1.,'branching_max_expected':1000.}
DEFAULTS = {**CORE_DEFAULTS, **FIELD_DEFAULTS, **SCENARIO_DEFAULTS}


def validate(config: dict) -> dict:
    unknown = set(config) - set(DEFAULTS) - {'_metadata'}
    if unknown:
        raise ValueError('Unknown configuration keys: ' + ', '.join(sorted(unknown)))
    c = {**DEFAULTS, **config}
    for name, value in c.items():
        if name == '_metadata':
            if not isinstance(value, dict): raise ValueError('_metadata must be an object')
            continue
        default = DEFAULTS[name]
        if isinstance(default, bool):
            if not isinstance(value, bool): raise ValueError(f'{name} must be boolean')
        elif isinstance(default, str):
            if not isinstance(value, str): raise ValueError(f'{name} must be a string')
        elif isinstance(default, list):
            if not isinstance(value, list) or len(value) != len(default):
                raise ValueError(f'{name} must have {len(default)} components')
            if not all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) for x in value):
                raise ValueError(f'{name} must contain finite numbers')
        else:
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f'{name} must be a finite number')
            if isinstance(default, int) and (value != int(value)):
                raise ValueError(f'{name} must be an integer')
            if isinstance(default, int): c[name] = int(value)
    positive = ['gap','flow_rate','jet_radius','density_polymer','density_solvent','viscosity',
                'relaxation_time','relative_permittivity','temperature','ambient_temperature',
                'latent_heat','molar_mass','psat_temperature','cp_polymer','cp_solvent',
                'air_viscosity','duration','dt_max','dt_min','injection_interval','output_interval',
                'max_beads','max_steps','emitter_radius','needle_inner_radius','needle_length','field_radius',
                'noise_correlation_time','random_seed']
    positive+=['gas_boundary_radius','vapor_diffusivity','air_thermal_conductivity','air_heat_capacity',
               'air_molar_mass','ambient_pressure','liquid_diffusivity','noise_sample_interval']
    for k in positive:
        if c[k] <= 0: raise ValueError(f'{k} must be > 0')
    nonnegative = ['surface_tension','conductivity','mass_transfer','heat_transfer',
                   'psat_ref','air_density','drag_coefficient','coulomb_softening',
                   'perturbation_amplitude','perturbation_frequency']
    for k in nonnegative:
        if c[k] < 0: raise ValueError(f'{k} must be >= 0')
    if c['voltage']*c['current'] < 0: raise ValueError('voltage and convective current must have matching signs')
    if c['perturbation_mode'] not in ('periodic','colored_noise'): raise ValueError('Invalid perturbation_mode')
    if c['constitutive_model'] not in ('maxwell','fene_p'): raise ValueError('Invalid constitutive_model')
    if c['drag_model'] not in ('effective','cylinder'): raise ValueError('Invalid drag_model')
    if c['transport_model'] not in ('constant','gas_film','crossflow'): raise ValueError('Invalid transport_model')
    if c['activity_model'] not in ('volume_fraction','flory_huggins'): raise ValueError('Invalid activity_model')
    if c['flory_huggins_chi']>0.5 or c['flory_huggins_chi']<0:
        raise ValueError('Use 0<=chi<=0.5 for the single-phase high-chain Flory-Huggins closure')
    if c['diffusion_concentration_exponent']<0:raise ValueError('Diffusion concentration exponent must be >=0')
    if c['transport_model']=='gas_film' and c['gas_boundary_radius']<=c['jet_radius']:
        raise ValueError('Gas boundary radius must exceed the initial jet radius')
    if c['transport_model']!='constant' and c['air_density']<=0:raise ValueError('Transport requires positive air density')
    if c['duration']/c['noise_sample_interval']>=1000000:raise ValueError('Noise history too large')
    if c['finite_extensibility'] <= 3: raise ValueError('finite_extensibility (L²) must exceed 3')
    if c['random_seed'] >= 2147483647: raise ValueError('random_seed must be < 2147483647')
    if c['branching_rate_s_inv']<0 or c['branching_threshold']<=0 or c['branching_max_expected']<=1:
        raise ValueError('Branching rate must be >=0, threshold >0 and maximum expected count >1')
    if c['coulomb_softening'] <= 0: raise ValueError('coulomb_softening must be > 0')
    if c['temperature'] <= 150 or c['ambient_temperature'] <= 150 or c['psat_temperature'] <= 150:
        raise ValueError('Absolute temperatures must exceed 150 K')
    if not 0 < c['injection_interval'] < c['duration']: raise ValueError('injection_interval must be < duration')
    inlet_length = c['flow_rate']/(math.pi*c['jet_radius']**2)*c['injection_interval']
    if inlet_length >= .2*c['gap']: raise ValueError('Reduce injection_interval: initial segment exceeds 20% of gap')
    if c['perturbation_amplitude'] >= .1*inlet_length: raise ValueError('Seed RMS must be < 10% of inlet segment length')
    if not 0 < c['polymer_fraction'] < 1: raise ValueError('polymer_fraction must lie strictly between 0 and 1')
    if not 0 <= c['humidity'] <= 1: raise ValueError('humidity must lie between 0 and 1')
    if c['needle_inner_radius'] >= c['emitter_radius']: raise ValueError('Needle inner radius must be smaller than outer radius')
    if c['jet_radius'] >= c['needle_inner_radius']: raise ValueError('Post-meniscus jet radius must be smaller than needle bore radius')
    if c['field_radius'] <= c['emitter_radius']: raise ValueError('Field domain must extend outside the needle')
    if c['needle_length'] <= c['needle_inner_radius']: raise ValueError('Needle length must exceed meniscus radius')
    if c['dt_min'] > c['dt_max']: raise ValueError('dt_min must not exceed dt_max')
    if not 24 <= c['field_nr'] <= 600 or not 48 <= c['field_nz'] <= 1000:
        raise ValueError('Field mesh requires 24<=nr<=600 and 48<=nz<=1000')
    if c['field_mode'] not in ('needle','uniform','file'): raise ValueError('field_mode must be needle, uniform or file')
    if c['field_mode'] == 'file' and not c['field_file']: raise ValueError('field_mode=file requires field_file')
    if c['max_beads'] < math.floor(c['duration']/c['injection_interval']) + 2:
        raise ValueError('max_beads too small for duration/injection_interval (deposited beads are retained)')
    return c


def load_config(path: Path | None = None, overrides: list[str] | None = None) -> dict:
    raw = json.loads(Path(path).read_text(encoding='utf-8')) if path else {}
    if not isinstance(raw, dict): raise ValueError('Configuration must be a JSON object')
    for item in overrides or []:
        key, sep, value = item.partition('=')
        if not sep: raise ValueError('--set requires name=JSON_value')
        try: raw[key] = json.loads(value)
        except json.JSONDecodeError as e: raise ValueError(f'Invalid JSON in --set {item}') from e
    return validate(raw)


def write_config(config: dict, path: Path) -> None:
    path.write_text(json.dumps(config, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def write_namelist(config: dict, path: Path) -> None:
    def literal(x):
        if isinstance(x, bool): return '.true.' if x else '.false.'
        if isinstance(x, str):
            if '\n' in x or '\r' in x: raise ValueError('Newlines forbidden in namelist paths')
            return "'" + x.replace("'", "''") + "'"
        if isinstance(x, list): return ', '.join(literal(v) for v in x)
        return str(x)
    path.write_text('&jet3d\n' + '\n'.join(f' {k} = {literal(config[k])}' for k in CORE_DEFAULTS) + '\n/\n', encoding='utf-8')
