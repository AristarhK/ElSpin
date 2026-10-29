"""Material records separate measured, tabulated and assumed properties."""
from __future__ import annotations
from copy import deepcopy
import json
from pathlib import Path
from .config import DEFAULTS, validate

MATERIAL_KEYS = [
    'polymer_fraction','density_polymer','density_solvent','viscosity','relaxation_time',
    'surface_tension','conductivity','relative_permittivity','temperature','molar_mass',
    'latent_heat','psat_ref','psat_temperature','cp_polymer','cp_solvent',
    'viscosity_exponent','activation_energy','mass_transfer','heat_transfer','finite_extensibility',
]
LEGACY_MATERIAL_KEYS=MATERIAL_KEYS.copy()
MATERIAL_KEYS+=['transport_model','activity_model','vapor_diffusivity','flory_huggins_chi',
                'liquid_diffusion','liquid_diffusivity','diffusion_concentration_exponent']


def pva_water() -> dict:
    c=deepcopy(DEFAULTS)
    c.update(polymer_fraction=.10,density_polymer=1260.,
             perturbation_mode='colored_noise',constitutive_model='fene_p',drag_model='cylinder',
             charge_transport=True,surface_contact=True,duration=.04,output_interval=.0004,max_beads=500,
             transport_model='gas_film',activity_model='flory_huggins',liquid_diffusion=True)
    c['_metadata']={
        'material_name':'ПВС–вода, 10 масс.% — стартовая гипотеза',
        'polymer':'ПВС (PVA)', 'solvent':'вода', 'filler':'нет',
        'molecular_weight_g_mol':None, 'hydrolysis_degree_percent':None,
        'property_status':'Не калибровано: реология, K, χ и жидкая диффузия — гипотезы. Газовый перенос рассчитывается по геометрии.',
        'source':'Свойства воды при 25 °C — справочные приближения; раствор — иллюстративный набор на основе пользовательского архива.',
        'note':'Измерить η, λ, K, γ, Mw и гидролиз. R газа=1 мм, D жидкости=1e-9 м²/с и χ=0.5 — гипотезы для анализа диапазонов. Начальная струя задана после неразрешённого мениска.',
    }
    return validate(c)


def material_presets() -> dict[str,dict]:
    """Ready-to-edit starting hypotheses, NEVER a measured solution database."""
    records={}
    for w,eta,relax in [(.08,.08,.004),(.10,.20,.010),(.12,.50,.020)]:
        c=pva_water();c.update(polymer_fraction=w,viscosity=eta,relaxation_time=relax)
        name=f'ПВС–вода · {w*100:.0f} масс.%'
        c['_metadata']['material_name']=name+' · предварительный набор'
        c['_metadata']['source']='Свойства воды при 25 °C — справочные приближения. η и λ выбраны как гипотезы для сравнения концентраций, не как измеренная зависимость.'
        records[name]=validate(c)
    c=pva_water();c.update(polymer_fraction=.04,density_polymer=1210.,viscosity=.25,relaxation_time=.025)
    c['_metadata'].update(material_name='ПЭО–вода · 4 масс.% · предварительный набор',polymer='ПЭО (PEO)',
        source='Водный раствор ПЭО: параметры реологии иллюстративные; молекулярная масса и измерения требуются.')
    records['ПЭО–вода · 4 масс.%']=validate(c)
    c=pva_water();c.update(polymer_fraction=.12,density_polymer=1200.,density_solvent=789.,
        viscosity=.12,relaxation_time=.008,surface_tension=.024,relative_permittivity=25.,vapor_diffusivity=1.2e-5,
        molar_mass=.046068,latent_heat=846000.,psat_ref=7870.,cp_solvent=2440.,humidity=.05,
        activation_energy=15000.)
    c['_metadata'].update(material_name='ПВП–этанол · 12 масс.% · предварительный набор',polymer='ПВП (PVP)',solvent='этанол',
        source='Свойства этанола при 25 °C — справочные приближения. Свойства раствора и перенос — гипотезы; η, λ, K, γ измерить.',
        note='Насыщенность воздуха относится к этанолу, а не к относительной влажности воды. Однокомпонентное испарение; смесь растворителей не моделируется.')
    records['ПВП–этанол · 12 масс.%']=validate(c)
    return records


def save_material(c: dict, name: str, source: str, note: str, path: Path) -> None:
    c=validate(c)
    if not name.strip(): raise ValueError('Укажите название материала')
    record={'schema_version':2,'name':name.strip(),'source':source.strip(),'note':note.strip(),
            'properties':{k:c[k] for k in MATERIAL_KEYS}}
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists(): raise FileExistsError('Файл материала уже существует. Выберите новое имя.')
    path.write_text(json.dumps(record,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def load_material(path: Path, base: dict) -> dict:
    record=json.loads(Path(path).read_text(encoding='utf-8'))
    if record.get('schema_version') not in (1,2): raise ValueError('Неизвестная версия материала')
    props=record['properties']
    expected=LEGACY_MATERIAL_KEYS if record['schema_version']==1 else MATERIAL_KEYS
    if set(props)!=set(expected): raise ValueError('Неполный или неизвестный набор свойств материала')
    c=deepcopy(base); c.update(props)
    if record['schema_version']==1:c.update(transport_model='constant',activity_model='volume_fraction',liquid_diffusion=False)
    c['_metadata']={'material_name':record['name'],'source':record.get('source',''),
                    'note':record.get('note',''),'property_status':'Пользовательский ввод; проверить происхождение свойств.'}
    return validate(c)
