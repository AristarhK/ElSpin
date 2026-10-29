"""Experimental dependencies for ONE fixed operating condition.

Filaments: a phenomenological saturation of N**(3/2), motivated by the
equal-split charge/capillary scaling. No verified jet birth law is implied.
Drying: bounded monotone saturation. Neither fit replaces the 3D equations.
Coordinate-group CV prevents leakage between replicates at identical z.
Uncertainty bands are conditional bootstrap bands, not guaranteed accuracy.
"""
from __future__ import annotations
import csv,io,json,math
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.interpolate import PchipInterpolator
from . import __version__

TARGETS={'filaments':'Дочерние струи: среднее число N',
         'evaporated_percent':'Испарилось исходного растворителя, %',
         'remaining_percent':'Осталось исходного растворителя, %'}


def validate_data(z,y,target,sigma=None):
    if target not in TARGETS:raise ValueError('Неизвестная измеряемая величина')
    z=np.asarray(z,float);y=np.asarray(y,float)
    if z.ndim!=1 or y.shape!=z.shape or not len(z):raise ValueError('Нужны пары координата–измерение')
    if not np.isfinite(z).all() or not np.isfinite(y).all() or np.any(z<0):raise ValueError('Координаты и значения должны быть конечными, z>=0')
    if target=='filaments' and np.any(y<1):raise ValueError('Для одной подающей струи N должно быть >=1. Среднее число допускается; отсутствие струи — другой режим.')
    if target!='filaments' and np.any((y<0)|(y>100)):raise ValueError('Доля растворителя должна быть в диапазоне 0–100%')
    if sigma is not None:
        sigma=np.asarray(sigma,float)
        if sigma.shape!=z.shape or not np.isfinite(sigma).all() or np.any(sigma<=0):raise ValueError('σ должна быть положительной и заданной для всех строк')
    return z,y,sigma


def parse_table(text):
    """Excel TSV, semicolon decimal-comma CSV, or ordinary decimal-point CSV."""
    text=text.lstrip('\ufeff');lines=[s for s in text.splitlines() if s.strip()]
    if not lines:raise ValueError('Таблица пуста')
    sample='\n'.join(lines[:15])
    delimiter='\t' if '\t' in sample else ';' if ';' in sample else ','
    rows=list(csv.reader(io.StringIO('\n'.join(lines)),delimiter=delimiter))
    def number(s):return float(s.strip().replace(',','.'))
    factor=1.;has_header=False
    try:number(rows[0][0]);number(rows[0][1])
    except (ValueError,IndexError):
        header=[s.strip().lower() for s in rows.pop(0)]
        if not header or header[0].replace(' ','') not in ('z','z_mm','z_m','distance_mm','distance_m','координата','координатаz,мм','z,мм'):
            raise ValueError('Первая строка должна содержать числа или заголовок z_mm / z_m')
        if header and header[0] in ('z_m','distance_m'):factor=1000.
        has_header=True
    z=[];y=[];sigma=[]
    for i,row in enumerate(rows,2 if has_header else 1):
        if len(row) not in (2,3):raise ValueError(f'Строка {i}: нужны 2 или 3 столбца. Для десятичной запятой используйте ; или табуляцию.')
        try:z.append(number(row[0])*factor);y.append(number(row[1]));sigma.append(number(row[2]) if len(row)==3 and row[2].strip() else None)
        except ValueError as exc:raise ValueError(f'Строка {i}: проверьте числовые значения') from exc
    if any(v is None for v in sigma) and not all(v is None for v in sigma):raise ValueError('Укажите σ для всех строк или оставьте весь столбец пустым')
    return np.asarray(z),np.asarray(y),(None if all(v is None for v in sigma) else np.asarray(sigma))


def evaluate(model,z):
    z=np.asarray(z,float)
    if model['method']=='pchip':
        return PchipInterpolator(model['knots_z_mm'],model['knots_y'],extrapolate=False)(z)
    x=np.maximum(z-model['onset_mm'],0)/model['length_mm']
    fraction=-np.expm1(-x)
    if model['target']=='filaments':return (1+(model['plateau']**1.5-1)*fraction)**(2/3)
    lost=model['plateau']*fraction
    return 100-lost if model['target']=='remaining_percent' else lost


def _fit(z,y,target,method,onset,fit_onset,sigma,cancelled=None):
    if method=='pchip':
        knots=np.unique(z);means=[float(np.average(y[z==x],weights=None if sigma is None else 1/sigma[z==x]**2)) for x in knots]
        return {'method':'pchip','target':target,'knots_z_mm':knots.tolist(),'knots_y':means},True,None
    span=max(float(np.ptp(z)),1.);weights=np.ones(len(y)) if sigma is None else sigma
    boundary=1. if target=='filaments' else 100. if target=='remaining_percent' else 0.
    if np.all(np.abs(y-boundary)<1e-12):
        return {'method':'saturation','target':target,'plateau':1. if target=='filaments' else 0.,
                'length_mm':span,'onset_mm':float(onset)},False,None
    cap=max(1000.,100*float(y.max())) if target=='filaments' else 100.
    lower=1. if target=='filaments' else 0.
    signal=y if target!='remaining_percent' else 100-y
    p0=min(cap*.99,max(lower+.01,float(signal.max())*1.3))
    best=None
    lo=[lower,math.log(span*1e-4)];hi=[cap,math.log(span*1e4)]
    if fit_onset:lo.append(0.);hi.append(float(z.max())-1e-9)
    def unpack(p):return {'method':'saturation','target':target,'plateau':float(p[0]),'length_mm':float(np.exp(p[1])),
                          'onset_mm':float(p[2] if fit_onset else onset)}
    def residual(p):
        if cancelled and cancelled():raise InterruptedError('Подбор остановлен; исходные данные сохранены в окне')
        return (evaluate(unpack(p),z)-y)/weights
    for scale in [.2,1.,5.]:
        initial=[p0,math.log(span*scale)]
        if fit_onset:initial.append(min(float(z.max())*.15,max(0.,onset)))
        fit=least_squares(residual,initial,bounds=(lo,hi),max_nfev=1500)
        if fit.success and np.isfinite(fit.cost) and (best is None or fit.cost<best.cost):best=fit
    if best is None:raise RuntimeError('Подбор не сошёлся; проверьте данные и выбранную зависимость')
    norms=np.linalg.norm(best.jac,axis=0)
    normalized=best.jac/np.maximum(norms,1e-20)
    condition=float(np.linalg.cond(normalized))
    identified=bool(np.all(norms>1e-8) and condition<1e4 and not np.any(best.active_mask) and np.ptp(y)>1e-8)
    return unpack(best.x),identified,condition


def predict(result,z,*,allow_extrapolation=False):
    z=np.atleast_1d(np.asarray(z,float))
    if not np.isfinite(z).all() or np.any(z<0):raise ValueError('Координата должна быть конечной и >=0')
    lo,hi=result['domain_mm'];outside=np.any((z<lo)|(z>hi))
    if outside and not allow_extrapolation:raise ValueError(f'Координата вне измеренного диапазона {lo:g}–{hi:g} мм. Экстраполяция выключена.')
    if outside and (result['method']=='pchip' or not result['parameters_identified']):
        raise ValueError('Экстраполяция недоступна: интерполяция или неидентифицируемые параметры. Нужны дополнительные измерения.')
    return evaluate(result['model'],z)


def fit_dependency(z,y,target='filaments',*,method='saturation',sigma=None,onset_mm=0.,fit_onset=False,
                   bootstrap=200,seed=42,cancelled=None,progress=None):
    z,y,sigma=validate_data(z,y,target,sigma)
    if method not in ('saturation','pchip'):raise ValueError('Неизвестный метод')
    if not math.isfinite(onset_mm) or onset_mm<0:raise ValueError('Начало зависимости должно быть >=0')
    unique=np.unique(z);minimum=2 if method=='pchip' else 6 if fit_onset else 5
    if len(unique)<minimum:raise ValueError(f'Нужно не менее {minimum} различных координат; сейчас {len(unique)}')
    if method!='pchip' and not fit_onset and onset_mm>=z.max():raise ValueError('Начало зависимости должно быть раньше последнего измерения')
    if not isinstance(bootstrap,int) or not 0<=bootstrap<=2000:raise ValueError('Число bootstrap-повторов должно быть от 0 до 2000')
    def check():
        if cancelled and cancelled():raise InterruptedError('Подбор остановлен; исходные данные сохранены в окне')
    check();model,identified,condition=_fit(z,y,target,method,onset_mm,fit_onset,sigma,cancelled)
    fitted=evaluate(model,z);residual=y-fitted
    grid=np.linspace(float(z.min()),float(z.max()),301);mean=evaluate(model,grid)
    warnings=[]
    if method=='saturation':warnings.append('Модель насыщения — эмпирическое замыкание для этой серии, а не подтверждённый закон рождения струй или высыхания.')
    if len(unique)<8:warnings.append('Меньше 8 различных координат: оценка ошибки чувствительна к отдельным точкам. Желательны дополнительные координаты и повторные опыты.')
    if not identified:warnings.append('Параметры недостаточно определены данными; экстраполяция заблокирована. Нужны точки начала роста и насыщения.')
    if np.ptp(y)<1e-8:warnings.append('Все значения одинаковы: изменение с координатой не обнаружено; поведение вне измерений неизвестно.')
    means=np.array([y[z==x].mean() for x in unique]);direction=-1 if target=='remaining_percent' else 1
    if method=='saturation' and np.any(direction*np.diff(means)<-.05*max(float(np.ptp(y)),1.)):
        warnings.append('Данные содержат заметное изменение против монотонной модели. Проверьте измерения, слияние/потерю нитей или выберите интерполяцию.')
    # Hold out whole coordinates: replicates at the test coordinate never train the model.
    folds=unique[np.unique(np.linspace(0,len(unique)-1,min(25,len(unique))).astype(int))]
    errors=[];failed=[];used=0;held_out=[]
    for x in folds:
        check();train=z!=x
        if method=='pchip' and x in (unique[0],unique[-1]):
            failed.append({'z_mm':float(x),'reason':'Интерполяция: крайняя координата вне обучающего диапазона'});continue
        try:
            m,_,_=_fit(z[train],y[train],target,method,onset_mm,fit_onset,None if sigma is None else sigma[train],cancelled)
            prediction=evaluate(m,z[~train]);errors.extend((y[~train]-prediction).tolist());used+=1
            held_out.extend({'z_mm':float(a),'observed':float(b),'predicted':float(p),'error':float(b-p)} for a,b,p in zip(z[~train],y[~train],prediction))
        except (ValueError,RuntimeError) as exc:failed.append({'z_mm':float(x),'reason':str(exc)})
        if progress:progress(int(25*(len(failed)+used)/len(folds)))
    cv={'method':'leave-one-coordinate-out','coordinate_groups_planned':len(folds),'coordinate_groups_used':used,
        'observations_predicted':len(errors),'rmse':float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
        'mae':float(np.mean(np.abs(errors))) if errors else None,'failed_or_unavailable':failed,'held_out_predictions':held_out}
    if used<len(folds):warnings.append(f'Проверено {used} из {len(folds)} исключённых координат; причины пропусков сохранены.')
    # Conditional mean-curve bands. No automatic interpretation as a physical prediction interval.
    curves=[];params=[];rng=np.random.default_rng(seed);failures=0
    effective_bootstrap=bootstrap
    if method=='pchip':
        effective_bootstrap=0
        warnings.append('Интерполяция: полоса неопределённости не строится. Для оценки параметров выберите модель с насыщением.')
    elif sigma is None and np.mean(residual**2)<1e-16*max(1.,float(np.mean(y**2))):
        effective_bootstrap=0
        warnings.append('Остатки практически нулевые, а σ не задана. Физическая неопределённость неизвестна; нулевая полоса не показывается.')
    for i in range(effective_bootstrap):
        check()
        if sigma is not None:yb=fitted+rng.normal(0,sigma)
        else:yb=fitted+residual*rng.choice([-1.,1.],len(y))
        try:
            m,_,_=_fit(z,yb,target,method,onset_mm,fit_onset,sigma,cancelled);curve=evaluate(m,grid)
            if not np.isfinite(curve).all():raise ValueError('Нечисловой bootstrap-прогноз')
            curves.append(curve)
            if method=='saturation':params.append([m['plateau'],m['length_mm'],m['onset_mm']])
        except (ValueError,RuntimeError):failures+=1
        if progress:progress(25+int(75*(i+1)/max(effective_bootstrap,1)))
    band=None;intervals=None
    if len(curves)>=20:
        quant=np.percentile(curves,[2.5,97.5],axis=0);band={'lower':quant[0].tolist(),'upper':quant[1].tolist()}
        if params:
            q=np.percentile(params,[2.5,97.5],axis=0)
            intervals={k:[float(q[0,i]),float(q[1,i])] for i,k in enumerate(['plateau','length_mm','onset_mm'])}
    elif effective_bootstrap:warnings.append('Успешных bootstrap-повторов меньше 20; полоса неопределённости не показана.')
    warnings.append('Зависимость относится только к условиям этой серии. Проверка исключённых точек не подтверждает перенос на другой раствор, U, Q или влажность.')
    warnings.append('Полоса — неопределённость средней кривой выбранной модели; не гарантированная ошибка отдельного опыта и не оценка всей физической погрешности.')
    variance=float(np.sum((y-y.mean())**2))
    if condition is not None and not np.isfinite(condition):condition=None
    return {'schema_version':1,'app_version':__version__,'created_at_utc':datetime.now(timezone.utc).isoformat(),
            'target':target,'method':method,'model':model,'domain_mm':[float(z.min()),float(z.max())],
            'parameters_identified':identified,'normalized_jacobian_condition':condition,
            'z_mm':z.tolist(),'observed':y.tolist(),'sigma':None if sigma is None else sigma.tolist(),
            'fitted_at_observations':fitted.tolist(),'residuals':residual.tolist(),
            'training':{'rmse':float(np.sqrt(np.mean(residual**2))),'mae':float(np.mean(np.abs(residual))),
                        'r_squared':float(1-np.sum(residual**2)/variance) if variance>0 else None},
            'validation':cv,'grid_z_mm':grid.tolist(),'mean':mean.tolist(),'mean_band':band,
            'parameter_intervals':intervals,'bootstrap_requested':bootstrap,'bootstrap_successful':len(curves),
            'bootstrap_failed':failures,'bootstrap_seed':seed,
            'band_assumption':'Gaussian observation noise from supplied sigma' if sigma is not None else 'Wild residual bootstrap; assumes residuals represent measurement/model noise',
            'experimental_accuracy_certified':False,'warnings':warnings}


def save_experiment(path,z,y,target,*,sigma=None,conditions=None,label='',result=None,settings=None):
    z,y,sigma=validate_data(z,y,target,sigma);path=Path(path)
    if path.exists():raise FileExistsError('Файл уже существует; выберите новое имя')
    payload={'schema_version':1,'kind':'elspin_experimental_series','label':label,'target':target,
             'z_mm':z.tolist(),'observed':y.tolist(),'sigma':None if sigma is None else sigma.tolist(),
             'conditions':conditions or {},'settings':settings or {},'fit':result,
             'data_origin':'User-entered; source and experimental status not independently verified'}
    path.write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def load_experiment(path):
    raw=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if raw.get('schema_version')!=1 or raw.get('kind')!='elspin_experimental_series':raise ValueError('Неизвестный формат серии')
    validate_data(raw['z_mm'],raw['observed'],raw['target'],raw.get('sigma'))
    if not isinstance(raw.get('conditions',{}),dict):raise ValueError('Неверный формат условий')
    if not isinstance(raw.get('settings',{}),dict):raise ValueError('Неверный формат настроек')
    # Stored numerical fits are review records; recompute from the loaded observations.
    raw['fit']=None
    return raw
