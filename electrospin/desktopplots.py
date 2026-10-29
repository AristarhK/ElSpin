"""Native Matplotlib plots. All coordinates and maps come from solver outputs."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from .analysis import read_table, coverage_map, deposited_segments, field_data, field_grid

PRIMARY_PLOTS=['Падение на коллектор — 2D','Высыхание по координате z','Филаменты: устойчивость и сценарии']
PLOT_DESCRIPTIONS={
    PRIMARY_PLOTS[0]:'Слева — полёт струи к коллектору, справа — осевшая нить при взгляде сверху. Ползунок меняет время.',
    PRIMARY_PLOTS[1]:'Доля испарившегося растворителя вдоль оси z от иглы. Показана история первого материального сегмента до выбранного времени.',
    PRIMARY_PLOTS[2]:'По умолчанию показаны условия неустойчивости, а не выдуманный прогноз N(z). Число дочерних струй требует отдельного закона рождения и данных; экспертный сценарий явно отделён.',
}
LEGACY_PLOTS=['Струя 3D','Проекции 2D','Диаметр и состав','Реология и температура',
            'Заряд и скорость','Коллектор','Рост изгиба','Электрическое поле','Баланс и шаг',
            'Высыхание','Филаменты и ветвление']
PLOT_NAMES=PRIMARY_PLOTS+LEGACY_PLOTS


class RunData:
    def __init__(self,path):
        self.path=Path(path)
        from .config import validate
        self.config=validate(json.loads((self.path/'config.json').read_text(encoding='utf-8')))
        self.summary=json.loads((self.path/'analysis.json').read_text(encoding='utf-8'))
        self.trajectory=read_table(self.path/'trajectory.csv')
        self.deposition=read_table(self.path/'deposition.csv')
        self.diagnostics=read_table(self.path/'diagnostics.csv')
        self.instability=read_table(self.path/'instability.csv')
        self.bending=read_table(self.path/'bending_profile.csv')
        self.branching=read_table(self.path/'branching_scenario.csv')
        self.branch_screen=read_table(self.path/'branching_screen.csv')
        self.field=field_data(self.path)
        self.times=np.unique(self.trajectory['time_s'])

    def frame(self,index=-1):
        if not len(self.times): return self.trajectory
        a=self.trajectory[self.trajectory['time_s']==self.times[index]]
        return np.sort(a,order='bead_id')


def empty(ax,text):
    ax.text(.5,.5,text,ha='center',va='center',transform=ax.transAxes,wrap=True)


def draw(fig: Figure, kind: str, run: RunData, index=-1, true_scale=True):
    fig.clear(); c=run.config; a=run.frame(index)
    live=a[a['deposited']<.5]
    time=(run.times[index]*1000 if len(run.times) else 0)
    if kind==PRIMARY_PLOTS[0]:
        # Two explicit views of the same process, rather than a report dashboard.
        ax=fig.add_subplot(121)
        if len(a)>1:
            connected=np.diff(a['bead_id'])==1
            points=np.column_stack([a['x_m'],a['z_m']])*1000
            seg=np.stack([points[:-1][connected],points[1:][connected]],axis=1)
            ax.add_collection(LineCollection(seg,colors='#146c86',linewidths=1.8))
        if len(a):ax.scatter(a['x_m']*1000,a['z_m']*1000,c='#146c86',s=3)
        ax.scatter([0],[0],c='#bd742c',s=40,zorder=5,label='Игла / начальное сечение')
        ax.axhline(c['gap']*1000,c='#404b50',lw=3,label='Плоскость коллектора')
        span=max(float(np.max(np.abs(a['x_m'])))*1000 if len(a) else 0,.1)*1.2
        ax.set(xlim=(-span,span),ylim=(c['gap']*1030,-c['gap']*30),xlabel='Поперечная координата x, мм',ylabel='Расстояние от иглы z, мм',title='Вид сбоку: падение струи (x–z)')
        ax.grid(alpha=.2);ax.legend(loc='upper right',fontsize=8)
        ax=fig.add_subplot(122)
        dep=run.deposition[run.deposition['time_s']<=time*.001+1e-12]
        segments,_,_=deposited_segments(dep)
        if len(segments):ax.add_collection(LineCollection(segments*1000,colors='#167e6b',linewidths=1.5));ax.autoscale()
        if len(dep):ax.scatter(dep['collector_x_m']*1000,dep['collector_y_m']*1000,c=dep['time_s']*1000,cmap='viridis',s=12)
        else:empty(ax,'К этому времени нить ещё не достигла коллектора.\nПереместите ползунок времени вправо.')
        ax.set(xlabel='x на коллекторе, мм',ylabel='y на коллекторе, мм',title='Вид сверху: нити на коллекторе (x–y)',aspect='equal');ax.grid(alpha=.2)
        fig.suptitle(f'Падение и укладка нити · t = {time:.2f} мс',fontsize=13)
    elif kind==PRIMARY_PLOTS[1]:
        from .drying import parcel_history,drying_values
        history,_=parcel_history(run.trajectory[run.trajectory['time_s']<=time*.001+1e-12])
        ax=fig.add_subplot(111)
        if len(history):
            retention,_=drying_values(history,c);evaporated=(1-retention)*100
            ax.plot(history['z_m']*1000,evaporated,color='#176eaa',lw=2.3,label='Испарившийся растворитель')
            ax.fill_between(history['z_m']*1000,evaporated,0,color='#176eaa',alpha=.1)
            ax.set_ylim(0,max(float(np.max(evaporated))*1.3,.01))
            ax.annotate(f'В последней рассчитанной точке:\nиспарилось {evaporated[-1]:.3f}% · осталось {retention[-1]*100:.3f}%',xy=(.04,.94),xycoords='axes fraction',va='top',bbox=dict(boxstyle='round,pad=.6',fc='#f0f6fa',ec='#a8bdcc'))
        else:empty(ax,'Нет материальных сегментов в выбранном кадре.')
        ax.set(xlim=(0,c['gap']*1000),xlabel='Расстояние от иглы z, мм',ylabel='Испарилось исходного растворителя, %',title='Высыхание раствора по координате z')
        ax.grid(alpha=.25)
    elif kind==PRIMARY_PLOTS[2]:
        from .branching import axial_profile
        ax=fig.add_subplot(111);b=run.branching
        screen=getattr(run,'branch_screen',[])
        if not c['branching_enabled'] and len(screen):
            ax.plot(screen['z_m']*1000,screen['electric_capillary_pressure_ratio'],lw=2,color='#78629d',label='Оценка электрического / капиллярного давления')
            ax.axhline(1.5,color='#a15c24',ls='--',label='Нейтральная граница m=2 в пределе проводящего цилиндра')
            ax.set(xlim=(0,c['gap']*1000),xlabel='Расстояние от иглы z, мм',ylabel='Электрическое / капиллярное давление',title='Условия деформации сечения — это не прогноз числа филаментов')
            ax.legend(loc='upper left',fontsize=8);ax.grid(alpha=.25)
            fraction=float(np.mean(screen['limit_model_in_domain']>=.5))
            fig.text(.5,.035,'N(z) не определено: закон рождения дочерних струй не идентифицирован. Критерий не доказывает ветвление.',ha='center',fontsize=9)
            fig.text(.5,.01,f'Условия применимости локальной предельной модели выполнены в {fraction*100:.1f}% срезов. См. branching_screen.csv.',ha='center',fontsize=9)
            fig.tight_layout(rect=(0,.07,1,1));return
        z,count=axial_profile(run.trajectory,b)
        if len(z):
            ax.plot(z*1000,count,color='#8a4a91',lw=2.3,label='Дочерние филаменты: сценарная оценка')
            upper=max(float(count.max())*1.4,2.);ax.set_ylim(0,upper)
            enabled=c['branching_enabled']
            text=('Сценарий ветвления выключен: N(z) = 1.\nДля сценарного расчёта включите «Ветвление» слева.' if not enabled else
                  f'Сценарий включён · β = {c["branching_rate_s_inv"]:g} 1/с\nN — ожидаемое число, закон ветвления требует калибровки.')
            if enabled and np.ptp(count)<1e-8:text+='\nВ этом расчёте увеличения числа филаментов нет.'
            ax.annotate(text,xy=(.04,.94),xycoords='axes fraction',va='top',bbox=dict(boxstyle='round,pad=.6',fc='#fcf6ed',ec='#ccac77'))
        else:empty(ax,'Прогноз N(z) отсутствует.\nДля реального числа дочерних струй нужны закон рождения и его проверка.\nНовый расчёт сохраняет физический скрининг устойчивости.\nВкладка «Ветвление» содержит только экспертный сценарий.')
        ax.set(xlim=(0,c['gap']*1000),xlabel='Расстояние от иглы z, мм',ylabel='Ожидаемое число филаментов N',title='Число дочерних филаментов в зависимости от z')
        ax.grid(alpha=.25)
    elif kind=='Струя 3D':
        ax=fig.add_subplot(111,projection='3d')
        if len(live)>1:
            pts=np.column_stack([live[k] for k in ('x_m','y_m','z_m')])*1000
            links=np.diff(live['bead_id'])==1
            lc=Line3DCollection(np.stack([pts[:-1][links],pts[1:][links]],axis=1),cmap='plasma',linewidths=1.6)
            lc.set_array(2*live['radius_m'][:-1][links]*1e6);ax.add_collection3d(lc)
            fig.colorbar(lc,ax=ax,shrink=.55,pad=.09,label='Эквивалентный диаметр, мкм')
            width=max(np.ptp(pts[:,0]),np.ptp(pts[:,1]),.02)*1.2
            center=pts[:,:2].mean(axis=0)
            ax.set_xlim(center[0]-width/2,center[0]+width/2);ax.set_ylim(center[1]-width/2,center[1]+width/2)
            if true_scale: ax.set_box_aspect([width,width,c['gap']*1000])
            else: ax.set_box_aspect([1,1,2])
        ax.scatter([0],[0],[0],c='#cf7a21',s=35,label='Начальное сечение')
        ax.set(zlabel='z, мм',xlabel='x, мм',ylabel='y, мм',zlim=(c['gap']*1000,0))
        ax.set_title(f't = {time:.2f} мс · '+('единый масштаб осей' if true_scale else 'поперечный масштаб увеличен'))
    elif kind=='Проекции 2D':
        for n,(x,y) in enumerate([('x_m','z_m'),('y_m','z_m'),('x_m','y_m')],1):
            ax=fig.add_subplot(1,3,n);ax.plot(a[x]*1e3,a[y]*1e3,lw=1)
            ax.set(xlabel=x[0]+', мм',ylabel=y[0]+', мм',title=f'{x[0]}–{y[0]}')
            if y=='z_m':ax.invert_yaxis()
            ax.grid(alpha=.2)
    elif kind in ('Диаметр и состав','Реология и температура','Заряд и скорость'):
        if not len(live):empty(fig.add_subplot(111),'В выбранном кадре нет свободной струи');return
        order=np.argsort(live['z_m']); s=live[order]; z=s['z_m']*1000
        if kind=='Диаметр и состав':
            values=[(2*s['radius_m']*1e6,'Диаметр, мкм'),(s['polymer_mass_kg']/s['mass_kg']*100,'Полимер, масс.%')]
        elif kind=='Реология и температура':
            values=[(s['stress_Pa']/1000,'Продольное напряжение, кПа'),(s['temperature_K']-273.15,'Температура, °C')]
        else:
            values=[(s['charge_C']/s['segment_length_m']*1e9,'Линейный заряд, нКл/м'),
                    (np.sqrt(s['vx_m_s']**2+s['vy_m_s']**2+s['vz_m_s']**2),'Скорость, м/с')]
        for j,(v,label) in enumerate(values,1):
            ax=fig.add_subplot(1,2,j);ax.plot(z,v,'.-',ms=3);ax.set(xlabel='z, мм',ylabel=label);ax.grid(alpha=.25)
        fig.suptitle('Значения сегментов; сортировка по z не задаёт связность струи')
    elif kind=='Коллектор':
        dep=run.deposition; seg,mass,length=deposited_segments(dep)
        ax=fig.add_subplot(221);ax.set(title='Укладка в координатах коллектора',xlabel='x, мм',ylabel='y, мм',aspect='equal')
        if len(seg): ax.add_collection(LineCollection(seg*1e3,colors='#157d82',linewidths=.8));ax.autoscale()
        else:empty(ax,'Связанных осевших сегментов пока нет')
        cov=coverage_map(dep);ax=fig.add_subplot(222)
        image=ax.pcolormesh(cov['xedges']*1000,cov['yedges']*1000,cov['areal_polymer_mass_kg_m2']*1000,cmap='viridis')
        fig.colorbar(image,ax=ax,label='Масса полимера, г/м²');ax.set(title='Консервативная карта массы',xlabel='x, мм',ylabel='y, мм',aspect='equal')
        ax=fig.add_subplot(223)
        if len(dep): ax.hist(dep['dry_radius_m']*2e9,bins=min(30,max(5,len(dep)//5)),color='#be7650')
        else:empty(ax,'Нет попаданий — диаметр осевшего волокна не определён')
        ax.set(xlabel='Сухой диаметр при неизменной длине, нм',ylabel='Число попаданий')
        ax=fig.add_subplot(224,projection='3d')
        x=(cov['xedges'][:-1]+cov['xedges'][1:])*500;y=(cov['yedges'][:-1]+cov['yedges'][1:])*500
        xx,yy=np.meshgrid(x,y);ax.plot_surface(xx,yy,cov['areal_polymer_mass_kg_m2']*1000,cmap='viridis',linewidth=0)
        ax.set(xlabel='x, мм',ylabel='y, мм',zlabel='г/м²',title='Поверхность массы, не толщина ковра')
    elif kind=='Рост изгиба':
        h=run.instability; b=run.bending
        ax=fig.add_subplot(221)
        if len(h):ax.semilogy(h['time_s']*1000,np.maximum(h['rms_transverse_m']*1e6,1e-10),label='RMS');ax.semilogy(h['time_s']*1000,np.maximum(h['centered_rms_m']*1e6,1e-10),label='без сдвига центра');ax.legend()
        ax.set(xlabel='t, мс',ylabel='Поперечный размер, мкм',title='Рост возмущений')
        ax=fig.add_subplot(222)
        if len(b):ax.semilogy(b['z_m']*1000,np.maximum(b['rms_transverse_m']*1e6,1e-10),'.-')
        ax.set(xlabel='z, мм',ylabel='RMS, мкм',title='Статистика по снимкам')
        dep=run.deposition;ax=fig.add_subplot(223)
        if len(dep)>8:
            dep=np.sort(dep,order='time_s');ts=np.linspace(dep['time_s'].min(),dep['time_s'].max(),len(dep))
            x=np.interp(ts,dep['time_s'],dep['x_m']);y=np.interp(ts,dep['time_s'],dep['y_m'])
            window=np.hanning(len(ts));power=(abs(np.fft.rfft((x-x.mean())*window))**2+abs(np.fft.rfft((y-y.mean())*window))**2)
            freq=np.fft.rfftfreq(len(ts),ts[1]-ts[0]);ax.plot(freq[1:],power[1:]);ax.set(xlabel='Частота, Гц',ylabel='Мощность, усл. ед.',title='Спектр координат попаданий (Hann)')
        else:empty(ax,'Для спектра нужно более 8 попаданий')
        ax=fig.add_subplot(224);ax.axis('off');empty(ax,'RMS и спектр — диагностика.\nДля вывода о неустойчивости сравните\nсилы Кулона вкл./выкл., несколько seed,\nшаг времени и длину сегмента.\nТраектория не задаётся спиралью.')
    elif kind=='Электрическое поле':
        if not len(run.field):empty(fig.add_subplot(111),'Выбрано однородное поле или нет файла поля');return
        r,z,phi=field_grid(run.field,'Phi_V');_,_,er=field_grid(run.field,'E_r_V_m');_,_,ez=field_grid(run.field,'E_z_V_m')
        ax=fig.add_subplot(121);v=ax.pcolormesh(r*1000,z*1000,np.hypot(er,ez)/1000,shading='auto',cmap='magma')
        ax.set(xlabel='r, мм',ylabel='z, мм',title='Модуль поля, кВ/м',xlim=(0,min(c['field_radius']*1000,15)));ax.invert_yaxis();fig.colorbar(v,ax=ax)
        ax=fig.add_subplot(122,projection='3d');rr,zz=np.meshgrid(r*1000,z*1000)
        ax.plot_surface(rr,zz,phi/1000,cmap='viridis',rcount=60,ccount=50)
        ax.set(xlabel='r, мм',ylabel='z, мм',zlabel='φ, кВ',title='Осесимметричный потенциал')
    elif kind=='Баланс и шаг':
        d=run.diagnostics; ax=fig.add_subplot(121)
        for name,label in [('mass_relative_error','масса'),('polymer_relative_error','полимер'),('charge_relative_error','заряд')]:
            ax.semilogy(d['time_s']*1000,np.maximum(d[name],1e-18),label=label)
        ax.legend();ax.set(xlabel='t, мс',ylabel='Относительная ошибка',title='Балансы')
        ax=fig.add_subplot(122);ax.semilogy(d['time_s']*1000,np.maximum(d['dt_s'],1e-18));ax.set(xlabel='t, мс',ylabel='Шаг, с',title='Фактический адаптивный шаг')
    elif kind=='Высыхание':
        from .drying import parcel_history,drying_values
        history,s=parcel_history(run.trajectory)
        retention,concentration=drying_values(history,c)
        for j,(x,label) in enumerate([(history['z_m']*1000,'Расстояние по оси z, мм'),(s*1000,'Пройденный путь частицы, мм')],1):
            ax=fig.add_subplot(2,2,j);ax.plot(x,retention*100,color='#176e78');ax.axhline(5,color='#b97739',ls='--',label='5% исходного растворителя')
            ax.set(xlabel=label,ylabel='Осталось растворителя, %',title='Первый материальный сегмент');ax.legend(fontsize=8);ax.grid(alpha=.2)
        ax=fig.add_subplot(223);ax.plot(s*1000,concentration*100);ax.set(xlabel='Пройденный путь, мм',ylabel='Полимер, масс.%');ax.grid(alpha=.2)
        ax=fig.add_subplot(224);ax.plot(s*1000,history['temperature_K']-273.15,color='#b97739');ax.set(xlabel='Пройденный путь, мм',ylabel='Температура, °C');ax.grid(alpha=.2)
        fig.suptitle('Материальная история: после осаждения испарение в основном ядре прекращается'+('' if 'path_length_m' in history.dtype.names else '\nПуть старого архива приближён хордами между снимками'))
    elif kind=='Филаменты и ветвление':
        ax=fig.add_subplot(221);planes=np.linspace(0,c['gap'],60)
        z=live['z_m'];lo=np.minimum(z[:-1],z[1:]);hi=np.maximum(z[:-1],z[1:])
        connected=np.diff(live['bead_id'])==1
        counts=[np.count_nonzero((lo<=p)&(hi>p)&connected) for p in planes]
        ax.step(planes*1000,counts,where='mid');ax.set(xlabel='z, мм',ylabel='Число пересечений',title='Пересечения одной цепочки плоскостью z')
        b=run.branching;ax=fig.add_subplot(222)
        if len(b):
            ax.plot(b['path_length_m']*1000,b['expected_filaments']);ax.set(xlabel='Пройденный путь, мм',ylabel='Ожидаемое N (не целое)',title='Отдельный сценарий дочерних струй')
        else:empty(ax,'Популяционный сценарий не рассчитан.\nРеальное N(z) не определено.')
        ax=fig.add_subplot(223)
        if len(b):
            ax.plot(b['path_length_m']*1000,2*b['radius_per_filament_m']*1e6,label='Диаметр каждой дочерней')
            ax.set(xlabel='Пройденный путь, мм',ylabel='Диаметр, мкм',title='Равное деление потоков в сценарии');ax.legend(fontsize=8)
        ax=fig.add_subplot(224);ax.axis('off')
        empty(ax,'Пересечения и дочерние струи — разные величины.\nСценарий ветвления: популяционный баланс,\nкоэффициент скорости требует калибровки.\nНет дочерних 3D-траекторий и обратной связи.\nСценарий '+('ВКЛЮЧЁН' if c['branching_enabled'] else 'ОТКЛЮЧЁН: N(z) не предсказывается')+'.')
    else:raise ValueError(kind)
    fig.tight_layout()


def export_report(run: RunData, directory: Path):
    from matplotlib.backends.backend_pdf import PdfPages
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    with PdfPages(directory/'figures.pdf') as pdf:
        for i,name in enumerate(PLOT_NAMES):
            fig=Figure(figsize=(12,7),dpi=120);draw(fig,name,run)
            fig.savefig(directory/f'{i+1:02d}.png',dpi=180);pdf.savefig(fig)
        for name,filename in zip(PRIMARY_PLOTS,['fall_2d','drying_z','filaments_z']):
            fig=Figure(figsize=(12,7),dpi=120);draw(fig,name,run)
            fig.savefig(directory/(filename+'.png'),dpi=180)


def draw_study(fig, rows, xkey, ykey, metric):
    fig.clear();good=[r for r in rows if r.get('status')=='COMPLETED' and r.get(metric) is not None]
    if not good:empty(fig.add_subplot(111),'Нет успешных случаев с выбранной метрикой');return
    x=np.array([r[xkey] for r in good]);y=np.array([r[ykey] for r in good]);v=np.array([r[metric] for r in good])
    ax=fig.add_subplot(121);im=ax.scatter(x,y,c=v,cmap='viridis',s=80,edgecolors='white');fig.colorbar(im,ax=ax,label=metric)
    ax.set(xlabel=xkey+' [СИ]',ylabel=ykey+' [СИ]',title='Рассчитанные точки')
    ax=fig.add_subplot(122,projection='3d');ax.scatter(x,y,v,c=v,cmap='viridis',s=35)
    # Surface only on a complete rectangular grid: failed cases are never filled.
    ux=np.unique(x);uy=np.unique(y)
    if len(ux)>1 and len(uy)>1 and len(good)==len(ux)*len(uy):
        zz=np.array([[v[(x==a)&(y==b)][0] for a in ux] for b in uy]);xx,yy=np.meshgrid(ux,uy)
        ax.plot_surface(xx,yy,zz,alpha=.35,cmap='viridis')
    ax.set(xlabel=xkey+' [СИ]',ylabel=ykey+' [СИ]',zlabel=metric,title='Отклик модели; узлы вычислены')
    fig.tight_layout()
