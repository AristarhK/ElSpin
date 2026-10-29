from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages
from electrospin.analysis import read_table
from electrospin.desktopplots import RunData,draw,PRIMARY_PLOTS
from electrospin.drying import summarize_drying
ROOT=Path.cwd();OUT=ROOT/'results/physics_report_v04';OUT.mkdir(exist_ok=True)
base=RunData(ROOT/'results/physics_v04_baseline');summary=summarize_drying(base.trajectory,base.config)
rows=[]
for folder in ['physics_audit_v04','physics_material_audit_v04']:
 for row in json.loads((ROOT/'results'/folder/'audit.json').read_text()):
  rows.append({**row,'dataset':folder})
(ROOT/'results/physics_report_v04/audit_combined.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
with PdfPages(OUT/'primary_graphs.pdf') as pdf:
 for name,file in zip(PRIMARY_PLOTS,['fall_2d','drying_z','filaments_screen']):
  fig=Figure(figsize=(12,7));draw(fig,name,base);fig.savefig(OUT/(file+'.png'),dpi=160);pdf.savefig(fig)
fig=Figure(figsize=(11,7));ax=fig.add_subplot(111);bounds=[];zgrid=np.linspace(.002,.1,200)
for row in rows:
 if row['group']!='uncertainty' or row['status']!='COMPLETED':continue
 a=read_table(ROOT/'results'/row['dataset']/row['name']/'drying_history.csv')
 if np.any(np.diff(a['z_m']) < -1e-9):continue
 bounds.append(np.interp(zgrid,a['z_m'],100*(1-a['remaining_solvent_fraction'])))
array=np.array(bounds);ax.fill_between(zgrid*1000,array.min(axis=0),array.max(axis=0),color='#247c9c',alpha=.17,label='Разведочный диапазон R газа и D жидкости')
a=read_table(ROOT/'results/physics_v04_baseline/drying_history.csv');ax.plot(a['z_m']*1000,100*(1-a['remaining_solvent_fraction']),color='#176eaa',lw=2,label='Новая базовая постановка')
a=read_table(ROOT/'results/physics_audit_v04/constant_legacy_transfer/drying_history.csv');ax.plot(a['z_m']*1000,100*(1-a['remaining_solvent_fraction']),color='#a56542',ls='--',label='Старое постоянное k = 0.01 м/с')
ax.set(xlabel='Расстояние от иглы z, мм',ylabel='Испарилось исходного растворителя, %',title='Перенос и неопределённость при заданных свойствах раствора',xlim=(0,100),ylim=(0,None));ax.grid(alpha=.25);ax.legend(loc='upper left')
fig.text(.5,.02,'Полоса — диапазон выбранных гипотез, не доверительный интервал и не экспериментальная погрешность.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.05,1,1));fig.savefig(OUT/'uncertainty.png',dpi=160)
fig=Figure(figsize=(11,7));ax=fig.add_subplot(111)
for radius in [3e-6,1e-6,.2e-6]:
 a=json.loads((ROOT/f'results/physics_audit_v04/reference_a{radius:g}_tol1e-09.json').read_text())
 ax.plot(np.array(a['z_m'])*1000,100*(1-np.array(a['remaining_solvent_fraction'])),lw=2,label=f'a₀={radius*1e6:g} мкм · Q₀={a["inlet_flow_ml_h"]:.4g} мл/ч')
ax.set(xlabel='Координата z, мм',ylabel='Испарилось исходного растворителя, %',title='Независимая проверочная задача: заданная скорость 5 м/с',xlim=(0,100),ylim=(0,100));ax.grid(alpha=.25);ax.legend()
fig.text(.5,.02,'Разные входные радиусы означают разные расходы. Механическая вытяжка и электрические силы в этой задаче не рассчитываются.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.05,1,1));fig.savefig(OUT/'radius_reference.png',dpi=160)
fig=Figure(figsize=(11,8))
for j,group in enumerate(['time','space']):
 p=[r for r in rows if r['group']==group];key='dt_max' if group=='time' else 'injection_interval'
 ax=fig.add_subplot(2,2,j+1);ax.semilogx([r[key]*1e6 for r in p],[100*(1-r['drying_remaining_solvent_fraction_last']) for r in p],'o-')
 ax.set(xlabel=('Шаг времени' if group=='time' else 'Интервал ввода')+', мкс',ylabel='Испарилось к 12 мс, %',title='Высыхание');ax.grid(alpha=.25)
 ax=fig.add_subplot(2,2,j+3);ax.loglog([r[key]*1e6 for r in p],[r['rms_transverse_gain'] for r in p],'o-',color='#9b5e36')
 ax.set(xlabel=('Шаг времени' if group=='time' else 'Интервал ввода')+', мкс',ylabel='Усиление RMS',title='Изгиб: пространственная сходимость не достигнута' if group=='space' else 'Изгиб');ax.grid(alpha=.25)
fig.tight_layout();fig.savefig(OUT/'refinement.png',dpi=160)
reference=json.loads((ROOT/'results/physics_audit_v04/reference_summary.json').read_text())
table='\n'.join(f"{r['radius0_m']*1e6:g} & {r['inlet_flow_ml_h']:.5f} & {100*(1-r['remaining_solvent_at_100mm']):.3f} \\\\" for r in reference if r['rtol']==1e-9)
loss=[100*(1-r['drying_remaining_solvent_fraction_at_first_impact']) for r in rows if r['group']=='uncertainty' and r.get('drying_remaining_solvent_fraction_at_first_impact') is not None]
tex=r'''\documentclass[12pt,a4paper]{article}
\usepackage{fontspec}\setmainfont{DejaVu Serif}\setmonofont{DejaVu Sans Mono}
\usepackage{polyglossia}\setdefaultlanguage{russian}
\usepackage{amsmath,amssymb,geometry,graphicx,booktabs}
\usepackage[unicode,colorlinks=true,linkcolor=black,urlcolor=blue]{hyperref}
\geometry{left=25mm,right=20mm,top=22mm,bottom=22mm}
\graphicspath{{../results/physics_report_v04/}}
\setlength{\parskip}{.4em}\setlength{\parindent}{1.25cm}\sloppy
\title{Электроспиннинг ПВС--вода:\\исправленная постановка высыхания\\\large Верификация, неопределённость и границы расчёта ветвления}
\author{Рабочая расчётная часть курсовой работы}\date{29 октября 2026 г.}
\begin{document}\maketitle
\begin{abstract}
Вместо постоянного коэффициента массоотдачи введены газовая диффузия,
теплопроводность, стефановский поток пара, активность растворителя и
ограничение жидкой диффузией. Для базового режима получена потеря
@@LOSS@@\% воды к первому попаданию. Это численное следствие выбранных
гипотез, а не экспериментальная аттестация состава. Аналитические пределы,
независимая задача переноса и уточнение численного решения проверены.
Сценарный закон рождения дочерних струй не выдаётся за физический прогноз
$N(z)$. Количественная пространственная сходимость изгиба остаётся нерешённой.
\end{abstract}
\section{Постановка и статус входных данных}
Аппаратная схема: неподвижная игла над плоским коллектором по мотивам
Douglas Miller; подача — пользовательский винтовой шприцевой насос по мотивам
MIT Critical Matter. Роботизированная рука не используется. Для ПВС не
предоставлены марка, молекулярная масса, степень гидролиза, измерения реологии
или проводимости. Поэтому свойства задаются гипотезами и диапазонами.
Диапазон чувствительности не является доверительным интервалом.

Прежняя потеря 0.16\% сама по себе не доказывала ошибку. Однако заданное
$k_m=0.01$ м/с не учитывало радиус и граничную задачу газа. При $a=3$ мкм,
$R_g=1$ мм и $D_g=2.5\cdot10^{-5}$ м$^2$/с новая газовая проводимость
даёт $k_m=1.4345$ м/с, $h=1503.4$ Вт/(м$^2$К), с последующим изменением по
мере вытяжки и испарения. Неполное высыхание остаётся допустимым результатом.

\section{Тепломассообмен}
\subsection{Газовая область и активность}
Для стационарной радиальной диффузии между поверхностью цилиндра $r=a$ и
заданной перемешанной границей $r=R_g$ точное решение имеет логарифмический
профиль. Поэтому
\[
 k_m=\frac{D_g}{a\ln(R_g/a)},\qquad
 h=\frac{k_g}{a\ln(R_g/a)}.
\]
$R_g$ является явной неопределённой границей плёночной аппроксимации;
сопряжённое нестационарное поле газа не решается. Установившийся профиль
вокруг движущейся и извивающейся струи — допущение, а не доказанный факт.

Парциальное давление на поверхности $p_s=a_s p_{sat}(T)$, давление в газе
$p_\infty=RH\,p_{sat}(T_\infty)$. Для массовой доли пара:
\[
 Y=\frac{M_s(p/P)}{M_s(p/P)+M_g(1-p/P)},\qquad
 j=\rho_g k_m\ln\frac{1-Y_\infty}{1-Y_s}.
\]
Отрицательный поток при конденсации не используется: текущая модель
описывает испарение. При насыщении и равной температуре поток равен нулю.
В высокополимерном пределе Флори--Хаггинса
\[
 \ln a_s=\ln\phi_s+(1-\phi_s)+\chi(1-\phi_s)^2.
\]
$\chi=0.5$ в базовой гипотезе не является измеренным параметром данного ПВС.
Давление пара вычисляется по Клаузиусу--Клапейрону с постоянной теплотой.

\subsection{Жидкая диффузия и температура}
Внутренний профиль концентрации не разрешается сеткой. Квазистационарное
параболическое замыкание связывает объёмную среднюю и поверхностную доли:
\[
 j_{gas}(\phi_{surface})=\frac{4\rho_s D_l}{a}
 (\phi_{bulk}-\phi_{surface}),\qquad
 D_l=D_{l0}(\phi_{bulk}/\phi_0)^n.
\]
Поверхностная доля находится из скалярного нелинейного уравнения. Это
ограничение сопротивлением диффузии, а не полное решение образования
оболочки, стеклования и кристаллизации. Показатель $n=3$ — гипотеза.

Для материального сегмента:
\[
 (m_p c_p+m_s c_s)\dot T=hS(T_\infty-T)-L\dot m_{ev}+P_J+P_{visc}.
\]
Радиус определяется из объёма и длины. Проверяются масса, масса полимера
и заряд. В файлах сохраняются действующие $k_m$, $h$, $Re$, поверхностная
доля растворителя и поток испарения.

Дополнительная корреляция Churchill--Bernstein задаёт тепло- и массообмен
при поперечном обтекании. Используется нормальная к струе скорость, а
выполнение $Re\,Pr$ и $Re\,Sc\geq0.2$ явно проверяется. В контрольном
расчёте около 86\% сохранённых активных отсчётов лежат вне этого диапазона;
корреляция не объявлена подтверждённой альтернативой для всего пролёта.

\section{Независимая верификация переноса}
Проверены аналитическая проводимость цилиндрического газового слоя,
нулевой поток при насыщении, предел сухого полимера и ограничение жидкой
диффузией. Для адиабатического испарения при постоянной массе полимера
проверено точное соотношение:
\[
 T-T_0=\frac{L}{c_s}\ln\frac{m_pc_p+m_sc_s}{m_pc_p+m_{s0}c_s}.
\]
Скомпилированное ядро также сравнено с независимым решением ОДУ переноса
на свободно движущемся сегменте без сил и вытяжки.

Отдельно рассчитан цилиндр с заданной скоростью 5 м/с, без электрических
сил и механического растяжения. Радиус уменьшается только из-за потери
растворителя. Полимерный поток постоянен; начальные расходы различаются:
\begin{center}\begin{tabular}{rrr}\toprule
$a_0$, мкм & $Q_0$, мл/ч & Испарилось на 100 мм, \%\\\midrule
@@TABLE@@
\bottomrule\end{tabular}\end{center}
Уточнение допуска интегратора от $10^{-5}$ до $10^{-9}$ изменяет конечную
долю растворителя менее чем на $3\cdot10^{-8}$ для этих трёх случаев.
Это верификация численного решения заданной задачи, не экспериментальная
валидация ПВС. Нельзя переносить тонкий случай на исходный расход 0.5 мл/ч,
не пересчитывая скорость и механическую вытяжку.
\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{radius_reference.png}
\caption{Отдельная проверочная задача. Разные радиусы и расходы.}
\end{figure}

\section{Связанный расчёт и диапазоны}
Базовая гипотеза: ПВС 10 масс.\%, 10 кВ, 0.5 мл/ч, зазор 100 мм, начальный
радиус 3 мкм; 40 мс. К первому попаданию потеряно @@LOSS@@\% растворителя.
Выбранный режим предсказывает влажную нить. В разведочном диапазоне
$R_g=0.2$--5 мм и $D_{l0}=10^{-12}$--$10^{-8}$ м$^2$/с потеря составляет
@@RANGE@@\%. Это диапазон выбранных гипотез; ни границы, ни базовые значения
не названы измеренными свойствами конкретной марки ПВС.

Отдельные прогоны меняют $\eta=0.05$--0.8 Па\,с, $\lambda=0.001$--0.1 с,
$K=2\cdot10^{-6}$--$2\cdot10^{-4}$ См/м. Они ограничены общим окном 12 мс;
непопадание на коллектор в этом окне не считается отказом интегратора и не
подменяется диаметром осевшего волокна. Полные конфигурации и статусы сохранены.
\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{uncertainty.png}
\caption{Влияние переноса при общих начальных условиях. Полоса не является вероятностным интервалом.}
\end{figure}

\section{Изгиб и число дочерних струй}
Литературные модели различают осесимметричную капиллярную, электрическую и
изгибную неустойчивости. Изгиб одной непрерывной струи не означает рождения
новых струй. Предыдущий сценарий $dN/dt=\beta(\Pi-\Pi_c)_+N$ не имеет
идентифицированного коэффициента для данного раствора. Поэтому при
отключённом сценарии больше не публикуется плоская линия $N(z)=1$ как
физический расчёт размножения: статус числа дочерних струй — не определено.

Сохраняется локальный скрининг деформации сечения. В пределе проводящего
ньютоновского цилиндра и стоксовского движения для азимутальной моды $m$:
\[
 B_e=\frac{\sigma_q^2a}{\epsilon_0\gamma},\qquad
 s_m=\frac{\gamma m}{2\eta a(m+1)}(B_e-m-1).
\]
При $m=2$ нейтральное отношение электрического и капиллярного давлений
равно 1.5. Это предельный критерий, а не число дочерних филаментов. В коде
проверяются большие $Oh$, малые числа Деборы и зарядовой релаксации,
а также малость продольного напряжения. Для базового пролёта полный набор
условий этого предела не выполняется. Скорость роста не интерпретируется
как скорость ветвления ПВС. Для количественного $N(z)$ нужен закон рождения,
проверенный высокоскоростной съёмкой, либо расчёт изменения топологии
свободной поверхности с соответствующей валидацией.

\section{Численная сходимость и ограничения}
Шум подачи теперь задан на единой физической временной сетке. При
уточнении интервала ввода используется тот же шумовой путь; это исключает
смену случайной реализации из такого сравнения. Уточнение временного шага
стабилизирует расчёт высыхания, но пространственная сходимость амплитуды
изгиба по-прежнему не достигнута. Количественные траектории и диаметры,
зависящие от изгиба, не объявляются аттестованными.
\begin{figure}[ht]\centering\includegraphics[width=\linewidth]{refinement.png}
\caption{Временная и пространственная чувствительность при общей затравке.}
\end{figure}

Для курсовой обоснованным результатом являются уравнения и пределы модели,
аналитическая верификация, независимый расчёт переноса и исследование
чувствительности. Экспериментальное подтверждение и достоверный закон
разветвления остаются отдельными задачами. Исторические результаты версии
0.3 сохранены; настоящий отчёт заменяет их интерпретацию для высыхания.

\section{Литература и доступ к источникам}
Метаданные и аннотации сверены через Crossref; полный текст этих статей
в данной проверке не использован как источник измеренных свойств ПВС.
\begin{enumerate}
\item Churchill, Bernstein. A Correlating Equation for Forced Convection From
Gases and Liquids to a Circular Cylinder in Crossflow (1977).
DOI: \href{https://doi.org/10.1115/1.3450685}{10.1115/1.3450685}.
\item Modeling of solvent evaporation from polymer jets in electrospinning (2011).
DOI: \href{https://doi.org/10.1063/1.3585148}{10.1063/1.3585148}.
В статье рассматривается PAN/DMF; это не база свойств ПВС--вода.
\item Hohman et al. Electrospinning and electrically forced jets. I. Stability
 theory (2001). DOI: 10.1063/1.1383791.
\item Hohman et al. Electrospinning and electrically forced jets. II. Applications
 (2001). DOI: 10.1063/1.1384013.
\item Reneker et al. Bending instability of electrically charged liquid jets of
 polymer solutions in electrospinning (2000). DOI: 10.1063/1.373532.
\end{enumerate}
\end{document}
'''
tex=tex.replace('@@LOSS@@',f"{100*(1-summary['drying_remaining_solvent_fraction_at_first_impact']):.2f}").replace('@@RANGE@@',f'{min(loss):.2f}--{max(loss):.2f}').replace('@@TABLE@@',table)
(ROOT/'docs/PHYSICS_AUDIT_V04_RU.tex').write_text(tex,encoding='utf-8')
print('Cases',len(rows),'loss envelope',min(loss),max(loss))
