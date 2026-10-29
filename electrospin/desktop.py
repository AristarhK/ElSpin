"""ElSpin desktop, Qt + Matplotlib, no web service or browser dependency."""
from __future__ import annotations
import csv
from copy import deepcopy
from datetime import datetime
import itertools
import hashlib
import json
import os
from pathlib import Path
import sys
import traceback

from PySide6.QtCore import Qt,QThread,Signal,QTimer,QStandardPaths
from PySide6.QtWidgets import (QApplication,QMainWindow,QWidget,QHBoxLayout,QVBoxLayout,
 QFormLayout,QLineEdit,QLabel,QComboBox,QPushButton,QTabWidget,QScrollArea,QSplitter,
 QFileDialog,QMessageBox,QDialog,QDialogButtonBox,QProgressBar,QSlider,QCheckBox,
 QPlainTextEdit,QTableWidget,QTableWidgetItem,QHeaderView,QSpinBox)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg,NavigationToolbar2QT
from matplotlib.figure import Figure
from .config import DEFAULTS,ROOT,validate,load_config,write_config
from .materials import pva_water,material_presets,MATERIAL_KEYS,save_material,load_material
from .desktopplots import PLOT_NAMES,PRIMARY_PLOTS,PLOT_DESCRIPTIONS,RunData,draw,draw_study,export_report,empty
from .cli import build,run_case

# Display scale: SI = displayed * scale + offset.
FIELDS={
 'voltage':('Напряжение, кВ',1000,0),'gap':('Зазор, см',.01,0),
 'flow_rate':('Расход, мл/ч',1e-6/3600,0),'jet_radius':('Радиус начального сечения, мкм',1e-6,0),
 'emitter_radius':('Наружный радиус иглы, мм',.001,0),'needle_inner_radius':('Внутренний радиус иглы, мм',.001,0),
 'needle_length':('Длина иглы, мм',.001,0),'field_radius':('Радиус камеры, см',.01,0),
 'polymer_fraction':('Полимер, масс.%',.01,0),'viscosity':('Вязкость раствора, Па·с',1,0),
 'relaxation_time':('Время релаксации, мс',.001,0),'surface_tension':('Поверхностное натяжение, мН/м',.001,0),
 'conductivity':('Проводимость, мкСм/м',1e-6,0),'relative_permittivity':('Относительная проницаемость',1,0),
 'current':('Заданный ток на входе, нА',1e-9,0),
 'temperature':('Температура раствора, °C',1,273.15),'ambient_temperature':('Температура воздуха, °C',1,273.15),
 'humidity':('Насыщенность парами растворителя, %',.01,0),
 'density_polymer':('Плотность полимера, кг/м³',1,0),'density_solvent':('Плотность растворителя, кг/м³',1,0),
 'duration':('Время моделирования, мс',.001,0),'dt_max':('Максимальный шаг, мкс',1e-6,0),
 'dt_min':('Минимальный шаг, с',1,0),'injection_interval':('Интервал ввода сегментов, мкс',1e-6,0),
 'output_interval':('Интервал снимков, мс',.001,0),'perturbation_amplitude':('RMS возмущения, мкм',1e-6,0),
 'noise_correlation_time':('Время корреляции шума, мс',.001,0),
 'coulomb_softening':('Регуляризация заряда, мкм',1e-6,0),
 'molar_mass':('Молярная масса растворителя, г/моль',.001,0),'latent_heat':('Теплота испарения, кДж/кг',1000,0),
 'psat_ref':('Давление насыщенного пара, Па',1,0),'psat_temperature':('Опорная температура пара, °C',1,273.15),
 'mass_transfer':('Массоотдача, м/с',1,0),'heat_transfer':('Теплоотдача, Вт/(м²·К)',1,0),
 'cp_polymer':('Теплоёмкость полимера, Дж/(кг·К)',1,0),'cp_solvent':('Теплоёмкость растворителя, Дж/(кг·К)',1,0),
 'finite_extensibility':('Предельный след конформации L²',1,0),
 'viscosity_exponent':('Показатель концентрационной зависимости η',1,0),
 'activation_energy':('Энергия активации вязкости, Дж/моль',1,0),
 'branching_enabled':('Экспертный сценарий N(z), НЕ проверенный прогноз',1,0),
 'branching_rate_s_inv':('Коэффициент ветвления (калибровать), 1/с',1,0),
 'branching_threshold':('Порог электрического/капиллярного давления',1,0),
 'branching_max_expected':('Предел области сценария: ожидаемое N',1,0),
 'transport_model':('Модель тепло- и массообмена',1,0),
 'activity_model':('Активность растворителя',1,0),
 'gas_boundary_radius':('Внешняя граница газовой диффузии, мм',.001,0),
 'vapor_diffusivity':('Диффузия пара в воздухе, м²/с',1,0),
 'liquid_diffusion':('Учитывать жидкую диффузию (параболический профиль)',1,0),
 'liquid_diffusivity':('Диффузия в растворе при c₀, м²/с',1,0),
 'diffusion_concentration_exponent':('Показатель снижения диффузии при высыхании',1,0),
 'flory_huggins_chi':('χ Флори–Хаггинса (гипотеза)',1,0),
 'ambient_pressure':('Давление воздуха, Па',1,0),
}
CHOICES={'field_mode':['needle','uniform','file'],'perturbation_mode':['colored_noise','periodic'],
         'constitutive_model':['fene_p','maxwell'],'drag_model':['cylinder','effective'],
         'transport_model':['gas_film','crossflow','constant'],'activity_model':['flory_huggins','volume_fraction']}


class ParameterForm(QWidget):
    def __init__(self,keys,c):
        super().__init__();self.keys=list(keys);self.editors={};form=QFormLayout(self)
        for key in self.keys:
            label,scale,offset=FIELDS.get(key,(key+' [СИ]',1,0));v=c[key]
            if key in CHOICES:
                e=QComboBox();e.addItems(CHOICES[key]);e.setCurrentText(v)
            elif isinstance(v,bool):e=QCheckBox();e.setChecked(v)
            else:
                text=json.dumps(v) if isinstance(v,list) else v if isinstance(v,str) else f'{(v-offset)/scale:.10g}'
                e=QLineEdit(text);e.setMinimumWidth(110)
            e.setToolTip(key+' · значения сохраняются в СИ');form.addRow(label,e);self.editors[key]=e

    def values(self):
        out={}
        for key,e in self.editors.items():
            default=DEFAULTS[key];_,scale,offset=FIELDS.get(key,(key,1,0))
            try:
                if key in CHOICES:out[key]=e.currentText()
                elif isinstance(default,bool):out[key]=e.isChecked()
                elif isinstance(default,list):out[key]=json.loads(e.text())
                elif isinstance(default,str):out[key]=e.text()
                else:
                    v=float(e.text().strip().replace(',','.'))*scale+offset
                    if isinstance(default,int):
                        if v!=int(v):raise ValueError('нужно целое число')
                        v=int(v)
                    out[key]=v
            except (ValueError,TypeError) as exc:raise ValueError(f'{key}: проверьте значение') from exc
        return out

    def set_values(self,c):
        for key,e in self.editors.items():
            v=c[key];_,scale,offset=FIELDS.get(key,(key,1,0))
            if key in CHOICES:e.setCurrentText(v)
            elif isinstance(v,bool):e.setChecked(v)
            elif isinstance(v,list):e.setText(json.dumps(v))
            elif isinstance(v,str):e.setText(v)
            else:e.setText(f'{(v-offset)/scale:.10g}')


def scrolled(widget):
    scroll=QScrollArea();scroll.setWidgetResizable(True);scroll.setWidget(widget);return scroll


class ParametersDialog(QDialog):
    def __init__(self,c,material=False,parent=None):
        super().__init__(parent);self.c=deepcopy(c);self.material=material;self.setWindowTitle('Новый материал' if material else 'Все параметры')
        self.resize(720,740);layout=QVBoxLayout(self)
        self.name=QLineEdit('Новый материал' if material else c.get('_metadata',{}).get('material_name','Пользовательский состав'))
        self.source=QLineEdit(c.get('_metadata',{}).get('source',''));self.note=QLineEdit(c.get('_metadata',{}).get('note',''))
        if material:
            f=QFormLayout();f.addRow('Название',self.name);f.addRow('Источник свойств / метод измерения',self.source);f.addRow('Полимер, Mw, гидролиз, добавка BN и её доля',self.note);layout.addLayout(f)
            info=QLabel('Для BN введите измеренные эффективные свойства готового раствора.\nОдно название добавки не меняет физические свойства. Все поля обязательны.');info.setWordWrap(True);layout.addWidget(info)
        self.form=ParameterForm(MATERIAL_KEYS if material else DEFAULTS.keys(),c);layout.addWidget(scrolled(self.form))
        self.buttons=QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel);self.buttons.accepted.connect(self.accept_checked);self.buttons.rejected.connect(self.reject);layout.addWidget(self.buttons)

    def accept_checked(self):
        try:
            self.c=validate({**self.c,**self.form.values()})
            if self.material and not self.name.text().strip():raise ValueError('Укажите название материала')
            self.accept()
        except Exception as exc:QMessageBox.warning(self,'Параметры',str(exc))


class PlotWidget(QWidget):
    def __init__(self):
        super().__init__();layout=QVBoxLayout(self);layout.setContentsMargins(2,2,2,2)
        self.figure=Figure(figsize=(10,7));self.canvas=FigureCanvasQTAgg(self.figure)
        layout.addWidget(NavigationToolbar2QT(self.canvas,self));layout.addWidget(self.canvas)


class Calculation(QThread):
    progress=Signal(int,str);result=Signal(str);study_result=Signal(object,str);error=Signal(str)
    def __init__(self,c,out,axes=None):
        super().__init__();self.c=deepcopy(c);self.out=Path(out);self.axes=axes

    def run(self):
        try:
            exe=build()
            if self.axes is None:
                run_case(self.c,self.out,executable=exe,report=False,
                         progress=lambda x:self.progress.emit(int(100*x),'Интегрирование струи'),cancelled=self.isInterruptionRequested)
                self.result.emit(str(self.out));return
            self.out.mkdir(parents=True,exist_ok=False);rows=[]
            cases=list(itertools.product(*self.axes.values()))
            for i,vals in enumerate(cases):
                if self.isInterruptionRequested():break
                changes=dict(zip(self.axes,vals));path=self.out/f'case_{i:03d}'
                row={'case':i,**changes,'status':'FAILED','directory':str(path)}
                try:
                    summary=run_case({**self.c,**changes},path,executable=exe,report=False,
                        progress=lambda x:self.progress.emit(int(100*(i+x)/len(cases)),f'Случай {i+1}/{len(cases)}'),cancelled=self.isInterruptionRequested)
                    row.update({k:v for k,v in summary.items() if v is None or isinstance(v,(str,int,float,bool))});row['status']='COMPLETED'
                except Exception as exc:row.update(status='CANCELLED' if isinstance(exc,InterruptedError) else 'FAILED',error=str(exc))
                rows.append(row)
                (self.out/'study.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
                keys=list(dict.fromkeys(k for r in rows for k in r))
                with (self.out/'study.csv').open('w',newline='',encoding='utf-8-sig') as f:
                    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
            self.study_result.emit(rows,str(self.out))
        except Exception as exc:self.error.emit(str(exc))


class MainWindow(QMainWindow):
    def __init__(self,initial=None):
        super().__init__();self.setWindowTitle('ElSpin 0.5 · электроспиннинг');self.resize(1480,960)
        self.config=pva_water();self.run_data=None;self.worker=None;self.study_rows=[];self.rendered=set();self.experiment_dialog=None
        central=QWidget();self.setCentralWidget(central);outer=QVBoxLayout(central)
        title=QLabel('ElSpin  /  ПВС–вода  /  пространственная динамика струи');title.setStyleSheet('font-size:19px;font-weight:600;color:#1a6f78;padding:8px');outer.addWidget(title)
        splitter=QSplitter();outer.addWidget(splitter,1);left=QWidget();left.setMinimumWidth(350);left.setMaximumWidth(475);l=QVBoxLayout(left)
        l.addWidget(QLabel('Готовые составы / сохранённые материалы'))
        self.material_options=material_presets()
        self.material_library=Path(os.environ.get('ELSPIN_DATA_DIR',QStandardPaths.writableLocation(QStandardPaths.AppDataLocation)))/'materials'
        self.library_errors=[]
        if self.material_library.is_dir():
            for path in sorted(self.material_library.glob('*.json')):
                try:
                    saved=load_material(path,self.config);name=saved['_metadata']['material_name']
                    self.material_options[name]=saved
                except Exception as exc:self.library_errors.append(path.name+': '+str(exc))
        self.material_menu=QComboBox();self.material_menu.addItems(self.material_options)
        self.material_menu.setCurrentText('ПВС–вода · 10 масс.%');self.material_menu.activated.connect(self.select_material);l.addWidget(self.material_menu)
        self.material_label=QLabel();self.material_label.setWordWrap(True);l.addWidget(self.material_label)
        params=QTabWidget();self.parameter_tabs=params;l.addWidget(params,1);self.forms=[]
        groups=[('Установка',['voltage','gap','flow_rate','current','jet_radius','emitter_radius','needle_inner_radius','needle_length','field_radius','ambient_temperature','humidity','field_mode','collector_velocity']),
                ('Раствор',['polymer_fraction','viscosity','relaxation_time','surface_tension','conductivity','relative_permittivity','temperature','density_polymer','density_solvent','constitutive_model','finite_extensibility','evaporation','charge_transport']),
                ('Расчёт',['duration','dt_max','injection_interval','output_interval','max_beads','field_nr','field_nz','perturbation_mode','perturbation_amplitude','noise_correlation_time','random_seed','coulomb_softening','coulomb','image_force','drag_model'])]
        for label,keys in groups:
            form=ParameterForm(keys,self.config);self.forms.append(form);params.addTab(scrolled(form),label)
        branch_form=ParameterForm(['branching_enabled','branching_rate_s_inv','branching_threshold','branching_max_expected'],self.config)
        self.forms.append(branch_form);params.addTab(scrolled(branch_form),'Ветвление')
        from .pump_ui import PumpWidget
        self.pump_widget=PumpWidget(self);params.addTab(scrolled(self.pump_widget),'Насос')
        thermal=ParameterForm(['transport_model','activity_model','gas_boundary_radius','vapor_diffusivity',
            'liquid_diffusion','liquid_diffusivity','diffusion_concentration_exponent','flory_huggins_chi','ambient_pressure'],self.config)
        self.forms.append(thermal);params.addTab(scrolled(thermal),'Испарение')
        for labels,slots in [(['Новый материал…','Импорт материала…'],[self.new_material,self.import_material]),(['Все параметры…','Сброс ПВС–вода'],[self.all_parameters,self.reset]),(['Загрузить JSON…','Сохранить JSON…'],[self.open_config,self.save_config])]:
            line=QHBoxLayout()
            for label,slot in zip(labels,slots):b=QPushButton(label);b.clicked.connect(slot);line.addWidget(b)
            l.addLayout(line)
        self.run_button=QPushButton('Рассчитать струю');self.run_button.setStyleSheet('background:#176e78;color:white;padding:11px;font-weight:600');self.run_button.clicked.connect(self.start_run);l.addWidget(self.run_button)
        self.cancel_button=QPushButton('Остановить');self.cancel_button.setEnabled(False);self.cancel_button.clicked.connect(self.cancel);l.addWidget(self.cancel_button)
        self.progress=QProgressBar();l.addWidget(self.progress);self.status_label=QLabel('Готово к вводу параметров');self.status_label.setWordWrap(True);l.addWidget(self.status_label)
        splitter.addWidget(left);right=QWidget();r=QVBoxLayout(right);splitter.addWidget(right);splitter.setStretchFactor(1,1)
        bar=QHBoxLayout();open_run=QPushButton('Открыть расчёт…');open_run.clicked.connect(self.open_run);bar.addWidget(open_run)
        example=QPushButton('Открыть пример');example.clicked.connect(self.open_example);bar.addWidget(example)
        self.experiment_button=QPushButton('Эксперимент и прогноз…');self.experiment_button.clicked.connect(self.open_experiment);bar.addWidget(self.experiment_button)
        export=QPushButton('Все графики → PNG / PDF…');export.clicked.connect(self.export);bar.addWidget(export)
        self.true_scale=QCheckBox('Единый масштаб 3D');self.true_scale.setChecked(True);self.true_scale.toggled.connect(self.refresh);bar.addWidget(self.true_scale);bar.addStretch();r.addLayout(bar)
        primary_bar=QHBoxLayout();self.primary_buttons={}
        labels=['1. Падение на коллектор\n2D: сбоку и сверху','2. Высыхание раствора\nпо координате z','3. Филаменты\nустойчивость и сценарии']
        for name,label in zip(PRIMARY_PLOTS,labels):
            button=QPushButton(label);button.setCheckable(True);button.setMinimumHeight(52)
            button.setStyleSheet('QPushButton {padding:6px;border:1px solid #a9b8c4;border-radius:5px;} QPushButton:checked {background:#176e78;color:white;border-color:#176e78;font-weight:600;}')
            button.clicked.connect(lambda _checked=False,n=name:self.show_plot(n));primary_bar.addWidget(button,1);self.primary_buttons[name]=button
        r.addLayout(primary_bar)
        secondary=QHBoxLayout();secondary.addWidget(QLabel('Другие графики:'))
        self.other_plot_names=[n for n in PLOT_NAMES if n not in PRIMARY_PLOTS]+['Параметрические расчёты']
        self.other_plots=QComboBox();self.other_plots.addItem('Выберите дополнительный график…');self.other_plots.addItems(self.other_plot_names)
        self.other_plots.activated.connect(lambda i:self.show_plot(self.other_plot_names[i-1]) if i else None);secondary.addWidget(self.other_plots,1)
        self.branch_settings=QPushButton('Настроить ветвление');self.branch_settings.clicked.connect(lambda:self.parameter_tabs.setCurrentIndex(3));secondary.addWidget(self.branch_settings);r.addLayout(secondary)
        self.view_title=QLabel();self.view_title.setStyleSheet('font-size:17px;font-weight:600;color:#154d5f;padding-top:4px');r.addWidget(self.view_title)
        self.view_description=QLabel();self.view_description.setWordWrap(True);r.addWidget(self.view_description)
        frames=QHBoxLayout();self.play_button=QPushButton('▶');self.play_button.setMaximumWidth(45);self.play_button.clicked.connect(self.toggle_animation);frames.addWidget(self.play_button)
        self.slider=QSlider(Qt.Horizontal);self.slider.setRange(0,0);self.slider.valueChanged.connect(self.frame_changed);frames.addWidget(self.slider,1);self.time_label=QLabel('Нет данных');frames.addWidget(self.time_label);r.addLayout(frames)
        self.tabs=QTabWidget();self.plot_widgets={}
        for name in PLOT_NAMES:
            p=PlotWidget();self.plot_widgets[name]=p;self.tabs.addTab(p,name)
        self.tabs.tabBar().hide();self.tabs.currentChanged.connect(self.refresh);r.addWidget(self.tabs,1)
        self.summary=QPlainTextEdit();self.summary.setReadOnly(True);self.summary.setMaximumHeight(105);r.addWidget(self.summary)
        self.study_tab=QWidget();self.tabs.addTab(self.study_tab,'Параметрические расчёты');self.build_study_ui()
        self.timer=QTimer(self);self.timer.setInterval(100);self.timer.timeout.connect(self.advance_frame)
        self.statusBar().showMessage('Свойства стартового раствора требуют измерений. Модель не заменяет эксперимент.')
        if self.library_errors:self.statusBar().showMessage('Не удалось прочитать материалы: '+'; '.join(self.library_errors))
        self.set_config(self.config)
        if initial:self.load_run(Path(initial))
        self.refresh()

    def current_config(self):
        c=deepcopy(self.config)
        for f in self.forms:c.update(f.values())
        return validate(c)

    def set_config(self,c):
        self.config=deepcopy(c)
        for f in self.forms:f.set_values(c)
        self.material_label.setText(c.get('_metadata',{}).get('material_name','Пользовательский раствор')+'\nСвойства: '+c.get('_metadata',{}).get('property_status','проверить источники'))

    def reset(self):self.set_config(pva_water())

    def select_material(self,index):
        try:
            # Material selection keeps machine settings; solvent-vapor saturation
            # is also changed because water RH is not ethanol saturation.
            c=self.current_config();preset=self.material_options[self.material_menu.itemText(index)]
            c.update({k:preset[k] for k in MATERIAL_KEYS})
            c['_metadata']=deepcopy(preset.get('_metadata',{}));c['humidity']=preset['humidity']
            self.set_config(validate(c))
        except Exception as exc:self.show_error(str(exc))

    def register_material(self,c,source_path=None):
        name=c.get('_metadata',{}).get('material_name','Пользовательский раствор')
        self.material_options[name]=deepcopy(c)
        if self.material_menu.findText(name)<0:self.material_menu.addItem(name)
        self.material_menu.setCurrentText(name);self.set_config(c)
        if source_path:
            try:
                content=Path(source_path).read_bytes();self.material_library.mkdir(parents=True,exist_ok=True)
                cached=self.material_library/(hashlib.sha256(content).hexdigest()+'.json')
                if not cached.exists():cached.write_bytes(content)
            except OSError:
                self.statusBar().showMessage('Материал сохранён в выбранном файле; автокаталог недоступен. После перезапуска используйте «Импорт материала».')

    def all_parameters(self):
        try:
            dlg=ParametersDialog(self.current_config(),parent=self)
            if dlg.exec():self.set_config(dlg.c)
        except Exception as exc:self.show_error(str(exc))

    def new_material(self):
        try:
            dlg=ParametersDialog(self.current_config(),material=True,parent=self)
            if not dlg.exec():return
            path,_=QFileDialog.getSaveFileName(self,'Сохранить материал','material.json','JSON (*.json)')
            if not path:return
            save_material(dlg.c,dlg.name.text(),dlg.source.text(),dlg.note.text(),Path(path))
            self.register_material(load_material(Path(path),dlg.c),Path(path))
        except Exception as exc:self.show_error(str(exc))

    def import_material(self):
        path,_=QFileDialog.getOpenFileName(self,'Материал','','JSON (*.json)')
        if path:
            try:self.register_material(load_material(Path(path),self.current_config()),Path(path))
            except Exception as exc:self.show_error(str(exc))

    def open_config(self):
        path,_=QFileDialog.getOpenFileName(self,'Конфигурация','','JSON (*.json)')
        if path:
            try:self.set_config(load_config(Path(path)))
            except Exception as exc:self.show_error(str(exc))

    def save_config(self):
        try:
            c=self.current_config();path,_=QFileDialog.getSaveFileName(self,'Конфигурация','config.json','JSON (*.json)')
            if path:write_config(c,Path(path))
        except Exception as exc:self.show_error(str(exc))

    def output_path(self,prefix):
        parent=QFileDialog.getExistingDirectory(self,'Папка для результатов (внутри создаётся новый каталог)')
        return Path(parent)/(prefix+'_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')) if parent else None

    def start_run(self):
        try:
            c=self.current_config();out=self.output_path('run')
            if out is not None:self.launch(c,out)
        except Exception as exc:self.show_error(str(exc))

    def launch(self,c,out,axes=None):
        if self.worker is not None and self.worker.isRunning():return
        self.timer.stop();self.play_button.setText('▶');self.progress.setValue(0)
        self.worker=Calculation(c,out,axes);self.worker.progress.connect(self.on_progress);self.worker.result.connect(self.load_run)
        self.worker.study_result.connect(self.on_study);self.worker.error.connect(self.show_error);self.worker.finished.connect(self.finished)
        self.run_button.setEnabled(False);self.sweep_button.setEnabled(False);self.cancel_button.setEnabled(True)
        self.status_label.setText('Подготовка электрического поля…');self.worker.start()

    def on_progress(self,value,text):self.progress.setValue(value);self.status_label.setText(text)
    def cancel(self):
        if self.worker:self.worker.requestInterruption();self.status_label.setText('Остановка; частичные данные сохраняются…')
    def finished(self):self.run_button.setEnabled(True);self.sweep_button.setEnabled(True);self.cancel_button.setEnabled(False)
    def show_error(self,text):self.status_label.setText('Операция не завершена');QMessageBox.warning(self,'ElSpin',text)

    def open_run(self):
        path=QFileDialog.getExistingDirectory(self,'Каталог завершённого расчёта')
        if path:self.load_run(path)

    def open_example(self):
        for path in [ROOT/'examples/baseline',ROOT/'results/physics_v04_baseline',ROOT/'results/release_baseline']:
            if (path/'analysis.json').is_file():self.load_run(path);self.show_plot(PRIMARY_PLOTS[0]);return
        self.show_error('Пример не включён в эту папку приложения. Откройте каталог завершённого расчёта или нажмите «Рассчитать струю».')

    def open_experiment(self):
        try:
            if self.experiment_dialog is None:
                from .experiment_ui import ExperimentDialog
                self.experiment_dialog=ExperimentDialog(self.current_config(),self)
            self.experiment_dialog.show();self.experiment_dialog.raise_();self.experiment_dialog.activateWindow()
        except Exception as exc:self.show_error(str(exc))

    def show_plot(self,name):
        for i in range(self.tabs.count()):
            if self.tabs.tabText(i)==name:self.tabs.setCurrentIndex(i);self.refresh();return
        raise ValueError('Unknown plot: '+name)

    def load_run(self,path):
        try:
            data=RunData(path);self.run_data=data;self.rendered.clear();self.timer.stop();self.play_button.setText('▶')
            self.set_config(data.config)
            self.slider.setRange(0,max(0,len(data.times)-1));self.slider.setValue(max(0,len(data.times)-1))
            s=data.summary;diam=s.get('mean_dry_diameter_same_length_at_hits_m');diam_text=f'{diam*1e9:.1f} нм' if diam is not None else 'нет попаданий'
            self.summary.setPlainText(f"Расчёт: {path}\nПопаданий: {s['deposited_parcels']} · сухой диаметр при неизменной длине: {diam_text}\n"
                f"Максимум ошибки массы: {s.get('maximum_mass_relative_error',0):.2e} · усиление RMS: {s.get('rms_transverse_gain')}\n"
                'Происхождение свойств: '+data.config.get('_metadata',{}).get('source','не указано')+'\n'
                'Нет экспериментальной калибровки и подтверждённой пространственной сходимости. Сценарий ветвления отделён от 3D-струи.\n'
                'Полное высыхание, радиальный распад и качество готового материала не подтверждаются одним завершённым расчётом.')
            self.progress.setValue(100);self.status_label.setText('Расчёт загружен');self.refresh()
        except Exception as exc:self.show_error(str(exc))

    def frame_changed(self):
        self.rendered.clear()
        if self.run_data is not None and len(self.run_data.times):self.time_label.setText(f'{self.run_data.times[self.slider.value()]*1000:.2f} мс')
        self.refresh()

    def refresh(self,*_):
        name=self.tabs.tabText(self.tabs.currentIndex())
        self.view_title.setText(name)
        self.view_description.setText(PLOT_DESCRIPTIONS.get(name,'Дополнительный график. К основным результатам можно вернуться тремя кнопками выше.'))
        for label,button in self.primary_buttons.items():button.setChecked(label==name)
        self.other_plots.blockSignals(True)
        self.other_plots.setCurrentIndex(self.other_plot_names.index(name)+1 if name in self.other_plot_names else 0)
        self.other_plots.blockSignals(False)
        self.branch_settings.setVisible(name==PRIMARY_PLOTS[2]);self.true_scale.setVisible(name=='Струя 3D')
        self.slider.setEnabled(name!=PRIMARY_PLOTS[2]);self.play_button.setEnabled(name!=PRIMARY_PLOTS[2])
        if name==PRIMARY_PLOTS[2] and hasattr(self,'timer'):self.timer.stop();self.play_button.setText('▶')
        if name not in self.plot_widgets:return
        if self.run_data is None:
            p=self.plot_widgets[name];p.figure.clear();ax=p.figure.add_subplot(111);ax.axis('off')
            empty(ax,'Нажмите «Открыть пример», «Открыть расчёт…»\nили запустите расчёт кнопкой слева.');p.canvas.draw_idle();return
        # Redraw on explicit changes; cached tabs avoid re-reading large CSVs.
        key=(name,self.slider.value(),self.true_scale.isChecked())
        if key in self.rendered:return
        p=self.plot_widgets[name]
        try:draw(p.figure,name,self.run_data,self.slider.value(),self.true_scale.isChecked());p.canvas.draw_idle();self.rendered.add(key)
        except Exception as exc:self.show_error(str(exc))

    def toggle_animation(self):
        if self.run_data is None:return
        if self.timer.isActive():self.timer.stop();self.play_button.setText('▶')
        else:self.timer.start();self.play_button.setText('Ⅱ')
    def advance_frame(self):self.slider.setValue((self.slider.value()+1)%(self.slider.maximum()+1))

    def export(self):
        if self.run_data is None:return
        out=self.output_path('figures')
        if out:
            try:export_report(self.run_data,out);self.status_label.setText('PNG и многостраничный PDF сохранены: '+str(out))
            except Exception as exc:self.show_error(str(exc))

    def build_study_ui(self):
        layout=QVBoxLayout(self.study_tab);line=QHBoxLayout();self.axis_boxes=[]
        for key,values in [('voltage','8, 10, 12'),('flow_rate','0.3, 0.5, 0.7')]:
            box=QComboBox();box.addItems(['voltage','flow_rate','viscosity','relaxation_time','current','gap','humidity','polymer_fraction','conductivity','random_seed']);box.setCurrentText(key)
            edit=QLineEdit(values);edit.setToolTip('Список чисел в единицах подписи; десятичный разделитель — точка');label=QLabel(FIELDS.get(key,(key,1,0))[0])
            box.currentTextChanged.connect(lambda text,l=label:l.setText(FIELDS.get(text,(text+' [СИ]',1,0))[0]))
            line.addWidget(box);line.addWidget(label);line.addWidget(edit);self.axis_boxes.append((box,edit))
        self.sweep_button=QPushButton('Рассчитать сетку');self.sweep_button.clicked.connect(self.start_sweep);line.addWidget(self.sweep_button);layout.addLayout(line)
        metrics=QHBoxLayout();metrics.addWidget(QLabel('Метрика:'));self.metric_box=QComboBox();self.metric_box.addItems(['mean_dry_diameter_same_length_at_hits_m','rms_impact_radius_material_m','final_max_radial_excursion_m','rms_transverse_gain','deposited_parcels']);self.metric_box.currentTextChanged.connect(self.redraw_study);metrics.addWidget(self.metric_box);metrics.addStretch();layout.addLayout(metrics)
        self.study_plot=PlotWidget();layout.addWidget(self.study_plot,1)
        self.study_table=QTableWidget(0,4);self.study_table.setHorizontalHeaderLabels(['Случай','Параметры [СИ]','Статус','Диагностика']);self.study_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch);self.study_table.setMaximumHeight(170);self.study_table.cellDoubleClicked.connect(self.open_study_case);layout.addWidget(self.study_table)
        layout.addWidget(QLabel('Неудачные случаи сохраняются с причиной и исключаются из поверхности. Двойной щелчок по строке открывает расчёт. До 100 случаев.'))

    def start_sweep(self):
        try:
            c=self.current_config();axes={}
            for box,edit in self.axis_boxes:
                k=box.currentText();_,scale,offset=FIELDS.get(k,(k,1,0))
                if k in axes:raise ValueError('Выберите два разных параметра')
                vals=[float(x.strip())*scale+offset for x in edit.text().split(',') if x.strip()]
                if not vals or len(set(vals))!=len(vals):raise ValueError('Нужны неповторяющиеся значения параметра')
                axes[k]=vals
            if __import__('math').prod(map(len,axes.values()))>100:raise ValueError('Максимум 100 случаев')
            for vals in itertools.product(*axes.values()):validate({**c,**dict(zip(axes,vals))})
            out=self.output_path('sweep')
            if out:self.study_axes=list(axes);self.launch(c,out,axes)
        except Exception as exc:self.show_error(str(exc))

    def on_study(self,rows,path):
        self.study_rows=rows;self.study_table.setRowCount(len(rows))
        for i,row in enumerate(rows):
            text=', '.join(f'{k}={row[k]:.5g}' for k in self.study_axes)
            for j,value in enumerate([row['case'],text,row['status'],row.get('error','')]):self.study_table.setItem(i,j,QTableWidgetItem(str(value)))
        self.redraw_study();self.status_label.setText(f'Серия сохранена: {path}');self.progress.setValue(100 if len(rows) else 0)

    def redraw_study(self,*_):
        if not self.study_rows:return
        draw_study(self.study_plot.figure,self.study_rows,*self.study_axes,self.metric_box.currentText());self.study_plot.canvas.draw_idle()
    def open_study_case(self,row,_column):
        if self.study_rows[row]['status']=='COMPLETED':self.load_run(self.study_rows[row]['directory']);self.tabs.setCurrentIndex(0)

    def closeEvent(self,event):
        if self.experiment_dialog is not None:self.experiment_dialog.close()
        if self.worker is not None and self.worker.isRunning():
            self.worker.requestInterruption();self.worker.wait()
        event.accept()


def main():
    app=QApplication(sys.argv);app.setApplicationName('ElSpin');app.setOrganizationName('ElSpin')
    window=MainWindow(sys.argv[1] if len(sys.argv)>1 else None);window.show();return app.exec()

if __name__=='__main__':raise SystemExit(main())
