"""Native experimental data entry, condition records, calibration and guarded prediction."""
from __future__ import annotations
from datetime import datetime
import csv,json
from pathlib import Path
import numpy as np
from PySide6.QtCore import Qt,QThread,Signal
from PySide6.QtGui import QKeySequence
from PySide6.QtWidgets import (QApplication,QDialog,QWidget,QVBoxLayout,QHBoxLayout,QFormLayout,
    QLabel,QLineEdit,QComboBox,QPushButton,QTableWidget,QTableWidgetItem,QHeaderView,
    QFileDialog,QSplitter,QCheckBox,QDoubleSpinBox,QSpinBox,QProgressBar,QPlainTextEdit)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg,NavigationToolbar2QT
from .calibration import TARGETS,parse_table,validate_data,fit_dependency,predict,save_experiment,load_experiment


class DataTable(QTableWidget):
    paste_requested=Signal()
    def keyPressEvent(self,event):
        if event.matches(QKeySequence.Paste):
            text=QApplication.clipboard().text()
            if '\t' in text or '\n' in text or ';' in text:self.paste_requested.emit()
            elif self.currentRow()>=0 and self.currentColumn()>=0:self.setItem(self.currentRow(),self.currentColumn(),QTableWidgetItem(text))
            return
        if event.key()==Qt.Key_Delete:
            for item in self.selectedItems():item.setText('')
            return
        super().keyPressEvent(event)


class FitWorker(QThread):
    progress=Signal(int);result=Signal(object);error=Signal(str)
    def __init__(self,z,y,target,options):
        super().__init__();self.z=z.copy();self.y=y.copy();self.target=target;self.options=options
    def run(self):
        try:self.result.emit(fit_dependency(self.z,self.y,self.target,**self.options,
            cancelled=self.isInterruptionRequested,progress=self.progress.emit))
        except Exception as exc:self.error.emit(str(exc))


class ExperimentDialog(QDialog):
    def __init__(self,config,parent=None):
        super().__init__(parent);self.setWindowTitle('ElSpin · Эксперимент и прогноз');self.resize(1320,900)
        self.result=None;self.worker=None;self._loading=False;self.last_query=None
        layout=QVBoxLayout(self)
        title=QLabel('Экспериментальные данные → зависимость → проверка прогноза');title.setStyleSheet('font-size:19px;font-weight:600;color:#176e78');layout.addWidget(title)
        note=QLabel('Одна серия соответствует одному раствору и одному режиму U, Q, влажности и температуры. Зависимость по z не переносится автоматически на другие условия.');note.setWordWrap(True);layout.addWidget(note)
        split=QSplitter();layout.addWidget(split,1);self.data_panel=QWidget();left=QVBoxLayout(self.data_panel);split.addWidget(self.data_panel)
        right=QWidget();r=QVBoxLayout(right);split.addWidget(right);split.setStretchFactor(1,1);self.data_panel.setMaximumWidth(460)
        self.label=QLineEdit();self.label.setPlaceholderText('Название серии и источник измерений')
        self.material=QLineEdit(config.get('_metadata',{}).get('material_name',''))
        self.protocol=QLineEdit();self.protocol.setPlaceholderText('Способ измерения, время усреднения, число кадров/повторов')
        self.target=QComboBox()
        for key,text in TARGETS.items():self.target.addItem(text,key)
        self.method=QComboBox();self.method.addItem('Модель с насыщением','saturation');self.method.addItem('Интерполяция PCHIP (только внутри данных)','pchip')
        form=QFormLayout();form.addRow('Серия / источник',self.label);form.addRow('Материал',self.material);form.addRow('Как измеряли',self.protocol);form.addRow('Величина',self.target);form.addRow('Зависимость',self.method);left.addLayout(form)
        self.conditions={}
        cond=QFormLayout()
        for key,text,value in [('voltage_kV','Напряжение, кВ',config['voltage']/1000),('flow_ml_h','Расход, мл/ч',config['flow_rate']*3600e6),
             ('polymer_percent','Полимер, масс.%',config['polymer_fraction']*100),('humidity_percent','Пары растворителя, %',config['humidity']*100),
             ('ambient_temperature_C','Температура воздуха, °C',config['ambient_temperature']-273.15),('gap_mm','Зазор, мм',config['gap']*1000)]:
            edit=QLineEdit(f'{value:.8g}');self.conditions[key]=edit;cond.addRow(text,edit)
        left.addLayout(cond)
        hint=QLabel('Повторные z допустимы. N — дочерние струи в сечении, не витки укладки. Высыхание — доля исходного растворителя. σ — стандартная неопределённость точки, необязательна.');hint.setWordWrap(True);left.addWidget(hint)
        self.table=DataTable(12,3);self.table.setHorizontalHeaderLabels(['z, мм','Значение','σ (необязательно)'])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch);self.table.paste_requested.connect(self.paste);left.addWidget(self.table,1)
        for labels,slots in [(['Добавить строку','Удалить строки'],[self.add_row,self.remove_rows]),
            (['Загрузить CSV…','Вставить из Excel'],[self.import_csv,self.paste]),
            (['Сохранить серию…','Открыть серию…'],[self.save_series,self.open_series])]:
            line=QHBoxLayout()
            for text,slot in zip(labels,slots):b=QPushButton(text);b.clicked.connect(slot);line.addWidget(b)
            left.addLayout(line)
        options=QFormLayout();self.onset=QDoubleSpinBox();self.onset.setRange(0,1e6);self.onset.setDecimals(3);self.onset.setSuffix(' мм')
        self.fit_onset=QCheckBox('Подбирать начало роста (нужно больше точек)')
        self.bootstrap=QSpinBox();self.bootstrap.setRange(0,1000);self.bootstrap.setValue(200)
        options.addRow('Начало зависимости z₀',self.onset);options.addRow('',self.fit_onset);options.addRow('Bootstrap-повторы',self.bootstrap);left.addLayout(options)
        self.fit_button=QPushButton('Построить зависимость и проверить');self.fit_button.setStyleSheet('background:#176e78;color:white;padding:10px;font-weight:600');self.fit_button.clicked.connect(self.start_fit);left.addWidget(self.fit_button)
        self.figure=Figure(figsize=(10,7));self.canvas=FigureCanvasQTAgg(self.figure);r.addWidget(NavigationToolbar2QT(self.canvas,self));r.addWidget(self.canvas,1)
        progress=QHBoxLayout();self.progress=QProgressBar();progress.addWidget(self.progress,1);self.cancel=QPushButton('Остановить');self.cancel.setEnabled(False);self.cancel.clicked.connect(self.stop);progress.addWidget(self.cancel);r.addLayout(progress)
        self.status=QLabel('Введите измерения, вставьте таблицу из Excel или откройте CSV.');self.status.setWordWrap(True);r.addWidget(self.status)
        self.summary=QPlainTextEdit();self.summary.setReadOnly(True);self.summary.setMaximumHeight(170);r.addWidget(self.summary)
        query=QHBoxLayout();query.addWidget(QLabel('Прогноз при z:'));self.query=QDoubleSpinBox();self.query.setRange(0,1e6);self.query.setDecimals(3);self.query.setSuffix(' мм');query.addWidget(self.query)
        self.extrapolate=QCheckBox('Разрешить экстраполяцию');self.extrapolate.setToolTip('По умолчанию прогноз только внутри диапазона измерений. Неидентифицируемые параметры блокируют экстраполяцию.');query.addWidget(self.extrapolate)
        self.predict_button=QPushButton('Рассчитать точку');self.predict_button.clicked.connect(self.query_prediction);self.predict_button.setEnabled(False);query.addWidget(self.predict_button);r.addLayout(query)
        self.prediction=QLabel('Прогноз появится после подбора.');self.prediction.setWordWrap(True);r.addWidget(self.prediction)
        self.export_button=QPushButton('Экспорт: модель JSON, прогноз CSV, график PNG/PDF…');self.export_button.clicked.connect(self.export);self.export_button.setEnabled(False);r.addWidget(self.export_button)
        for signal in [self.table.itemChanged,self.target.currentIndexChanged,self.method.currentIndexChanged,self.onset.valueChanged,self.fit_onset.toggled,self.bootstrap.valueChanged,
                       self.label.textChanged,self.material.textChanged,self.protocol.textChanged,*[edit.textChanged for edit in self.conditions.values()]]:signal.connect(self.invalidate)
        self.method.currentIndexChanged.connect(self.update_options);self.fit_onset.toggled.connect(self.update_options);self.update_options()
        self._show_empty()

    def _show_empty(self):
        self.figure.clear();ax=self.figure.add_subplot(111);ax.axis('off')
        ax.text(.5,.5,'Сначала введите экспериментальные точки.\nСинтетические данные не подставляются автоматически.',ha='center',va='center',transform=ax.transAxes);self.canvas.draw_idle()

    def invalidate(self,*_):
        if self._loading:return
        self.result=None;self.last_query=None;self.predict_button.setEnabled(False);self.export_button.setEnabled(False)
        self.status.setText('Данные или настройки изменены. Нажмите «Построить зависимость и проверить».')
        self.prediction.setText('Текущий прогноз не рассчитан.');self.summary.clear();self._show_empty()

    def add_row(self):self.table.insertRow(self.table.rowCount())

    def update_options(self,*_):
        saturation=self.method.currentData()=='saturation'
        self.onset.setEnabled(saturation and not self.fit_onset.isChecked());self.fit_onset.setEnabled(saturation);self.bootstrap.setEnabled(saturation)

    def settings_record(self):
        return {'method':self.method.currentData(),'onset_mm':self.onset.value(),'fit_onset':self.fit_onset.isChecked(),'bootstrap':self.bootstrap.value()}
    def remove_rows(self):
        for row in sorted({i.row() for i in self.table.selectedIndexes()},reverse=True):self.table.removeRow(row)
        self.invalidate()

    def set_data(self,z,y,sigma=None):
        self._loading=True;self.table.setRowCount(max(len(z),12));self.table.clearContents()
        for i,(x,value) in enumerate(zip(z,y)):
            for j,v in enumerate([x,value,None if sigma is None else sigma[i]]):
                if v is not None:self.table.setItem(i,j,QTableWidgetItem(f'{v:.10g}'))
        self._loading=False;self.invalidate()

    def data(self):
        lines=[]
        for i in range(self.table.rowCount()):
            row=[self.table.item(i,j).text().strip() if self.table.item(i,j) else '' for j in range(3)]
            if any(row):lines.append('\t'.join(row))
        z,y,sigma=parse_table('\n'.join(lines));return validate_data(z,y,self.target.currentData(),sigma)

    def conditions_record(self):
        result={'material':self.material.text(),'measurement_protocol':self.protocol.text()}
        for key,edit in self.conditions.items():
            try:value=float(edit.text().replace(',','.'))
            except ValueError as exc:raise ValueError('Проверьте условие: '+key) from exc
            if not np.isfinite(value):raise ValueError('Неконечное условие: '+key)
            result[key]=value
        if result['flow_ml_h']<=0 or result['gap_mm']<=0 or not 0<result['polymer_percent']<100 or not 0<=result['humidity_percent']<=100:
            raise ValueError('Проверьте расход, зазор, долю полимера и влажность')
        if result['ambient_temperature_C']<=-273.15:raise ValueError('Температура должна быть выше абсолютного нуля')
        return result

    def paste(self):
        try:self.set_data(*parse_table(QApplication.clipboard().text()))
        except Exception as exc:self.status.setText(str(exc))

    def import_csv(self):
        path,_=QFileDialog.getOpenFileName(self,'Экспериментальные точки','','CSV / TSV (*.csv *.tsv *.txt)')
        if not path:return
        try:self.set_data(*parse_table(Path(path).read_text(encoding='utf-8-sig')))
        except Exception as exc:self.status.setText(str(exc))

    def save_series(self):
        try:
            z,y,sigma=self.data();conditions=self.conditions_record()
            path,_=QFileDialog.getSaveFileName(self,'Сохранить экспериментальную серию','experiment.json','JSON (*.json)')
            if path:
                save_experiment(path,z,y,self.target.currentData(),sigma=sigma,conditions=conditions,label=self.label.text(),result=self.result,settings=self.settings_record())
                self.status.setText('Серия сохранена: '+path)
        except Exception as exc:self.status.setText(str(exc))

    def open_series(self):
        path,_=QFileDialog.getOpenFileName(self,'Экспериментальная серия','','JSON (*.json)')
        if not path:return
        try:
            raw=load_experiment(path);self.target.setCurrentIndex(self.target.findData(raw['target']));self.label.setText(raw.get('label',''))
            c=raw.get('conditions',{});self.material.setText(c.get('material',''))
            self.protocol.setText(c.get('measurement_protocol',''))
            for key,edit in self.conditions.items():
                if key in c:edit.setText(str(c[key]))
            settings=raw.get('settings',{})
            method=self.method.findData(settings.get('method','saturation'))
            if method<0:raise ValueError('Неизвестный метод в файле серии')
            self.method.setCurrentIndex(method);self.onset.setValue(float(settings.get('onset_mm',0)))
            self.fit_onset.setChecked(bool(settings.get('fit_onset',False)));self.bootstrap.setValue(int(settings.get('bootstrap',200)))
            self.set_data(raw['z_mm'],raw['observed'],raw.get('sigma'))
        except Exception as exc:self.status.setText(str(exc))

    def start_fit(self):
        if self.worker is not None and self.worker.isRunning():return
        try:
            z,y,sigma=self.data();self.fit_conditions=self.conditions_record();self.fit_label=self.label.text()
            self.result=None;self.last_query=None;self.prediction.setText('Выполняется подбор…');self.summary.clear();self._show_empty()
            self.predict_button.setEnabled(False);self.export_button.setEnabled(False)
            options={'method':self.method.currentData(),'sigma':sigma,'onset_mm':self.onset.value(),
                     'fit_onset':self.fit_onset.isChecked(),'bootstrap':self.bootstrap.value()}
            self.worker=FitWorker(z,y,self.target.currentData(),options);self.worker.progress.connect(self.progress.setValue)
            self.worker.result.connect(self.show_result);self.worker.error.connect(self.status.setText);self.worker.finished.connect(self.fit_finished)
            self.progress.setValue(0);self.data_panel.setEnabled(False);self.cancel.setEnabled(True)
            self.status.setText('Подбор, проверка исключённых координат и оценка неопределённости…');self.worker.start()
        except Exception as exc:self.status.setText(str(exc))

    def stop(self):
        if self.worker:self.worker.requestInterruption()
    def fit_finished(self):self.data_panel.setEnabled(True);self.cancel.setEnabled(False)

    def show_result(self,result):
        self.result={**result,'conditions':self.fit_conditions,'series_label':self.fit_label,'settings':self.settings_record()}
        self.progress.setValue(100);self.predict_button.setEnabled(True);self.export_button.setEnabled(True)
        self.query.setValue(float(np.mean(result['domain_mm'])));self.extrapolate.setChecked(False)
        self.status.setText('Зависимость построена. Проверяйте ошибку исключённых точек и предупреждения ниже.')
        training=result['training'];cv=result['validation']
        fmt=lambda x:'не определена' if x is None else f'{x:.5g}'
        text=f"Обучение: RMSE={fmt(training['rmse'])}; MAE={fmt(training['mae'])}\n"
        text+=f"Исключённые координаты: {cv['coordinate_groups_used']}/{cv['coordinate_groups_planned']}; RMSE={fmt(cv['rmse'])}; MAE={fmt(cv['mae'])}\n"
        text+=f"Измеренный диапазон: {result['domain_mm'][0]:g}–{result['domain_mm'][1]:g} мм\n"
        if result['method']=='saturation':
            m=result['model'];text+=f"Параметры: предел={m['plateau']:.5g}; длина={m['length_mm']:.5g} мм; начало={m['onset_mm']:.5g} мм\n"
        text+=f"Bootstrap: {result['bootstrap_successful']}/{result['bootstrap_requested']} успешных\n"+'\n'.join(result['warnings'])
        self.summary.setPlainText(text);self.draw_result()

    def draw_result(self,query=None):
        result=self.result
        if result is None:return
        self.figure.clear();ax=self.figure.add_subplot(211);res=self.figure.add_subplot(212)
        x=np.asarray(result['grid_z_mm']);y=np.asarray(result['mean'])
        ax.errorbar(result['z_mm'],result['observed'],yerr=result['sigma'],fmt='o',ms=4,color='#285b69',capsize=3,label='Введённые измерения')
        ax.plot(x,y,color='#aa5d36',lw=2,label='Средняя зависимость')
        band=result['mean_band']
        if band:ax.fill_between(x,band['lower'],band['upper'],color='#cc9876',alpha=.25,label='Условная bootstrap-полоса средней кривой')
        if query is not None:
            q,value=query
            if q>x.max() or q<x.min():
                ext=np.linspace(min(q,x.min()),max(q,x.max()),301);values=predict(result,ext,allow_extrapolation=True)
                outside=(ext<x.min())|(ext>x.max());ax.plot(ext[outside],values[outside],'--',color='#a64147',label='Экстраполяция (гипотеза)')
            ax.scatter([q],[value],s=60,marker='D',color='#a64147',zorder=5,label='Запрошенная точка')
        ax.set(xlabel='Расстояние от иглы z, мм',ylabel=TARGETS[result['target']],title='Эксперимент и подобранная зависимость');ax.legend(fontsize=8);ax.grid(alpha=.25)
        res.axhline(0,color='gray',lw=1);res.scatter(result['z_mm'],result['residuals'],s=16,color='#719397',label='Остатки на обучении')
        held=result['validation']['held_out_predictions']
        if held:res.scatter([h['z_mm'] for h in held],[h['error'] for h in held],marker='x',color='#b45145',label='Ошибка исключённых координат')
        res.set(xlabel='z, мм',ylabel='Измерение − прогноз',title='Проверка ошибки');res.grid(alpha=.25);res.legend(fontsize=8)
        label=result.get('series_label','')
        if label:self.figure.suptitle(label[:120],fontsize=10)
        self.figure.tight_layout();self.canvas.draw_idle()

    def query_prediction(self):
        if self.result is None:return
        try:
            q=self.query.value();value=float(predict(self.result,[q],allow_extrapolation=self.extrapolate.isChecked())[0])
            outside=q<self.result['domain_mm'][0] or q>self.result['domain_mm'][1]
            text=f'z={q:g} мм → {TARGETS[self.result["target"]]}: {value:.5g}'
            if outside:text+=' · ЭКСТРАПОЛЯЦИЯ: точность вне измерений не проверена'
            elif self.result['mean_band']:
                band=self.result['mean_band'];x=self.result['grid_z_mm'];low=np.interp(q,x,band['lower']);high=np.interp(q,x,band['upper'])
                text+=f' · условная полоса средней кривой [{low:.5g}; {high:.5g}]'
            self.last_query={'z_mm':q,'mean':value,'extrapolated':outside}
            self.prediction.setText(text);self.draw_result((q,value))
        except Exception as exc:self.last_query=None;self.prediction.setText(str(exc));self.draw_result()

    def export(self):
        if self.result is None:return
        parent=QFileDialog.getExistingDirectory(self,'Папка для экспорта (создаётся новый каталог)')
        if not parent:return
        try:
            out=Path(parent)/('experiment_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));out.mkdir()
            (out/'model_and_validation.json').write_text(json.dumps(self.result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
            save_experiment(out/'experimental_series.json',self.result['z_mm'],self.result['observed'],self.result['target'],
                sigma=self.result['sigma'],conditions=self.result['conditions'],label=self.result['series_label'],result=self.result,settings=self.result['settings'])
            if self.last_query:(out/'queried_point.json').write_text(json.dumps(self.last_query,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
            with (out/'forecast.csv').open('w',newline='',encoding='utf-8-sig') as f:
                w=csv.writer(f);w.writerow(['z_mm','mean','mean_band_lower','mean_band_upper'])
                band=self.result['mean_band']
                for i,(z,y) in enumerate(zip(self.result['grid_z_mm'],self.result['mean'])):w.writerow([z,y,band['lower'][i] if band else '',band['upper'][i] if band else ''])
            with (out/'validation.csv').open('w',newline='',encoding='utf-8-sig') as f:
                w=csv.DictWriter(f,fieldnames=['z_mm','observed','predicted','error']);w.writeheader();w.writerows(self.result['validation']['held_out_predictions'])
            self.figure.savefig(out/'experimental_fit.png',dpi=180);self.figure.savefig(out/'experimental_fit.pdf')
            (out/'summary.txt').write_text(self.summary.toPlainText(),encoding='utf-8');self.status.setText('Экспорт сохранён: '+str(out))
        except Exception as exc:self.status.setText(str(exc))

    def closeEvent(self,event):
        if self.worker is not None and self.worker.isRunning():self.worker.requestInterruption();self.worker.wait()
        event.accept()
