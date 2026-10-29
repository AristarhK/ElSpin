from __future__ import annotations
from dataclasses import asdict
import json
from PySide6.QtWidgets import QWidget,QVBoxLayout,QFormLayout,QLineEdit,QComboBox,QLabel,QPushButton,QPlainTextEdit
from .pump import Pump

LABELS={'syringe_inner_diameter_mm':'Внутренний диаметр шприца, мм (измерить)',
 'screw_lead_mm_rev':'Ход винта за оборот, мм (не шаг резьбы!)',
 'motor_steps_rev':'Полных шагов двигателя на оборот','microsteps':'Импульсов STEP на полный шаг',
 'flow_ml_h':'Средний заданный расход, мл/ч','needle_inner_radius_mm':'Внутренний радиус иглы, мм',
 'needle_length_mm':'Длина канала иглы, мм','tubing_inner_radius_mm':'Внутренний радиус трубки, мм',
 'tubing_length_mm':'Длина трубки, мм (0 — без трубки)','viscosity_Pa_s':'Вязкость, Па·с (ньютоновская оценка)',
 'mechanical_efficiency':'КПД винтовой передачи (оценка)'}


class PumpWidget(QWidget):
    def __init__(self,window):
        super().__init__();self.window=window;layout=QVBoxLayout(self)
        label=QLabel('Основа установки: Miller / OpenESpin, игла над плоским коллектором.\nПодача: винтовой шприцевой насос с кареткой по мотивам MIT.\nРазмеры на рисунке не заменяют измерений вашего насоса.');label.setWordWrap(True);layout.addWidget(label)
        form=QFormLayout();self.inputs={}
        for k,v in asdict(Pump()).items():
            if k=='driver':e=QComboBox();e.addItems(['DRV8825','A4988','Другой (проверить паспорт)'])
            else:e=QLineEdit(str(v))
            self.inputs[k]=e;form.addRow('Драйвер' if k=='driver' else LABELS[k],e)
        layout.addLayout(form)
        sync=QPushButton('Взять раствор, расход и иглу из модели');sync.clicked.connect(self.sync);layout.addWidget(sync)
        calc=QPushButton('Рассчитать привод и применить расход к струе');calc.clicked.connect(self.apply);layout.addWidget(calc)
        self.result=QPlainTextEdit();self.result.setReadOnly(True);layout.addWidget(self.result)
        warning=QLabel('Микрошаг задаёт среднюю подачу, но не гарантирует плавный расход.\nМодель струи использует средний расход; пульсации, упругость шприца и трубки требуют измерений.\nΔp и момент — оценки Пуазейля без трения уплотнения; не паспорт двигателя.');warning.setWordWrap(True);layout.addWidget(warning)
        layout.addStretch()

    def sync(self):
        try:
            c=self.window.current_config()
            for k,v in {'flow_ml_h':c['flow_rate']*3600e6,'viscosity_Pa_s':c['viscosity'],
                        'needle_inner_radius_mm':c['needle_inner_radius']*1000,'needle_length_mm':c['needle_length']*1000}.items():self.inputs[k].setText(f'{v:.8g}')
        except Exception as exc:self.window.show_error(str(exc))

    def apply(self):
        try:
            raw={}
            for k,e in self.inputs.items():
                raw[k]=e.currentText() if k=='driver' else float(e.text().replace(',','.'))
            pump=Pump(**raw);p=pump.calculate();c=self.window.current_config()
            c.update(flow_rate=p['flow_rate_m3_s'],needle_inner_radius=raw['needle_inner_radius_mm']*.001,
                     needle_length=raw['needle_length_mm']*.001)
            c.setdefault('_metadata',{})['hardware']={'architecture':'Miller stationary emitter + flat collector; MIT-inspired lead-screw feed',
                'pump_inputs':raw,'pump_result':p,'feed_assumption':'Mean flow only; hydraulic pulsations not coupled.'}
            from .config import validate
            self.window.set_config(validate(c))
            self.result.setPlainText('Средний расход применён к модели.\n'+json.dumps(p,ensure_ascii=False,indent=2)+
                '\nВ прошивке MIT указаны 64 микрошаговых импульса с DRV8825; этот драйвер поддерживает максимум 32. Проверяйте фактическую настройку.')
        except Exception as exc:self.window.show_error(str(exc))
