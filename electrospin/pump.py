"""Lead-screw syringe feed, independent of any robot arm or network firmware."""
from __future__ import annotations
from dataclasses import dataclass,asdict
import math


@dataclass
class Pump:
    syringe_inner_diameter_mm:float=20. # editable hypothesis, NOT measured from CAD
    screw_lead_mm_rev:float=2.
    motor_steps_rev:int=200
    microsteps:int=16
    driver:str='DRV8825'
    flow_ml_h:float=.5
    needle_inner_radius_mm:float=.2
    needle_length_mm:float=10.
    tubing_inner_radius_mm:float=.5
    tubing_length_mm:float=0.
    viscosity_Pa_s:float=.2
    mechanical_efficiency:float=.3

    def calculate(self):
        values=asdict(self)
        for key,value in values.items():
            if key=='driver':continue
            if not math.isfinite(value) or value<0:raise ValueError(key+': требуется конечное неотрицательное значение')
        for key in ['syringe_inner_diameter_mm','screw_lead_mm_rev','motor_steps_rev','microsteps','flow_ml_h','needle_inner_radius_mm','viscosity_Pa_s','mechanical_efficiency']:
            if values[key]<=0:raise ValueError(key+': должно быть >0')
        if self.motor_steps_rev!=int(self.motor_steps_rev) or self.microsteps!=int(self.microsteps):raise ValueError('Число шагов и делитель должны быть целыми')
        limits={'DRV8825':32,'A4988':16,'Другой (проверить паспорт)':256}
        if self.driver not in limits:raise ValueError('Неизвестный драйвер')
        if self.microsteps not in [1,2,4,8,16,32,64,128,256] or self.microsteps>limits[self.driver]:
            raise ValueError(f'{self.driver}: недопустимый делитель микрошагов. Для DRV8825 максимум 32, A4988 — 16.')
        if self.mechanical_efficiency>1:raise ValueError('КПД должен быть <=1')
        if self.tubing_length_mm>0 and self.tubing_inner_radius_mm<=0:raise ValueError('Укажите положительный радиус трубки')
        area=math.pi*self.syringe_inner_diameter_mm**2/4
        step_mm=self.screw_lead_mm_rev/(self.motor_steps_rev*self.microsteps)
        quantum_ml=area*step_mm/1000
        frequency=self.flow_ml_h/(3600*quantum_ml)
        Q=self.flow_ml_h*1e-6/3600
        resistance=8*self.viscosity_Pa_s*(self.needle_length_mm*.001)/(math.pi*(self.needle_inner_radius_mm*.001)**4)
        if self.tubing_length_mm>0:resistance+=8*self.viscosity_Pa_s*(self.tubing_length_mm*.001)/(math.pi*(self.tubing_inner_radius_mm*.001)**4)
        pressure=resistance*Q;force=pressure*area*1e-6
        return {'flow_rate_m3_s':Q,'plunger_speed_mm_s':self.flow_ml_h*1000/(3600*area),
                'step_travel_mm':step_mm,'volume_per_step_ml':quantum_ml,
                'step_frequency_Hz':frequency,'step_interval_s':1/frequency,
                'pulse_half_period_us':1e6/(2*frequency),'pressure_drop_Pa':pressure,
                'plunger_hydraulic_force_N':force,
                'screw_torque_N_m':force*self.screw_lead_mm_rev*.001/(2*math.pi*self.mechanical_efficiency)}
