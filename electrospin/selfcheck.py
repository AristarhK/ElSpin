"""Functional desktop/solver check, usable with embedded or frozen Python."""
from pathlib import Path
import tempfile
import json
import os
import time


def main():
    os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
    from PySide6.QtWidgets import QApplication
    from .desktop import MainWindow
    from .cli import run_case,build
    from .materials import pva_water
    from .desktopplots import export_report
    from .calibration import predict
    import numpy as np
    app=QApplication.instance() or QApplication([])
    with tempfile.TemporaryDirectory(prefix='ElSpin-check-') as d:
        c=pva_water();c.update(duration=.008,gap=.03,field_nr=32,field_nz=64)
        out=Path(d)/'run';s=run_case(c,out,executable=build(),report=False)
        if not s['deposited_parcels'] or s['maximum_mass_relative_error']>1e-8:
            raise RuntimeError('Functional solver check failed')
        w=MainWindow();w.load_run(out)
        for i in range(len(w.plot_widgets)):
            w.tabs.setCurrentIndex(i);w.refresh();app.processEvents()
            if not w.plot_widgets[w.tabs.tabText(i)].figure.axes:raise RuntimeError('Empty native plot')
        if w.material_menu.count()<5:raise RuntimeError('Missing material presets')
        w.open_experiment();editor=w.experiment_dialog
        if editor is None or editor.result is not None:raise RuntimeError('Experimental editor must open without invented measurements')
        # Synthetic recovery test only: this is NOT an experimental calibration.
        z=np.linspace(0,100,11);dry=80*(1-np.exp(-z/30))
        editor.target.setCurrentIndex(editor.target.findData('evaporated_percent'))
        editor.set_data(z,dry);editor.bootstrap.setValue(0);editor.start_fit()
        deadline=time.monotonic()+15
        while editor.worker.isRunning() and time.monotonic()<deadline:
            app.processEvents();time.sleep(.005)
        app.processEvents()
        if editor.worker.isRunning() or editor.result is None:
            editor.close();raise RuntimeError('Experimental fitting check failed: '+editor.status.text())
        value=float(predict(editor.result,[50])[0])
        if abs(value-80*(1-np.exp(-50/30)))>1e-5 or len(editor.figure.axes)!=2:
            raise RuntimeError('Synthetic recovery/plot check failed')
        w.close()
        result={'status':'PASSED','deposited_parcels':s['deposited_parcels'],
                'maximum_mass_relative_error':s['maximum_mass_relative_error'],'native_plots':len(w.plot_widgets),'material_presets':5,
                'experimental_editor':'PASSED','calibration_check_data':'synthetic_only',
                'synthetic_recovery_validation_rmse':editor.result['validation']['rmse']}
        print(json.dumps(result));return result
