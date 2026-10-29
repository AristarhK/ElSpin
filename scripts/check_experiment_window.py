"""Reproducible GUI check. All generated measurements are explicitly SYNTHETIC."""
from pathlib import Path
import json
import os
import sys
import time

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from PySide6.QtWidgets import QApplication
from electrospin.experiment_ui import ExperimentDialog
from electrospin.materials import pva_water


def main():
    app = QApplication.instance() or QApplication([])
    out = ROOT/'results/experiment_v05'
    out.mkdir(exist_ok=True)
    editor = ExperimentDialog(pva_water())
    try:
        editor.show();app.processEvents()
        editor.grab().save(str(out/'empty_window.png'))
        z = np.linspace(0, 120, 13)
        y = (1+(7**1.5-1)*(1-np.exp(-z/38)))**(2/3)
        y += np.array([0, .13, -.16, .2, -.1, .11, -.18, .12, -.09, .14, -.1, .06, -.04])
        editor.label.setText('СИНТЕТИЧЕСКИЕ ДАННЫЕ — проверка интерфейса, не эксперимент')
        editor.protocol.setText('Сгенерированная тестовая кривая и заданная погрешность 0.15')
        editor.set_data(z, y, np.full(len(z), .15));editor.bootstrap.setValue(80)
        editor.start_fit()
        deadline = time.monotonic()+45
        while editor.worker.isRunning() and time.monotonic()<deadline:
            app.processEvents();time.sleep(.005)
        app.processEvents()
        if editor.worker.isRunning() or editor.result is None:
            raise RuntimeError('Synthetic GUI check failed: '+editor.status.text())
        editor.query.setValue(75);editor.query_prediction();app.processEvents()
        editor.grab().save(str(out/'synthetic_fit_window.png'))
        editor.figure.savefig(out/'synthetic_fit.png', dpi=180)
        editor.figure.savefig(out/'synthetic_fit.pdf')
        (out/'SYNTHETIC_CHECK_ONLY.json').write_text(
            json.dumps(editor.result, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        print(json.dumps({'test_data': 'SYNTHETIC_ONLY_NOT_EXPERIMENTAL',
                          'bootstrap_successful': editor.result['bootstrap_successful'],
                          'withheld_coordinate_rmse': editor.result['validation']['rmse'],
                          'out': str(out)}, ensure_ascii=False))
    finally:
        editor.close();app.processEvents()


if __name__ == '__main__':
    main()
