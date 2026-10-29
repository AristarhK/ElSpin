import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from electrospin.desktop import MainWindow
from electrospin.experiment_ui import ExperimentDialog
from electrospin.materials import pva_water


class ExperimentWindowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def wait_for_fit(self, dialog):
        deadline = time.monotonic()+30
        while dialog.worker.isRunning() and time.monotonic() < deadline:
            self.app.processEvents()
            QTest.qWait(5)
        self.app.processEvents()
        self.assertFalse(dialog.worker.isRunning(), 'Fit thread did not finish')
        self.assertIsNotNone(dialog.result, dialog.status.text())

    def test_main_window_opens_empty_editor_and_preserves_unsaved_series(self):
        w = MainWindow()
        try:
            w.show();self.app.processEvents()
            self.assertTrue(w.experiment_button.isVisible())
            w.experiment_button.click();self.app.processEvents()
            d = w.experiment_dialog
            self.assertTrue(d.isVisible())
            self.assertIsNone(d.result)
            self.assertFalse(d.predict_button.isEnabled())
            self.assertIsNone(d.table.item(0, 0))
            d.fit_button.click()
            self.assertIn('пуста', d.status.text())
            d.set_data([0, 10], [1, 2]);d.close()
            w.experiment_button.click()
            self.assertIs(w.experiment_dialog, d)
            self.assertEqual(d.table.item(1, 1).text(), '2')
        finally:
            w.close();self.app.processEvents()

    def test_excel_fit_query_export_reload_and_condition_invalidation(self):
        d = ExperimentDialog(pva_water())
        try:
            d.show();self.app.processEvents()
            z = np.linspace(0, 100, 11)
            y = (1+(6**1.5-1)*(1-np.exp(-z/30)))**(2/3)
            table = 'z_mm\tN\n'+'\n'.join(f'{a:.10g}\t{b:.10g}' for a,b in zip(z,y))
            QApplication.clipboard().setText(table)
            d.table.setFocus();QTest.keyClick(d.table, Qt.Key_V, Qt.ControlModifier)
            self.assertEqual(len(d.data()[0]), 11)
            d.label.setText('СИНТЕТИЧЕСКИЕ ДАННЫЕ — тест интерфейса')
            d.protocol.setText('Численные тестовые точки, не эксперимент')
            d.bootstrap.setValue(0);d.fit_button.click();self.wait_for_fit(d)
            self.assertEqual(len(d.figure.axes), 2)
            self.assertIn('RMSE=', d.summary.toPlainText())
            d.query.setValue(50);d.predict_button.click()
            self.assertAlmostEqual(d.last_query['mean'], y[5], places=5)
            d.query.setValue(110);d.predict_button.click()
            self.assertIn('Экстраполяция выключена', d.prediction.text())
            d.extrapolate.setChecked(True);d.predict_button.click()
            self.assertTrue(d.last_query['extrapolated'])
            with tempfile.TemporaryDirectory() as tmp:
                with patch('electrospin.experiment_ui.QFileDialog.getExistingDirectory', return_value=tmp):
                    d.export_button.click()
                out = next(Path(tmp).glob('experiment_*'))
                for name in ['model_and_validation.json', 'experimental_series.json', 'forecast.csv',
                             'validation.csv', 'experimental_fit.png', 'experimental_fit.pdf',
                             'summary.txt', 'queried_point.json']:
                    self.assertTrue((out/name).is_file(), name)
                stored = json.loads((out/'model_and_validation.json').read_text())
                self.assertEqual(stored['conditions']['measurement_protocol'], d.protocol.text())
                path = str(out/'experimental_series.json')
                with patch('electrospin.experiment_ui.QFileDialog.getOpenFileName', return_value=(path, '')):
                    d.open_series()
                self.assertIsNone(d.result)  # Reload observations, recompute numerical results.
                self.assertEqual(d.bootstrap.value(), 0)
                self.assertEqual(len(d.data()[0]), 11)
            d.fit_button.click();self.wait_for_fit(d)
            d.conditions['voltage_kV'].setText('12')
            self.assertIsNone(d.result)
            self.assertFalse(d.export_button.isEnabled())
            self.assertFalse(d.predict_button.isEnabled())
        finally:
            d.close();self.app.processEvents()

    def test_interpolation_disables_inapplicable_saturation_settings(self):
        d = ExperimentDialog(pva_water())
        try:
            d.method.setCurrentIndex(d.method.findData('pchip'))
            self.assertFalse(d.bootstrap.isEnabled())
            self.assertFalse(d.fit_onset.isEnabled())
            self.assertFalse(d.onset.isEnabled())
        finally:
            d.close()
