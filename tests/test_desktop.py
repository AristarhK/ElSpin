import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
from pathlib import Path
import tempfile
import unittest
from PySide6.QtWidgets import QApplication
from electrospin.desktop import MainWindow,ParameterForm,ParametersDialog,Calculation
from electrospin.materials import pva_water
from electrospin.desktopplots import PLOT_NAMES,PRIMARY_PLOTS
from unittest.mock import patch


class DesktopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.app=QApplication.instance() or QApplication([])

    def test_window_units_material_dialog_and_every_plot(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'case';c=pva_water();c.update(duration=.008,gap=.03,field_nr=32,field_nz=64)
            from electrospin.cli import run_case
            run_case(c,out,report=False)
            w=MainWindow();w.load_run(out)
            self.assertAlmostEqual(w.current_config()['voltage'],10000.)
            self.assertEqual(len(w.plot_widgets),len(PLOT_NAMES))
            for i in range(len(w.plot_widgets)):
                w.tabs.setCurrentIndex(i);w.refresh();self.app.processEvents()
                self.assertGreater(len(w.plot_widgets[w.tabs.tabText(i)].figure.axes),0)
            dlg=ParametersDialog(c,material=True,parent=w)
            dlg.form.editors['viscosity'].setText('0,45');dlg.accept_checked()
            self.assertEqual(dlg.c['viscosity'],.45)
            w.close()

    def test_primary_navigation_is_visible_without_horizontal_tabs(self):
        w=MainWindow();w.resize(1280,850);w.show();self.app.processEvents()
        self.assertTrue(w.tabs.tabBar().isHidden())
        for name in PRIMARY_PLOTS:
            button=w.primary_buttons[name];self.assertTrue(button.isVisible());button.click();self.app.processEvents()
            self.assertEqual(w.tabs.tabText(w.tabs.currentIndex()),name)
            self.assertEqual(w.view_title.text(),name);self.assertTrue(button.isChecked())
            self.assertEqual(sum(b.isChecked() for b in w.primary_buttons.values()),1)
        self.assertTrue(w.branch_settings.isVisible())
        w.branch_settings.click();self.assertEqual(w.parameter_tabs.currentIndex(),3)
        w.show_plot('Струя 3D');self.assertGreater(w.other_plots.currentIndex(),0)
        self.assertEqual(sum(b.isChecked() for b in w.primary_buttons.values()),0)
        w.close()

    def test_cancelled_worker_leaves_manifest_and_no_running_child(self):
        with tempfile.TemporaryDirectory() as d:
            from electrospin.cli import run_case
            out=Path(d)/'cancelled'
            with self.assertRaises(InterruptedError):run_case(pva_water(),out,report=False,cancelled=lambda:True)
            import json
            self.assertEqual(json.loads((out/'manifest.json').read_text())['status'],'CANCELLED')

    def test_material_menu_changes_solvent_and_retains_user_material_after_restart(self):
        from electrospin.materials import save_material,load_material
        with tempfile.TemporaryDirectory() as d,patch.dict(os.environ,{'ELSPIN_DATA_DIR':d}):
            w=MainWindow();self.assertEqual(w.material_menu.count(),5)
            w.material_menu.setCurrentText('ПВП–этанол · 12 масс.%');w.select_material(w.material_menu.currentIndex())
            self.assertAlmostEqual(w.current_config()['molar_mass'],.046068)
            self.assertAlmostEqual(w.current_config()['density_solvent'],789.)
            p=Path(d)/'user.json';c=w.current_config();c['viscosity']=.66
            save_material(c,'Мой состав','Измерение','BN в растворе',p)
            w.register_material(load_material(p,c),p);w.close()
            again=MainWindow();self.assertGreaterEqual(again.material_menu.count(),6)
            again.material_menu.setCurrentText('Мой состав');again.select_material(again.material_menu.currentIndex())
            self.assertAlmostEqual(again.current_config()['viscosity'],.66);again.close()
