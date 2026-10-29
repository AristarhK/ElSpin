"""User-facing regressions: physical z, visible small drying changes, time-filtered mat."""
from types import SimpleNamespace
import unittest
import numpy as np
from matplotlib.figure import Figure
from matplotlib.collections import LineCollection
from electrospin.desktopplots import draw,PRIMARY_PLOTS
from electrospin.materials import pva_water
from electrospin.branching import axial_profile


class PrimaryPlots(unittest.TestCase):
    def setUp(self):
        dtype=[(k,float) for k in ['time_s','bead_id','x_m','y_m','z_m','deposited','path_length_m','solvent_mass_kg','polymer_mass_kg']]
        data=np.array([(.01,1,0,0,.02,0,.03,8,2),(.02,1,.001,0,.05,0,.08,4,2)],dtype=dtype)
        dep=np.array([(.01,1,.001,.002,1),(.03,2,.003,.004,1)],dtype=[(k,float) for k in ['time_s','bead_id','collector_x_m','collector_y_m','polymer_mass_kg']])
        branch=np.array([(.03,1),(.08,3)],dtype=[('path_length_m',float),('expected_filaments',float)])
        c=pva_water();c.update(polymer_fraction=.2,branching_enabled=True)
        self.run=SimpleNamespace(config=c,trajectory=data,deposition=dep,branching=branch,times=np.array([.01,.02]),frame=lambda index:data[index:index+1] if index>=0 else data[-1:])

    def test_drying_is_evaporated_percent_against_axial_z(self):
        fig=Figure();draw(fig,PRIMARY_PLOTS[1],self.run)
        ax=fig.axes[0];np.testing.assert_allclose(ax.lines[0].get_xdata(),[20,50])
        np.testing.assert_allclose(ax.lines[0].get_ydata(),[0,50])
        self.assertIn('z, мм',ax.get_xlabel());self.assertIn('Испарилось',ax.get_ylabel())

    def test_drying_excludes_future_history_at_selected_time(self):
        fig=Figure();draw(fig,PRIMARY_PLOTS[1],self.run,index=0)
        np.testing.assert_allclose(fig.axes[0].lines[0].get_xdata(),[20])

    def test_branch_count_uses_axial_coordinate_for_legacy_csv(self):
        fig=Figure();draw(fig,PRIMARY_PLOTS[2],self.run)
        np.testing.assert_allclose(fig.axes[0].lines[0].get_xdata(),[20,50])
        np.testing.assert_allclose(fig.axes[0].lines[0].get_ydata(),[1,3])
        self.assertIn('z, мм',fig.axes[0].get_xlabel())

    def test_new_branch_csv_uses_explicit_z(self):
        b=np.array([(.03,1,.01),(.08,3,.04)],dtype=[('path_length_m',float),('expected_filaments',float),('z_m',float)])
        z,_=axial_profile(self.run.trajectory,b);np.testing.assert_allclose(z,[.01,.04])

    def test_fall_view_has_plane_and_does_not_show_future_impacts(self):
        fig=Figure();draw(fig,PRIMARY_PLOTS[0],self.run)
        side,top=fig.axes
        self.assertIn('Вид сбоку',side.get_title());self.assertIn('Вид сверху',top.get_title())
        np.testing.assert_allclose(side.lines[0].get_ydata(),[100,100])
        self.assertFalse(any(isinstance(c,LineCollection) for c in top.collections))
        self.assertEqual(len(top.collections[0].get_offsets()),1)

    def test_disabled_branching_is_explained_on_the_plot(self):
        self.run.config['branching_enabled']=False
        self.run.branching['expected_filaments']=1
        fig=Figure();draw(fig,PRIMARY_PLOTS[2],self.run)
        self.assertTrue(any('выключен' in t.get_text() and 'N(z) = 1' in t.get_text() for t in fig.axes[0].texts))

    def test_new_results_show_stability_screen_instead_of_flat_count(self):
        self.run.config['branching_enabled']=False
        self.run.branch_screen=np.array([(.02,.3,0),(.05,.7,0)],dtype=[('z_m',float),('electric_capillary_pressure_ratio',float),('limit_model_in_domain',float)])
        fig=Figure();draw(fig,PRIMARY_PLOTS[2],self.run)
        self.assertIn('не прогноз',fig.axes[0].get_title())
        self.assertTrue(any('N(z) не определено' in text.get_text() for text in fig.texts))
