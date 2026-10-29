"""Numerical and data-integrity checks; synthetic data are not experiments."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from electrospin.calibration import (
    evaluate, fit_dependency, load_experiment, parse_table, predict,
    save_experiment, validate_data,
)


class CalibrationTests(unittest.TestCase):
    def test_excel_decimal_comma_and_metre_csv(self):
        z, y, sigma = parse_table('z_mm;N;sigma\n0;1,0;0,1\n10,5;2,3;0,2\n')
        np.testing.assert_allclose(z, [0, 10.5])
        np.testing.assert_allclose(y, [1, 2.3])
        np.testing.assert_allclose(sigma, [.1, .2])
        z, y, sigma = parse_table('\ufeffz_m,remaining_percent\n0.01,95\n0.02,90\n')
        np.testing.assert_allclose(z, [10, 20])
        self.assertIsNone(sigma)
        np.testing.assert_allclose(parse_table('0\t1\t\n10\t2\t')[1], [1, 2])

    def test_invalid_rows_cannot_silently_disappear(self):
        for text in ['oops;1\n10;2', 'z_mm;N\n1;2;3;4', '1;2;.2\n2;3;', '']:
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_table(text)
        for z, y, target, sigma in [([-1], [1], 'filaments', None),
                                     ([0], [0], 'filaments', None),
                                     ([0], [101], 'remaining_percent', None),
                                     ([0], [np.nan], 'filaments', None),
                                     ([0], [1], 'filaments', [0])]:
            with self.subTest(y=y), self.assertRaises(ValueError):
                validate_data(z, y, target, sigma)

    def test_recovery_and_group_validation_on_exact_synthetic_curve(self):
        z = np.linspace(0, 120, 13)
        known = {'method': 'saturation', 'target': 'filaments',
                 'plateau': 8., 'length_mm': 32., 'onset_mm': 0.}
        result = fit_dependency(z, evaluate(known, z), bootstrap=0)
        self.assertAlmostEqual(result['model']['plateau'], 8, places=6)
        self.assertAlmostEqual(result['model']['length_mm'], 32, places=5)
        self.assertLess(result['validation']['rmse'], 1e-6)
        self.assertEqual(result['validation']['coordinate_groups_used'], len(z))
        self.assertFalse(result['experimental_accuracy_certified'])
        self.assertIsNone(result['mean_band'])  # No fake zero-width uncertainty.
        with self.assertRaisesRegex(ValueError, 'Экстраполяция выключена'):
            predict(result, [121])
        self.assertTrue(1 < predict(result, [200], allow_extrapolation=True)[0] < 8)

    def test_replicates_at_held_out_coordinate_never_train_pchip(self):
        result = fit_dependency([0, 0, 10, 10, 20, 20], [1, 1, 10, 12, 3, 3],
                                method='pchip', sigma=[1, 1, 1, 2, 1, 1], bootstrap=20)
        self.assertAlmostEqual(predict(result, [10])[0], 10.4)
        cv = result['validation']
        self.assertEqual(cv['coordinate_groups_used'], 1)
        self.assertEqual(cv['observations_predicted'], 2)
        np.testing.assert_allclose([p['predicted'] for p in cv['held_out_predictions']], [2, 2])
        self.assertEqual(len(cv['failed_or_unavailable']), 2)
        self.assertIsNone(result['mean_band'])
        with self.assertRaisesRegex(ValueError, 'Экстраполяция недоступна'):
            predict(result, [30], allow_extrapolation=True)

    def test_fit_onset_and_drying_complement_remain_bounded(self):
        z = np.linspace(0, 150, 16)
        lost = 87 * (1 - np.exp(-np.maximum(z-17, 0)/43))
        result = fit_dependency(z, lost, target='evaporated_percent', fit_onset=True, bootstrap=0)
        self.assertAlmostEqual(result['model']['onset_mm'], 17, places=4)
        self.assertAlmostEqual(result['model']['plateau'], 87, places=4)
        remaining = fit_dependency(z, 100-lost, target='remaining_percent', fit_onset=True, bootstrap=0)
        np.testing.assert_allclose(np.array(remaining['mean']) + result['mean'], 100, atol=1e-6)
        values = predict(result, [0, 1e6], allow_extrapolation=True)
        self.assertTrue(np.all((values >= 0) & (values <= 100)))

    def test_constant_data_cannot_identify_future_birth_or_uncertainty(self):
        result = fit_dependency(np.linspace(0, 100, 10), np.ones(10), bootstrap=30)
        self.assertFalse(result['parameters_identified'])
        self.assertIsNone(result['mean_band'])
        with self.assertRaisesRegex(ValueError, 'неидентифицируемые'):
            predict(result, [101], allow_extrapolation=True)
        json.dumps(result, allow_nan=False)

    def test_noisy_bootstrap_reproducible_and_validation_is_separate(self):
        z = np.linspace(0, 120, 12)
        y = 70 * (1 - np.exp(-z/45)) + np.array([0, .8, -1, 1, -.6, .5, -1.2, .5, .6, -.5, .4, -.3])
        options = dict(target='evaporated_percent', sigma=np.ones(len(z)), bootstrap=24, seed=17)
        result = fit_dependency(z, y, **options)
        repeat = fit_dependency(z, y, **options)
        self.assertEqual(result['mean_band'], repeat['mean_band'])
        self.assertEqual(result['bootstrap_successful'], 24)
        self.assertGreater(result['validation']['rmse'], result['training']['rmse'])
        self.assertTrue(np.all(np.array(result['mean_band']['lower']) <= result['mean_band']['upper']))
        self.assertGreater(result['parameter_intervals']['length_mm'][1], result['parameter_intervals']['length_mm'][0])

    def test_sample_count_cancellation_and_persisted_series(self):
        with self.assertRaisesRegex(ValueError, 'не менее 5'):
            fit_dependency([0, 1, 2], [1, 2, 3])
        with self.assertRaises(InterruptedError):
            fit_dependency(np.arange(5), np.arange(5)+1, cancelled=lambda: True)
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'Серия ПВС.json'
            save_experiment(path, [0, 10], [1, 2], 'filaments',
                            conditions={'voltage_kV': 10}, label='Synthetic test',
                            settings={'method': 'pchip'}, result={'arbitrary_stored_fit': True})
            loaded = load_experiment(path)
            self.assertIsNone(loaded['fit'])
            self.assertEqual(loaded['conditions']['voltage_kV'], 10)
            self.assertEqual(loaded['settings']['method'], 'pchip')
            with self.assertRaises(FileExistsError):
                save_experiment(path, [0, 10], [1, 2], 'filaments')
