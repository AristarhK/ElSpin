import csv
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from electrospin.field import _solve, solve_field


class FieldTests(unittest.TestCase):
    def test_voltage_polarity_and_charge_balance(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            config = dict(field_nr=48, field_nz=80)
            solve_field(config, root / "positive")
            solve_field({**config, "voltage": -3500}, root / "negative")
            p = np.genfromtxt(root / "positive/potential_rz.csv", delimiter=",", names=True)
            n = np.genfromtxt(root / "negative/potential_rz.csv", delimiter=",", names=True)
            np.testing.assert_allclose(p["Phi_V"], -n["Phi_V"], atol=1e-10)
            self.assertGreater(np.max(p["E_z_V_m"]), 0)
            s = json.loads((root / "positive/field_summary.json").read_text())
            self.assertLess(s["relative_linear_residual"], 1e-9)
            self.assertLess(s["relative_boundary_charge_imbalance"], 1e-8)
            self.assertGreater(s["electrode_charge_C"]["needle"], 0)
            with (root / "positive/electrode_charge.csv").open() as stream:
                material = {r["material"] for r in csv.DictReader(stream)}
            self.assertTrue({"steel", "liquid"}.issubset(material))

    def test_zero_voltage(self):
        with tempfile.TemporaryDirectory() as temp:
            solve_field(dict(voltage=0, field_nr=48, field_nz=80), Path(temp))
            p = np.genfromtxt(Path(temp)/"potential_rz.csv", delimiter=",", names=True)
            self.assertTrue(np.all(p["Phi_V"] == 0))
            self.assertTrue(np.all(p["E_z_V_m"] == 0))

    def test_coaxial_analytical_limit_and_refinement(self):
        # An axial-invariant conductor fills the inner cylinder. End at zmax
        # is grounded in this solver, so use a long cylinder and compare far
        # from that end, where analytic coaxial potential is approached.
        errors = []
        for nr in (40, 80):
            # Geometric radial grid includes electrode boundary exactly.
            rf = np.geomspace(.01, 1., nr)
            rf = np.concatenate([[0.], rf])
            zf = np.linspace(0, 20, 41)
            r = .5*(rf[:-1]+rf[1:])
            mask = np.tile(r < .01, (40, 1))
            result = _solve(rf, zf, mask, 1.)
            exact = np.log(1/r)/np.log(100)
            free = ~mask[0]
            errors.append(np.max(abs(result["phi"][0, free]-exact[free])))
        self.assertLess(errors[1], errors[0]*.6)
        self.assertLess(errors[1], 5e-4)

    def test_invalid_geometry(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(ValueError):
                solve_field(dict(needle_inner_radius=.001), Path(temp))


if __name__ == "__main__":
    unittest.main()
