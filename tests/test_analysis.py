"""Independent checks of map conservation and bead connectivity."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np

from electrospin.analysis import (coverage_map, deposited_segments, field_grid,
                                  read_table, summarize_run)


NAMES = ("time_s", "bead_id", "x_m", "y_m", "collector_x_m", "collector_y_m",
         "radius_m", "dry_radius_m", "polymer_mass_kg", "solvent_mass_kg", "charge_C", "speed_m_s")


def deposits(ids=(1, 2, 3), xy=((0, 0), (.001, 0), (.001, .001))):
    result = np.zeros(len(ids), dtype=[(name, float) for name in NAMES])
    result["bead_id"] = ids
    result["time_s"] = np.arange(len(ids)) * .01
    result["collector_x_m"], result["collector_y_m"] = np.array(xy).T
    result["radius_m"] = 1e-6
    result["dry_radius_m"] = .5e-6
    result["polymer_mass_kg"] = 2e-10
    return result


class AnalysisTests(unittest.TestCase):
    def test_conservative_mass_and_length_map(self):
        data = deposits()
        for bins in (4, 32, 64):
            coverage = coverage_map(data, bins=bins)
            self.assertAlmostEqual(coverage["mass_kg"].sum(), 4e-10, places=20)
            self.assertAlmostEqual(coverage["length_m"].sum(), .002, places=15)
            integral = coverage["areal_polymer_mass_kg_m2"].sum() * coverage["cell_area_m2"]
            self.assertAlmostEqual(integral, 4e-10, places=20)

    def test_never_connect_missing_ids_or_sort_by_time(self):
        data = deposits(ids=(3, 1, 4), xy=((0, 0), (.01, .01), (.001, 0)))
        segments, masses, lengths = deposited_segments(data)
        self.assertEqual(len(segments), 1)
        np.testing.assert_allclose(segments[0], [[0, 0], [.001, 0]])
        self.assertAlmostEqual(float(lengths[0]), .001)

    def test_empty_and_coincident_landings(self):
        empty = deposits(ids=(), xy=np.empty((0, 2)))
        self.assertEqual(float(coverage_map(empty)["mass_kg"].sum()), 0)
        data = deposits(ids=(1, 2), xy=((0, 0), (0, 0)))
        coverage = coverage_map(data)
        self.assertAlmostEqual(float(coverage["mass_kg"].sum()), 2e-10, places=20)
        self.assertEqual(float(coverage["length_m"].sum()), 0)

    def test_field_grid_shuffled_and_incomplete(self):
        field = np.array([(1, 2, 12), (0, 1, 1), (1, 1, 11), (0, 2, 2)],
                         dtype=[("r_m", float), ("z_m", float), ("Phi_V", float)])
        r, z, values = field_grid(field, "Phi_V")
        np.testing.assert_allclose(values, [[1, 11], [2, 12]])
        with self.assertRaises(ValueError):
            field_grid(field[:-1], "Phi_V")

    def test_header_only_csv_and_nonfinite_rejection(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "test.csv"
            path.write_text("a,b\n")
            self.assertEqual(len(read_table(path)), 0)
            path.write_text("a,b\n1,nan\n")
            with self.assertRaises(ValueError):
                read_table(path)
            path.write_text("a,b\n1,2,\n")
            with self.assertRaisesRegex(ValueError, "expected 2 CSV fields, found 3"):
                read_table(path)

    def test_summary_distinguishes_delivered_from_resolved_mass(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with (root / "deposition.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(NAMES)
                for row in deposits():
                    writer.writerow(row.tolist())
            (root / "config.json").write_text("{}")
            result = summarize_run(root)
            self.assertEqual(result["connected_deposited_chords"], 2)
            self.assertAlmostEqual(result["mapped_polymer_mass_kg"], 4e-10, places=20)
            self.assertAlmostEqual(result["delivered_polymer_mass_kg"], 6e-10, places=20)
            self.assertAlmostEqual(result["orientation_nematic_order_0_to_1"], 0, places=10)
            self.assertTrue((root / "coverage.csv").exists())
            self.assertEqual(json.loads((root / "analysis.json").read_text())["deposited_parcels"], 3)

    def test_vtp_exports_final_snapshot_and_only_consecutive_edges(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            names = ["time_s", "bead_id", "x_m", "y_m", "z_m", "radius_m", "charge_C",
                     "temperature_K", "stress_Pa", "deposited"]
            with (root / "trajectory.csv").open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(names)
                writer.writerow([0, 1, 0, 0, 0, 1e-6, 0, 298, 0, 0])
                for bead in (1, 2, 4):
                    writer.writerow([1, bead, bead * .01, 0, .1, 1e-6, 0, 298, 0, 1])
            summarize_run(root)
            tree = ET.parse(root / "trajectory.vtp")
            piece = tree.find("./PolyData/Piece")
            self.assertEqual(piece.attrib["NumberOfPoints"], "3")
            self.assertEqual(piece.attrib["NumberOfLines"], "1")
            connectivity = piece.find("./Lines/DataArray[@Name='connectivity']").text
            self.assertEqual(connectivity, "0 1")


if __name__ == "__main__":
    unittest.main()
