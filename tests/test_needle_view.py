"""Geometry checks for visualization of computed electrode facets."""
import unittest

import numpy as np

from electrospin.needle_view import facet_mesh


class NeedleViewTests(unittest.TestCase):
    def test_cylinder_and_annulus_areas_converge_to_supplied_fv_area(self):
        facets = [
            dict(r_m=.00035, z_m=-.001, normal_r=1, normal_z=0,
                 area_m2=2 * np.pi * .00035 * .0002, sigma_C_m2=3e-5),
            dict(r_m=.000275, z_m=-.0002, normal_r=0, normal_z=1,
                 area_m2=np.pi * (.00035 ** 2 - .0002 ** 2), sigma_C_m2=7e-5),
        ]
        exact = sum(face["area_m2"] for face in facets)
        errors = []
        for sectors in (24, 48):
            points, triangles, values = facet_mesh(facets, sectors)
            vertices = points[triangles]
            area = np.linalg.norm(np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]), axis=1).sum() / 2
            errors.append(abs(area - exact))
            self.assertAlmostEqual(float(points[:, 0].max()), .00035)
            self.assertAlmostEqual(float(points[:, 2].min()), -.0011)
            self.assertAlmostEqual(float(points[:, 2].max()), -.0002)
            self.assertEqual(len(values), len(triangles))
            self.assertEqual(set(values), {3e-5, 7e-5})
        self.assertLess(errors[1], .26 * errors[0])

    def test_invalid_facet_rejected(self):
        with self.assertRaises(ValueError):
            facet_mesh([dict(r_m=0, z_m=0, area_m2=1, normal_r=1, normal_z=0, sigma_C_m2=0)])


if __name__ == "__main__":
    unittest.main()
