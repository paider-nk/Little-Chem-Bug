import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from analysis_core import (Cube, difference_density, energy_difference, matching_grid,
                           parse_fchk, parse_log, solvation_difference, write_report)
from gaussian_core import InputError


class AnalysisTests(unittest.TestCase):
    def test_real_gaussian_water_log(self):
        path = ROOT.parents[1] / "work" / "integration-probe" / "water.log"
        if not path.exists():
            self.skipTest("integration probe is only present in the development workspace")
        result = parse_log(path)
        self.assertTrue(result.normal)
        self.assertAlmostEqual(result.energy_hartree, -74.9456282263, places=8)
        self.assertEqual(len(result.atoms), 3)
        self.assertEqual(len(result.mulliken), 3)
        self.assertAlmostEqual(result.orbitals["alpha"]["gap_ev"], 27.9115, places=3)
        self.assertAlmostEqual(result.dipole_debye, 1.7192, places=3)

    def test_html_report(self):
        path = ROOT.parents[1] / "work" / "integration-probe" / "water.log"
        if not path.exists():
            self.skipTest("integration probe is only present in the development workspace")
        with tempfile.TemporaryDirectory() as tmp:
            report = write_report(parse_log(path), tmp)
            self.assertTrue(report.is_file())
            self.assertTrue((Path(tmp) / "orbital_levels.svg").is_file())
            self.assertTrue((Path(tmp) / "mulliken_charges.csv").is_file())

    def test_real_gaussian_fchk(self):
        path = ROOT.parents[1] / "work" / "integration-probe" / "water.fchk"
        if not path.exists():
            self.skipTest("formatted checkpoint is only present in the development workspace")
        result = parse_fchk(path)
        self.assertEqual(result["Number of atoms"], 3)
        self.assertAlmostEqual(result["Total Energy"], -74.9456282263, places=8)
        self.assertAlmostEqual(result["orbitals"]["alpha"]["gap_ev"], 27.9115, places=3)

    def test_density_difference_shared_grid_and_geometry(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            atoms_ab = [[8, 0, 0, 0, 0], [1, 0, 0, 0, 1], [1, 0, 0, 0, -1]]
            atoms_a, atoms_b = atoms_ab[:1], atoms_ab[1:]
            def save(name, atoms, values):
                from analysis_core import write_cube
                write_cube(Cube([0, 0, 0], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], atoms, values), tmp / name)
            import numpy as np
            save("ab.cube", atoms_ab, np.array([[[5.0]]]))
            save("a.cube", atoms_a, np.array([[[2.0]]]))
            save("b.cube", atoms_b, np.array([[[1.0]]]))
            result = difference_density(tmp / "ab.cube", tmp / "a.cube", tmp / "b.cube", tmp / "delta.cube")
            self.assertEqual(result["integrated_delta_electrons"], 2.0)

    def test_grid_mismatch_rejected(self):
        import numpy as np
        a = Cube(np.zeros(3), np.eye(3), np.empty((0, 5)), np.zeros((1, 1, 1)))
        b = Cube(np.ones(3), np.eye(3), np.empty((0, 5)), np.zeros((1, 1, 1)))
        with self.assertRaises(InputError):
            matching_grid(a, b)


if __name__ == "__main__":
    unittest.main()
