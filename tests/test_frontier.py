import csv
import tempfile
import unittest
from pathlib import Path

from PIL import Image

from analysis_core import HARTREE_EV, parse_log
from frontier_levels import generate_frontier_levels
from gaussian_core import InputError


class FrontierTests(unittest.TestCase):
    def test_gaussian_uses_only_last_block_and_keeps_spin_channels(self):
        with tempfile.TemporaryDirectory() as temp:
            log = Path(temp) / 'open-shell.log'
            latest = (' Alpha  occ. eigenvalues -- -8.0 -0.40\n'
                      ' Alpha virt. eigenvalues --  0.20\n'
                      ' Beta   occ. eigenvalues -- -0.50\n'
                      ' Beta  virt. eigenvalues --  0.10\n')
            log.write_text(' Alpha  occ. eigenvalues -- -1.0\n'
                           ' Alpha virt. eigenvalues --  1.5\n'
                           ' Other SCF step output\n' + latest + ' Normal termination of Gaussian 16\n')
            result = generate_frontier_levels(log)
            alpha, beta = result['rows']
            self.assertEqual((alpha['spin'], beta['spin']), ('alpha', 'beta'))
            self.assertAlmostEqual(alpha['gap_ev'], .60 * HARTREE_EV)
            self.assertAlmostEqual(beta['gap_ev'], .60 * HARTREE_EV)
            self.assertAlmostEqual(beta['lumo_ev'], .10 * HARTREE_EV)
            with result['summary'].open(encoding='utf-8-sig') as stream:
                self.assertEqual([row['spin'] for row in csv.DictReader(stream)], ['alpha', 'beta'])
            with Image.open(result['png']) as figure:
                figure.verify()
            self.assertTrue((result['folder'] / 'orbital_levels.svg').is_file())
            log.write_text(latest + ' Normal termination of Gaussian 16\n'
                           ' Alpha  occ. eigenvalues -- -0.80\n'
                           ' Error termination\n')
            orbitals = parse_log(log).orbitals
            self.assertNotIn('beta', orbitals)
            self.assertIsNone(orbitals['alpha']['lumo_ev'])

    def test_orca_uses_explicit_occupations_and_ignores_other_tables(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'sample.out'
            path.write_text('FINAL SINGLE POINT ENERGY -7.1\nORBITAL ENERGIES\n-----\n'
                            'SPIN UP ORBITALS\n NO OCC E(Eh) E(eV)\n'
                            ' 0 1.0 -0.4 -10.88\n 1 0.0 0.3 8.16\n 2 1.0 -0.2 -5.44\n\n'
                            'SPIN DOWN ORBITALS\n NO OCC E(Eh) E(eV)\n'
                            ' 0 0.0 0.12 3.265\n\n'
                            'Other table\n 1 1.0 0.4 10.88\n'
                            '****ORCA TERMINATED NORMALLY****\n')
            orbitals = parse_log(path).orbitals
            self.assertEqual((orbitals['alpha']['homo_index'], orbitals['alpha']['lumo_index']), (3, 2))
            self.assertIsNone(orbitals['beta']['homo_ev'])
            self.assertAlmostEqual(orbitals['beta']['lumo_ev'], .12 * HARTREE_EV)
            path.write_text('FINAL SINGLE POINT ENERGY -7.1\n Other table\n 0 2.0 -0.4 -10.88\n')
            with self.assertRaises(InputError):
                generate_frontier_levels(path)
