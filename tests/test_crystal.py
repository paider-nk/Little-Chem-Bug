import sys
from pathlib import Path
import unittest
import tempfile
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gaussian_core import Settings, generate, InputError
from analysis_core import parse_log

class CrystalInputTests(unittest.TestCase):
    def setUp(self):
        self.cif = Path(__file__).resolve().parents[1] / 'examples' / 'crystal_water.cif'
    def test_cross_cell_molecule_is_unwrapped(self):
        result = generate(str(self.cif), 'cif', Settings(job='sp', method='HF', basis='STO-3G', processors=1))
        self.assertIn('H2O', result.formula)
        rows = [x.split() for x in result.text.splitlines() if x.split() and x.split()[0] in ('O','H')]
        self.assertEqual(len(rows), 3)
        self.assertAlmostEqual(float(rows[1][1])-float(rows[0][1]), 0.957, places=2)
    def test_periodic_cell_has_tv_and_rejects_orca(self):
        result = generate(str(self.cif), 'cif', Settings(job='sp', method='HF', basis='STO-3G', processors=1, cif_scope='periodic'))
        self.assertEqual(result.geometry_method, 'CIF 周期晶胞')
        self.assertEqual(sum(1 for x in result.text.splitlines() if x.startswith('Tv')), 3)
        with self.assertRaises(InputError):
            generate(str(self.cif), 'cif', Settings(job='opt', engine='orca', cif_scope='periodic'))
    def test_vmd_opens_final_xyz_from_orca_output(self):
        from vmd_support import locate_vmd, open_vmd
        source = self.cif.with_name('orca_water_excerpt.out')
        def fake_popen(args, cwd, env):
            status = Path(bytes.fromhex(env['WORKFLOW_VMD_STATUS']).decode('utf-8'))
            status.write_text('0 3 1', encoding='ascii')
            self.assertTrue((Path(cwd) / args[-1]).is_file())
            class Process: pid = 1234
            return Process()
        with tempfile.TemporaryDirectory() as directory, patch('vmd_support.subprocess.Popen', side_effect=fake_popen):
            exe = Path(directory)/'vmd.exe'
            exe.touch()
            output = Path(directory)/'water.xyz'
            result = open_vmd(source, vmd=exe, output=output)
            self.assertEqual(output.read_text().splitlines()[0], '3')
            self.assertEqual(result['atoms'], 3)

    def test_reject_disorder_and_infinite_chain(self):
        text = self.cif.read_text(encoding='ascii')
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'specimen.cif'
            source.write_text(text.replace('0.3000 1','0.3000 0.5'), encoding='ascii')
            with self.assertRaisesRegex(InputError,'占位'):
                generate(str(source))
            source.write_text(text[:text.index('O1 O')].replace('10.000(2)','1.5') + 'C1 C 0 0 0 1\n', encoding='ascii')
            with self.assertRaisesRegex(InputError,'无限连接'):
                generate(str(source))

    def test_orca_output_parser(self):
        out = self.cif.with_name('orca_water_excerpt.out')
        a = parse_log(out)
        self.assertTrue(a.normal); self.assertEqual(len(a.atoms), 3); self.assertIsNotNone(a.energy_hartree)

if __name__ == '__main__': unittest.main()

