"""Regression coverage for the optional preprocessing and complex builder."""
import json
import os
from pathlib import Path
import sys
import subprocess
import unittest
from unittest.mock import patch

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gaussian_core import InputError, Settings, generate, _geometry, _trigonal_boron_dihedrals, _xtb_optimize
from complex_core import PRESETS, validate_spec
from gaussian_input import cli

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


class PreprocessingTests(unittest.TestCase):
    def test_preprocess_method_choices_and_external_xtb_dispatch(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        with patch("gaussian_core._xtb_optimize") as optimize:
            cid, method, energy, notes = _geometry(mol, 1, "xtb", 0, 1)
        optimize.assert_called_once_with(mol, cid, "xtb", 0, 1, 1)
        self.assertIn("GFN2", method)
        self.assertIsNone(energy)
        self.assertTrue(any("xTB" in note for note in notes))
        cid, method, energy, notes = _geometry(Chem.AddHs(Chem.MolFromSmiles("CCO")), 1, "none", 0, 1)
        self.assertIn("ETKDGv3", method)
        self.assertTrue(any("跳过力场" in note for note in notes))

    def test_explicit_force_field_missing_parameters_is_not_silently_replaced(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
        with patch("gaussian_core.AllChem.MMFFHasAllMoleculeParams", return_value=False):
            with self.assertRaisesRegex(InputError, "MMFF94s"):
                _geometry(mol, 1, "mmff94s", 0, 1)
        with patch("gaussian_core.AllChem.UFFHasAllMoleculeParams", return_value=False):
            with self.assertRaisesRegex(InputError, "UFF"):
                _geometry(mol, 1, "uff", 0, 1)

    def test_xtb_selection_requires_an_executable(self):
        with patch.dict(os.environ, {"XTB_EXE": ""}), \
                patch("gaussian_core._user_xtb_executable", return_value=""), \
                patch("gaussian_core.shutil.which", return_value=None), \
                self.assertRaisesRegex(InputError, "找不到 xtb"):
            generate("CCO", settings=Settings(conformers=1, preprocess_method="xtb"))

    def test_tblite_cli_and_template_use_selected_backend(self):
        import tempfile
        with patch("gaussian_core._xtb_optimize") as optimize, tempfile.TemporaryDirectory() as directory:
            output, meta = Path(directory) / "ethanol.gjf", Path(directory) / "ethanol.json"
            self.assertEqual(cli(["CCO", "--preprocess-method", "tblite", "--processors", "3",
                                  "--conformers", "1", "-o", str(output), "--metadata", str(meta)]), 0)
            self.assertEqual(optimize.call_args.args[2:], ("tblite", 0, 1, 3))
            data = json.loads(meta.read_text(encoding="utf-8"))
            self.assertEqual(data["settings"]["preprocess_method"], "tblite")
            self.assertIn("tblite", data["geometry_method"])
        with patch("gaussian_core._xtb_optimize") as optimize:
            result = generate("ferrocene", "ferrocene", Settings(preprocess_method="tblite", processors=2))
            self.assertEqual(optimize.call_args.args[2:], ("tblite", 0, 1, 2))
            self.assertIn("tblite", result.geometry_method)
            self.assertIn("Geom=Connectivity", result.text)

    def test_tblite_command_threads_charge_spin_and_convergence(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("[CH3]"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        def run(command, **kwargs):
            self.assertEqual(command[2:], ["--gfn", "2", "--opt", "tight", "--tblite",
                                          "--chrg", "0", "--uhf", "1", "--parallel", "3"])
            for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
                self.assertEqual(kwargs["env"][variable], "3")
            folder = Path(kwargs["cwd"])
            (folder / "xtbopt.xyz").write_text((folder / "input.xyz").read_text())
            (folder / ".xtboptok").touch()
            return subprocess.CompletedProcess(command, 0, "optimization converged", "")
        with patch("gaussian_core._xtb_executable", return_value="xtb.exe"), \
                patch("gaussian_core.subprocess.run", side_effect=run) as runner:
            _xtb_optimize(mol, 0, "tblite", 0, 2, 3)
            runner.assert_called_once()

    def test_tblite_rejects_unconverged_or_invalid_coordinates_without_mutation(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("O"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        before = mol.GetConformer().GetPositions().copy()
        for bad in ("unconverged", "not_converged", "reordered", "nan", "truncated"):
            def run(command, **kwargs):
                folder = Path(kwargs["cwd"])
                text = (folder / "input.xyz").read_text()
                if bad == "reordered":
                    text = text.replace("O ", "C ")
                elif bad in ("nan", "truncated"):
                    rows = text.splitlines()
                    rows[-1] = "H nan 0 0" if bad == "nan" else "H 0 0"
                    text = "\n".join(rows)
                (folder / "xtbopt.xyz").write_text(text)
                if bad != "unconverged":
                    (folder / ".xtboptok").touch()
                if bad == "not_converged":
                    (folder / "NOT_CONVERGED").touch()
                return subprocess.CompletedProcess(command, 0, "", "")
            with self.subTest(bad=bad), patch("gaussian_core._xtb_executable", return_value="xtb.exe"), \
                    patch("gaussian_core.subprocess.run", side_effect=run), self.assertRaises(InputError):
                _xtb_optimize(mol, 0, "tblite", 0, 1)
            np.testing.assert_array_equal(before, mol.GetConformer().GetPositions())

    def test_tblite_unavailable_does_not_fall_back(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("O"))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        failure = subprocess.CompletedProcess([], 1, "", "not compiled with tblite support")
        with patch("gaussian_core._xtb_executable", return_value="xtb.exe"), \
                patch("gaussian_core.subprocess.run", return_value=failure) as runner, \
                self.assertRaisesRegex(InputError, "不支持 tblite"):
            _xtb_optimize(mol, 0, "tblite", 0, 1)
        runner.assert_called_once()

    def test_skip_preserves_coordinates_without_planarization_force_field(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1B(O)O"))
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=42), 0)
        before = mol.GetConformer().GetPositions().copy()
        notes = []
        with patch("gaussian_core.AllChem.UFFGetMoleculeForceField", side_effect=AssertionError("must skip UFF")):
            constraints = _trigonal_boron_dihedrals(mol, mol.GetConformer(), "skip", notes=notes)
        np.testing.assert_array_equal(before, mol.GetConformer().GetPositions())
        self.assertEqual(len(constraints), 2)
        self.assertEqual(constraints[0][-1], 180.0)
        self.assertTrue(any("已跳过" in note for note in notes))

    def test_user_macrocycle_skip_keeps_every_ring_bond_and_freezes_both_boron_centres(self):
        source = json.loads((EXAMPLES / "boron_macrocycle_skip.json").read_text(encoding="utf-8"))
        result = generate(source["input"], "smiles", Settings(conformers=1, boron_preprocess="skip"))
        self.assertEqual((result.formula, result.atom_count), ("C102H84B2N2", 190))
        parts = result.text.strip().split("\n\n")
        self.assertEqual(len(parts), 5)  # header/title/coordinates/connectivity/torsions
        self.assertIn("Geom=Connectivity", parts[0])
        self.assertIn("Opt=(ModRedundant)", parts[0])
        self.assertEqual(len(parts[2].splitlines()) - 1, 190)
        actual = {}
        for row in parts[3].splitlines():
            fields = row.split()
            for other, order in zip(fields[1::2], fields[2::2]):
                pair = (int(fields[0]) - 1, int(other) - 1)
                self.assertNotIn(pair, actual)
                actual[pair] = float(order)
        mol = Chem.AddHs(Chem.MolFromSmiles(source["input"]))
        expected = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))): b.GetBondTypeAsDouble()
                    for b in mol.GetBonds()}
        self.assertEqual(actual, expected)
        torsions = parts[4].splitlines()
        self.assertEqual(len(torsions), 8)
        for value, frozen in zip(torsions[::2], torsions[1::2]):
            self.assertTrue(value.startswith("D "))
            self.assertEqual(value.rsplit(" ", 1)[0] + " F", frozen)
        self.assertNotIn("\nA ", result.text)

    def test_macrocycle_projection_preserves_bonds_and_makes_both_boron_centers_planar(self):
        source = json.loads((EXAMPLES / "boron_macrocycle_tblite.json").read_text(encoding="utf-8"))
        value = source["input"]
        def sections(mode):
            return generate(value, "smiles", Settings(conformers=1, preprocess_method="none",
                            boron_preprocess=mode)).text.strip().split("\n\n")
        before, after = sections("skip"), sections("projected")
        self.assertEqual(before[3:], after[3:])
        self.assertNotIn("\nA ", after[4])
        def positions(parts):
            return np.array([[float(x) for x in row.split()[1:]] for row in parts[2].splitlines()[1:]])
        xyz0, xyz = positions(before), positions(after)
        mol = Chem.AddHs(Chem.MolFromSmiles(value))
        edges = np.array([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])
        self.assertEqual((len(xyz), len(edges)), (190, 207))
        a, b = edges.T
        drift = np.abs(np.linalg.norm(xyz[a] - xyz[b], axis=1) - np.linalg.norm(xyz0[a] - xyz0[b], axis=1))
        self.assertLessEqual(drift.max(), 0.03)
        conf = Chem.Conformer(len(xyz))
        for i, pos in enumerate(xyz):
            conf.SetAtomPosition(i, pos)
        for line in after[4].splitlines()[::2]:
            _, a, b, c, d, target = line.split()
            measured = rdMolTransforms.GetDihedralDeg(conf, *(int(x) - 1 for x in (a, b, c, d)))
            self.assertLessEqual(abs((measured - float(target) + 180) % 360 - 180), 0.5)
        centers = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 5 and a.GetDegree() == 3]
        self.assertEqual(len(centers), 2)
        for center in centers:
            a, b, c = xyz[[n.GetIdx() for n in center.GetNeighbors()]]
            normal = np.cross(b - a, c - a)
            self.assertGreater(np.linalg.norm(normal), 1e-8)
            height = abs(np.dot(xyz[center.GetIdx()] - a, normal) / np.linalg.norm(normal))
            self.assertLessEqual(height, 0.001)

    def test_invalid_preprocessing_options_fail_before_embedding(self):
        for settings in (Settings(boron_preprocess="unknown"), Settings(boron_bond_tolerance=float("nan")),
                         Settings(boron_bond_tolerance=0), Settings(boron_bond_tolerance=1.01)):
            with self.subTest(settings=settings), self.assertRaises(InputError):
                generate("O", settings=settings)

    def test_cli_passes_preprocessing_controls(self):
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            output, meta = Path(directory) / "b.gjf", Path(directory) / "b.json"
            self.assertEqual(cli(["CCB(O)O", "--boron-preprocess", "skip", "--boron-bond-tolerance", "0.2",
                                  "--conformers", "1", "-o", str(output), "--metadata", str(meta)]), 0)
            data = json.loads(meta.read_text(encoding="utf-8"))
            self.assertEqual(data["settings"]["boron_preprocess"], "skip")
            self.assertEqual(data["settings"]["boron_bond_tolerance"], 0.2)
            self.assertIn("180.0", output.read_text())

    def test_relaxed_cli_works_in_legacy_windows_console(self):
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "relaxed.gjf"
            run = subprocess.run([sys.executable, str(EXAMPLES.parent / "gaussian_input.py"),
                                  "c1ccccc1B(O)O", "--boron-preprocess", "relaxed",
                                  "--conformers", "1", "-o", str(output)],
                                 env={**os.environ, "PYTHONIOENCODING": "gbk"}, capture_output=True, timeout=30)
            self.assertEqual(run.returncode, 0, run.stderr.decode("gbk", errors="replace"))
            self.assertIn("Opt=(ModRedundant)", output.read_text())
            self.assertIn("\\xc5", run.stdout.decode("gbk"))


class ComplexTests(unittest.TestCase):
    def test_custom_multidentate_ligand_donors_are_required_and_validated(self):
        spec = {"metal": "Fe", "oxidation_state": 2, "geometry": "oct",
                "ligands": [{"kind": "smiles", "value": "NCCN", "count": 3, "donors": [1, 4]}]}
        self.assertEqual(validate_spec(spec), spec)
        for donors in ([], [0, 4], [1, 5], [1, 1], [True]):
            with self.subTest(donors=donors), self.assertRaises(InputError):
                validate_spec({**spec, "ligands": [{**spec["ligands"][0], "donors": donors}]})

    def test_saved_backend_topology_has_correct_coordination_and_original_ligand_bonds(self):
        for stem, preset, count in (("cu_ammonia4", "四氨合铜(II)", 4),
                                    ("ferrocene_molsimplify", "二茂铁（π 配位）", 10)):
            with self.subTest(stem=stem):
                spec, charge, spin = PRESETS[preset]
                # The fixture contains real molSimplify output, so tests need no optional runtime.
                text = (EXAMPLES / f"{stem}.gjf").read_text()
                parts = text.strip().split("\n\n")
                atoms = [row.split() for row in parts[2].splitlines()[1:]]
                metal = next(i for i, atom in enumerate(atoms) if atom[0] == spec["metal"])
                bonds = []
                for row in parts[3].splitlines():
                    fields = row.split()
                    bonds.extend((int(fields[0])-1, int(j)-1, float(order))
                                 for j, order in zip(fields[1::2], fields[2::2]))
                payload = {"version": "2.0.0", "metal_index": metal,
                           "atoms": [{"symbol": a[0], "xyz": list(map(float, a[1:]))} for a in atoms],
                           "bonds": bonds}
                with patch("complex_core.build_geometry", return_value=payload):
                    result = generate(json.dumps(spec), "complex", Settings(charge=charge, multiplicity=spin, basis="Def2SVP"))
                self.assertEqual(result.text.strip().split("\n\n")[3], parts[3])
                self.assertEqual(len([b for b in bonds if metal in b[:2]]), count)
                self.assertTrue(all(b[2] == 0.1 for b in bonds if metal in b[:2]))
                self.assertEqual((result.charge, result.multiplicity), (charge, spin))
                with patch("complex_core.build_geometry", return_value=payload), \
                        patch("gaussian_core._xtb_optimize") as optimize:
                    refined = generate(json.dumps(spec), "complex", Settings(charge=charge, multiplicity=spin,
                                       basis="Def2SVP", preprocess_method="tblite", processors=2))
                self.assertEqual(optimize.call_args.args[2:], ("tblite", charge, spin, 2))
                self.assertIn("tblite", refined.geometry_method)
                self.assertEqual(refined.text.strip().split("\n\n")[3], parts[3])
                with patch("complex_core.build_geometry", return_value=payload), self.assertRaisesRegex(InputError, "不相容"):
                    generate(json.dumps(spec), "complex", Settings(charge=charge, multiplicity=spin + 1))
                with patch("complex_core.build_geometry", side_effect=AssertionError("missing spin must fail first")):
                    with self.assertRaisesRegex(InputError, "明确指定"):
                        generate(json.dumps(spec), "complex", Settings())


if __name__ == "__main__":
    unittest.main()
