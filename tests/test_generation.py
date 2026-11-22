import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rdkit import Chem
from gaussian_core import (InputError, Settings, _geometry, _trigonal_boron_dihedrals,
                           build_polymer, default_checkpoint,
                           generate, resolve_molecule)
from rdkit.Chem import AllChem, rdMolDescriptors
from gaussian_input import cli


def atoms_from_input(text):
    sections = text.strip().split("\n\n")
    return [(line.split()[0], tuple(map(float, line.split()[1:])))
            for line in sections[2].splitlines()[1:]]


class PolymerTests(unittest.TestCase):
    def test_repeat_count_and_hydrogen_caps(self):
        for degree in (1, 2, 5, 100):
            with self.subTest(degree=degree):
                mol, _, _ = build_polymer("*CC*", degree)
                self.assertEqual(mol.GetNumAtoms(), 2 * degree)
                self.assertEqual(rdMolDescriptors.CalcMolFormula(mol), f"C{2 * degree}H{4 * degree + 2}")
                self.assertEqual(len(Chem.GetMolFrags(mol)), 1)
                self.assertFalse(any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()))

    def test_asymmetric_unit_orientation_and_methyl_caps(self):
        for cap, expected in (("hydrogen", "CCOCCOCCO"), ("methyl", "CCCOCCOCCOC")):
            mol, _, _ = build_polymer("[*:1]CCO[*:2]", 3, cap)
            self.assertEqual(Chem.MolToSmiles(mol), Chem.MolToSmiles(Chem.MolFromSmiles(expected)))
        mol, _, _ = build_polymer("*CCO*", 1, "methyl")
        self.assertEqual(rdMolDescriptors.CalcMolFormula(mol), "C4H10O")
        # Reverse both port labels: the capped chain is the same graph.
        reversed_mol, _, _ = build_polymer("[*:2]CCO[*:1]", 3, "methyl")
        forward, _, _ = build_polymer("*CCO*", 3, "methyl")
        self.assertEqual(Chem.MolToSmiles(reversed_mol), Chem.MolToSmiles(forward))

    def test_branch_charge_and_isotope_are_copied(self):
        branched, _, notes = build_polymer("*CC(C)*", 4)
        self.assertEqual(rdMolDescriptors.CalcMolFormula(branched), "C12H26")
        self.assertTrue(any("未指定" in w for w in notes))
        ionic, _, _ = build_polymer("*C[NH+](C)C*", 3)
        self.assertEqual(Chem.GetFormalCharge(ionic), 3)
        labeled, _, _ = build_polymer("*[13CH2]C*", 3)
        self.assertEqual(sum(a.GetIsotope() == 13 for a in labeled.GetAtoms()), 3)

    def test_port_stereochemistry_is_retained(self):
        for cap, expected in (("hydrogen", "[H][C@](F)(Cl)CO[C@](F)(Cl)CO[C@](F)(Cl)CO"),
                              ("methyl", "C[C@](F)(Cl)CO[C@](F)(Cl)CO[C@](F)(Cl)COC")):
            with self.subTest(cap=cap):
                mol, _, _ = build_polymer("[*:1][C@](F)(Cl)CO[*:2]", 3, cap)
                self.assertEqual(Chem.MolToSmiles(mol, isomericSmiles=True),
                                 Chem.MolToSmiles(Chem.MolFromSmiles(expected), isomericSmiles=True))

    def test_invalid_repeat_and_limits_are_rejected(self):
        invalid = ("CC", "*CC", "*C(*)C*", "*CC*.[Na+]", "C*C*", "*=CC*",
                   "[*:1]CC*", "[*:1]CC[*:1]", "[*:3]CC[*:4]", "[1*]CC*", "[H]*CC*", "*CC*\nC")
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(InputError):
                build_polymer(value, 2)
        for degree in (0, -1, 101, 2.5, True, "3"):
            with self.subTest(degree=degree), self.assertRaises(InputError):
                build_polymer("*CC*", degree)
        for value, degree, cap in (("*CCC*", 67, "hydrogen"), ("*CC*", 100, "methyl"),
                                    ("*CC*", 2, "unknown")):
            with self.subTest(value=value, degree=degree, cap=cap), self.assertRaises(InputError):
                build_polymer(value, degree, cap)

    def test_cli_polymer_geometry_and_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            output, metadata = Path(directory) / "poly.gjf", Path(directory) / "poly.json"
            args = ["*CC*", "--mode", "polymer", "--degree", "3", "--end-group", "methyl",
                    "--conformers", "2", "-o", str(output), "--metadata", str(metadata)]
            self.assertEqual(cli(args), 0)
            data = json.loads(metadata.read_text(encoding="utf-8"))
            self.assertEqual((data["formula"], data["degree"], data["repeat_unit"], data["end_group"]),
                             ("C8H18", 3, "*CC*", "methyl"))
            atoms = atoms_from_input(output.read_text(encoding="ascii"))
            self.assertEqual(len(atoms), 26)
            mol, _, _ = build_polymer("*CC*", 3, "methyl")
            for bond in Chem.AddHs(mol).GetBonds():
                distance = math.dist(atoms[bond.GetBeginAtomIdx()][1], atoms[bond.GetEndAtomIdx()][1])
                self.assertTrue(0.8 < distance < 1.7, distance)
            self.assertEqual(cli(["CC", "--degree", "3", "-o", str(output), "--force"]), 2)
        checkpoints = {default_checkpoint("*CC*", "polymer", n, cap)
                       for n in (2, 3) for cap in ("hydrogen", "methyl")}
        self.assertEqual(len(checkpoints), 4)


class ChemicalGenerationTests(unittest.TestCase):
    def test_ethanol_format_geometry_and_bond_lengths(self):
        result = generate("CCO", "smiles", Settings(conformers=4))
        self.assertEqual((result.formula, result.charge, result.multiplicity, result.atom_count),
                         ("C2H6O", 0, 1, 9))
        self.assertTrue(result.text.startswith("%chk=molecule.chk\n%mem=4GB\n%nprocshared=4\n#p "))
        self.assertIn("\n\nMolecule\n\n0 1\n", result.text)
        self.assertTrue(result.text.endswith("\n\n"))
        self.assertTrue(result.text.isascii())
        atoms = atoms_from_input(result.text)
        self.assertEqual(len(atoms), 9)
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            length = math.dist(atoms[i][1], atoms[j][1])
            self.assertTrue(0.8 < length < 1.7, (i, j, length))
        self.assertIsNotNone(result.energy)

    def test_charge_radical_and_high_spin(self):
        for smiles, charge, mult in (("CC(=O)[O-]", -1, 1), ("[NH4+]", 1, 1),
                                     ("[CH3]", 0, 2), ("[O][O]", 0, 3), ("[H]", 0, 2)):
            with self.subTest(smiles=smiles):
                result = generate(smiles, "smiles", Settings(conformers=2))
                self.assertEqual((result.charge, result.multiplicity), (charge, mult))

    def test_electron_spin_parity_rejected(self):
        with self.assertRaisesRegex(InputError, "不相容"):
            generate("CCO", settings=Settings(multiplicity=2))

    def test_manual_charge_changes_electrons_without_hydrogens(self):
        result = generate("CCO", settings=Settings(charge=1, multiplicity=2, conformers=1))
        self.assertEqual((result.atom_count, result.formula, result.charge), (9, "C2H6O", 1))
        self.assertTrue(any("不改变氢数" in w for w in result.warnings))

    def test_stereochemistry_survives_coordinates(self):
        for smiles in ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"):
            with self.subTest(smiles=smiles):
                result = generate(smiles, "smiles", Settings(conformers=3))
                expected = Chem.MolFromSmiles(smiles)
                expected_tag = Chem.FindMolChiralCenters(expected, includeUnassigned=True)
                probe = Chem.AddHs(expected)
                conf = Chem.Conformer(probe.GetNumAtoms())
                conf.Set3D(True)
                for i, (_, xyz) in enumerate(atoms_from_input(result.text)):
                    conf.SetAtomPosition(i, xyz)
                probe.AddConformer(conf, assignId=True)
                Chem.RemoveStereochemistry(probe)
                Chem.AssignAtomChiralTagsFromStructure(probe)
                self.assertEqual(Chem.FindMolChiralCenters(Chem.RemoveHs(probe), includeUnassigned=True), expected_tag)

    def test_unspecified_stereo_is_reported(self):
        result = generate("NC(C)C(=O)O", settings=Settings(conformers=1))
        self.assertTrue(any("未指定" in w for w in result.warnings))

    def test_isotope_mass_is_written(self):
        result = generate("[13CH4]", settings=Settings(conformers=1))
        self.assertIn("C(Iso=13)", result.text)

    def test_offline_chinese_name(self):
        with patch("gaussian_core._pubchem_name", side_effect=AssertionError("unexpected network")):
            result = generate("乙醇", "name", Settings(conformers=1))
        self.assertEqual((result.smiles, result.source), ("CCO", "内置常见名称"))

    def test_pubchem_response_preserves_isomeric_smiles(self):
        payload = {"PropertyTable": {"Properties": [{"CID": 5950, "SMILES": "N[C@@H](C)C(=O)O"}]}}
        with patch("gaussian_core.urllib.request.urlopen") as opened:
            opened.return_value.__enter__.return_value.read.return_value = json.dumps(payload).encode()
            result = generate("L-alanine", "name", Settings(conformers=1))
        self.assertIn("@", result.smiles)
        self.assertEqual(result.source, "PubChem CID 5950")

    def test_pubchem_ambiguous_name_is_rejected(self):
        payload = {"PropertyTable": {"Properties": [{"CID": 1, "SMILES": "C"}, {"CID": 2, "SMILES": "CC"}]}}
        with patch("gaussian_core.urllib.request.urlopen") as opened:
            opened.return_value.__enter__.return_value.read.return_value = json.dumps(payload).encode()
            with self.assertRaisesRegex(InputError, "多个结构"):
                generate("ambiguous", "name")

    def test_network_failure_is_readable(self):
        import urllib.error
        with patch("gaussian_core.urllib.request.urlopen", side_effect=urllib.error.URLError("offline")):
            with self.assertRaisesRegex(InputError, "无法连接"):
                generate("a remote chemical name", "name")

    def test_invalid_or_unsupported_inputs(self):
        for value in ("", "not a SMILES", "CCO extra text", "C(C", "*", "c1cccc1", "CCO\n%mem=1GB"):
            with self.subTest(value=value), self.assertRaises(InputError):
                generate(value, "smiles")

    def test_metals_require_manual_spin(self):
        with self.assertRaisesRegex(InputError, "手动填写"):
            generate("[Fe+2]", "smiles")

    def test_coordination_smiles_metal_complex_with_manual_spin(self):
        result = generate("N->[Cu+2](<-N)(<-N)<-N", "smiles",
                          Settings(multiplicity=2, conformers=1))
        self.assertEqual(result.formula, "H12CuN4+2")
        self.assertEqual((result.charge, result.multiplicity), (2, 2))
        self.assertIn("Cu", result.text)
        self.assertTrue(any("配位 SMILES" in warning for warning in result.warnings))

    def test_rdkit_dot_disconnected_metal_and_ligands_are_accepted_without_overlap(self):
        result = generate("[Fe+2].[cH-]1cccc1.[cH-]1cccc1", "smiles",
                          Settings(multiplicity=1, conformers=1))
        self.assertEqual(result.charge, 0)
        self.assertEqual(result.formula, "C10H10Fe")
        self.assertTrue(any("多个独立片段" in warning for warning in result.warnings))
        rows = []
        started = False
        for line in result.text.splitlines():
            if line.strip() == "0 1":
                started = True
                continue
            if started:
                if not line.strip():
                    break
                fields = line.split()
                rows.append((fields[0], np.array(list(map(float, fields[1:4])))))
        iron = next(point for symbol, point in rows if symbol == "Fe")
        ring_atoms = [point for symbol, point in rows if symbol in ("C", "H")]
        self.assertGreater(min(np.linalg.norm(point - iron) for point in ring_atoms), 2.5)

    def test_ferrocene_template_adds_fe_carbon_connectivity_for_gaussian(self):
        result = generate("ferrocene", "ferrocene", Settings(conformers=1))
        route = next(line for line in result.text.splitlines() if line.startswith("#"))
        self.assertIn("Geom=Connectivity", route)
        self.assertEqual((result.formula, result.charge, result.multiplicity), ("C10H10Fe", 0, 1))
        self.assertTrue(any("10 个碳" in warning for warning in result.warnings))
        lines = result.text.splitlines()
        charge_line = lines.index("0 1")
        coordinate_end = next(i for i in range(charge_line + 1, len(lines)) if not lines[i].strip())
        coordinates = [line.split() for line in lines[charge_line + 1:coordinate_end]]
        connectivity = [line.split() for line in lines[coordinate_end + 1:] if line.strip()]
        iron = next(i for i, row in enumerate(coordinates, 1) if row[0] == "Fe")
        fe_neighbors = set()
        for row in connectivity:
            atom = int(row[0])
            for neighbor, order in zip(row[1::2], row[2::2]):
                if atom == iron and coordinates[int(neighbor) - 1][0] == "C":
                    fe_neighbors.add(int(neighbor))
                if int(neighbor) == iron and coordinates[atom - 1][0] == "C":
                    fe_neighbors.add(atom)
        self.assertEqual(len(fe_neighbors), 10)
        iron_xyz = np.array(list(map(float, coordinates[iron - 1][1:4])))
        for idx in fe_neighbors:
            atom_xyz = np.array(list(map(float, coordinates[idx - 1][1:4])))
            self.assertAlmostEqual(np.linalg.norm(atom_xyz - iron_xyz), 2.02, delta=0.05)

    def test_trigonal_boron_and_three_neighbors_are_planar_and_dihedral_constrained(self):
        result = generate("c1ccccc1B(O)O", "smiles", Settings(conformers=1))
        route = next(line for line in result.text.splitlines() if line.startswith("#"))
        self.assertIn("Opt=(ModRedundant)", route)
        constraints = [line for line in result.text.splitlines() if line.startswith("D ")]
        self.assertEqual(len(constraints), 4)
        self.assertTrue(any(line.endswith("180.0") for line in constraints))
        for value_line, freeze_line in zip(constraints[::2], constraints[1::2]):
            self.assertEqual(value_line.rsplit(" ", 1)[0] + " F", freeze_line)
        coordinates = {}
        reading = False
        for line in result.text.splitlines():
            if line.strip() == "0 1":
                reading = True
                continue
            if reading:
                if not line.strip():
                    break
                fields = line.split()
                coordinates[len(coordinates)] = tuple(map(float, fields[1:4]))
        atoms = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1B(O)O"))
        b = next(a.GetIdx() for a in atoms.GetAtoms() if a.GetAtomicNum() == 5)
        neighbors = [a.GetIdx() for a in atoms.GetAtomWithIdx(b).GetNeighbors()]
        origin = np.array(coordinates[b])
        vectors = np.array([np.array(coordinates[i]) - origin for i in neighbors])
        self.assertLess(abs(np.linalg.det(vectors)), 1e-4)
        self.assertTrue(any("180°" in warning for warning in result.warnings))

    def test_bpin_boronate_keeps_closed_ring_bonds_and_uses_dihedral_constraints(self):
        smiles = "CC1(C)OB(c2coc3ccccc23)OC1(C)C"
        result = generate(smiles, "smiles", Settings(conformers=1))
        self.assertEqual(result.formula, "C14H17BO3")
        self.assertEqual(len([line for line in result.text.splitlines() if line.startswith("D ")]), 4)
        self.assertNotIn("\nA ", result.text)
        atoms = Chem.AddHs(Chem.MolFromSmiles(smiles))
        coordinates = []
        reading = False
        for line in result.text.splitlines():
            if line.strip() == "0 1":
                reading = True
                continue
            if reading:
                if not line.strip():
                    break
                fields = line.split()
                coordinates.append(np.array(list(map(float, fields[1:4]))))
        self.assertEqual(len(coordinates), atoms.GetNumAtoms())
        for bond in atoms.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            distance = np.linalg.norm(coordinates[i] - coordinates[j])
            self.assertGreater(distance, 0.8, f"bond {i + 1}-{j + 1} collapsed")
            self.assertLess(distance, 1.9, f"bond {i + 1}-{j + 1} is disconnected")

    def test_boron_planarization_preserves_all_bond_lengths_in_attached_groups(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1B(O)O"))
        cid, *_ = _geometry(mol, 1)
        conf = mol.GetConformer(cid)
        AllChem.UFFOptimizeMolecule(mol, confId=cid)
        before = {}
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            before[(i, j)] = np.linalg.norm(np.asarray(conf.GetAtomPosition(i)) -
                                            np.asarray(conf.GetAtomPosition(j)))
        constraints = _trigonal_boron_dihedrals(mol, conf)
        self.assertEqual(len(constraints), 2)
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            after = np.linalg.norm(np.asarray(conf.GetAtomPosition(i)) -
                                   np.asarray(conf.GetAtomPosition(j)))
            self.assertAlmostEqual(before[(i, j)], after, places=7)

    def test_boron_sp_does_not_add_optimization_constraint(self):
        result = generate("CCB(O)O", "smiles", Settings(job="sp", conformers=1))
        self.assertNotIn("ModRedundant", result.text)
        self.assertFalse(any(line.startswith("D ") for line in result.text.splitlines()))

    def test_special_routes_and_header_injection_rejected(self):
        for settings in (Settings(basis="GenECP"), Settings(extra_keywords="Geom=AllCheck"),
                         Settings(extra_keywords="Opt=ModRedundant"), Settings(extra_keywords="Guess=Read"),
                         Settings(extra_keywords="\n%chk=bad"), Settings(checkpoint="../bad"),
                         Settings(title="bad\n\n0 1"), Settings(conformers=0), Settings(memory_gb=-1)):
            with self.subTest(settings=settings), self.assertRaises(InputError):
                generate("O", "smiles", settings)

    def test_cli_writes_gjf_metadata_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "water.gjf"
            metadata = Path(directory) / "water.json"
            args = ["O", "--mode", "smiles", "-o", str(output), "--metadata", str(metadata), "--conformers", "1"]
            self.assertEqual(cli(args), 0)
            self.assertIn("\n0 1\n", output.read_text(encoding="ascii"))
            self.assertEqual(json.loads(metadata.read_text(encoding="utf-8"))["formula"], "H2O")
            self.assertEqual(cli(args), 2)
            self.assertEqual(cli(args + ["--force"]), 0)

    def test_deterministic_initial_coordinates(self):
        settings = Settings(conformers=2)
        self.assertEqual(generate("CCCO", settings=settings).text, generate("CCCO", settings=settings).text)


if __name__ == "__main__":
    unittest.main()
