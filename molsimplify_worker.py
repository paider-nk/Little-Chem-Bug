"""JSON-in/JSON-out worker. Run with the optional molSimplify interpreter."""
from __future__ import annotations

import importlib.metadata
import json
from pathlib import Path
import random
import sys
import traceback


def install_windows_compatibility():
    """Equivalent to 2.0.0's converter, using memory to avoid locked XYZ files."""
    from openbabel import openbabel
    from molSimplify.Classes.mol3D import mol3D

    def convert(self, force_clean=False, ignoreX=False):
        matrix = None
        if self.OBMol and not force_clean:
            matrix = self.populateBOMatrix()
        elif self.bo_mat.size != 0 and not force_clean:
            matrix = self.bo_mat
        conversion = openbabel.OBConversion()
        conversion.SetInFormat("xyz")
        molecule = openbabel.OBMol()
        if not conversion.ReadString(molecule, self.writexyz("memory", symbsonly=True, ignoreX=ignoreX, writestring=True)):
            raise ValueError("Open Babel 无法读取构建中的 XYZ 结构。")
        self.OBMol = molecule
        if matrix is not None and not force_clean:
            self.cleanBonds()
            for i in range(self.natoms):
                for j in range(self.natoms):
                    if matrix[i][j] > 0:
                        self.OBMol.AddBond(i + 1, j + 1, int(matrix[i][j]))

    mol3D.convert2OBMol = convert


def build(request):
    import numpy as np
    from molSimplify.Scripts.generator import startgen_pythonic
    from molSimplify.Scripts.io import getlicores
    from molSimplify.Scripts import structgen, distgeom
    from molSimplify.Classes.mol3D import mol3D
    from openbabel import openbabel

    install_windows_compatibility()
    # 2.0.0's default generator discards the assembled bond dictionary.
    # Capture each *placed* ligand's original graph and explicit donor atoms;
    # never infer coordination or ligand bond orders from XYZ distances.
    placements = []
    active = {}
    original_init = structgen.init_ligand
    original_combine = mol3D.roland_combine

    def init_ligand(*args, **kwargs):
        ligand, remove_pi, pi_atoms = original_init(*args, **kwargs)
        active.update(remove_pi=remove_pi, pi_atoms=pi_atoms)
        return ligand, remove_pi, pi_atoms

    def capture(self, ligand, donors, *args, **kwargs):
        count = ligand.natoms - int(active["remove_pi"])
        actual_donors = list(active["pi_atoms"] if active["remove_pi"] else donors)
        placements.append({"symbols": [a.symbol() for a in ligand.atoms[:count]],
                           "bonds": [(i, j, order) for (i, j), order in ligand.bo_dict.items() if i < count and j < count],
                           "donors": actual_donors})
        return original_combine(self, ligand, donors, *args, **kwargs)

    structgen.init_ligand = init_ligand
    mol3D.roland_combine = capture
    random.seed(20261117)
    np.random.seed(20261117)
    spec = request["spec"]
    library = getlicores()
    coordination = {"li": 2, "tpl": 3, "sqp": 4, "thd": 4, "spy": 5,
                    "tbp": 5, "oct": 6, "tpr": 6, "pbp": 7, "sqap": 8}[spec["geometry"]]
    sites = 0
    for ligand in spec["ligands"]:
        if ligand["kind"] == "library":
            if ligand["value"] not in library:
                raise ValueError(f"配体库中不存在 {ligand['value']}；可使用 ammonia、water、bipy、en、chloride、cyclopentadienyl，或选择自定义 SMILES。")
            cats = library[ligand["value"]][2]
            denticity = 1 if isinstance(cats, str) or "pi" in cats else len(cats)
        else:
            denticity = len(ligand["donors"])
        sites += denticity * ligand["count"]
    if sites != coordination:
        raise ValueError(f"配体共占用 {sites} 个配位位点，但所选构型要求 {coordination} 个；π 配体按一个方向位点计。")
    ligand_inputs = []
    for index, ligand in enumerate(spec["ligands"]):
        if ligand["kind"] == "smiles":
            # A controlled full path avoids upstream fuzzy library-name matching
            # (e.g. interpreting a valid SMILES as a similarly named library ligand).
            path = Path.cwd() / f"explicit_user_ligand_{index:02d}.mol"
            molecule = mol3D().getOBMol(ligand["value"], "smistring", True)
            conversion = openbabel.OBConversion()
            conversion.SetOutFormat("mol")
            if not conversion.WriteFile(molecule, str(path)):
                raise ValueError("无法准备自定义配体结构。")
            conversion.CloseOutFile()
            ligand_inputs.append(path.as_posix())
        else:
            ligand_inputs.append(ligand["value"])
    arguments = {
        "-core": spec["metal"].lower(), "-oxstate": str(spec["oxidation_state"]),
        "-geometry": spec["geometry"], "-coord": str(coordination),
        "-lig": ",".join(ligand_inputs),
        "-ligocc": ",".join(str(ligand["count"]) for ligand in spec["ligands"]),
        "-spin": str(request["multiplicity"]), "-charge": str(request["charge"]),
        "-skipANN": "True", "-ffoption": "N", "-keepHs": "yes",
        "-rundir": str(Path.cwd() / "Runs"),
    }
    custom = [ligand for ligand in spec["ligands"] if ligand["kind"] == "smiles"]
    if custom:
        arguments["-smicat"] = ",".join("[" + ",".join(map(str, ligand["donors"])) + "]" for ligand in custom)
    explicit_denticities = {name: len(ligand["donors"]) for name, ligand in zip(ligand_inputs, spec["ligands"]) if ligand["kind"] == "smiles"}
    original_shape_loader = distgeom.lig_load

    def load_for_shape(name, *args, **kwargs):
        ligand, error = original_shape_loader(name, *args, **kwargs)
        if ligand and name in explicit_denticities:
            ligand.denticity = explicit_denticities[name]
        return ligand, error

    # Upstream findshape otherwise assumes every file ligand is monodentate.
    distgeom.lig_load = load_for_shape
    try:
        result = startgen_pythonic(arguments, flag=True, write=False)
    finally:
        structgen.init_ligand = original_init
        mol3D.roland_combine = original_combine
        distgeom.lig_load = original_shape_loader
    if not isinstance(result, tuple) or len(result) != 3:
        raise ValueError(str(result))
    _, error, diagnostic = result
    if error:
        raise ValueError(str(error))
    if not diagnostic or not diagnostic.mol:
        raise ValueError("构建器未返回有效的 mol3D 结构。")
    mol = diagnostic.mol
    if diagnostic.sanity_is_set and diagnostic.sanity:
        raise ValueError("构建器检测到原子碰撞；请调整配体 / 构型。")
    atoms = [{"symbol": atom.symbol(), "xyz": [float(v) for v in atom.coords()]} for atom in mol.atoms]
    metals = [i for i, atom in enumerate(atoms) if atom["symbol"] == spec["metal"]]
    if len(metals) != 1:
        raise ValueError("当前仅支持单金属中心配合物。")
    bonds = []
    offset = 1
    expected_symbols = [spec["metal"]]
    if metals[0] != 0:
        raise ValueError("构建器改变了金属中心编号，无法映射配体连接。")
    for placed in placements:
        expected_symbols.extend(placed["symbols"])
        for i, j, value in placed["bonds"]:
            order = 1.5 if str(value).lower() in ("ar", "aro", "aromatic") else (1.0 if str(value).lower() == "am" else float(value))
            if order > 0:
                bonds.append([int(i) + offset, int(j) + offset, order])
        bonds.extend([0, int(donor) + offset, 1.0] for donor in placed["donors"])
        offset += len(placed["symbols"])
    if expected_symbols != [atom["symbol"] for atom in atoms]:
        raise ValueError("构建器改变了配体原子顺序或数量，已阻止不可靠的拓扑映射。")
    if not bonds:
        raise ValueError("构建器没有返回键连表，不能安全导出配合物。")
    return {"atoms": atoms, "bonds": bonds, "metal_index": metals[0],
            "version": importlib.metadata.version("molSimplify")}


if __name__ == "__main__":
    source, destination = map(Path, sys.argv[1:3])
    try:
        response = build(json.loads(source.read_text(encoding="utf-8")))
    except Exception as exc:
        response = {"error": str(exc)}
        traceback.print_exc()
    destination.write_text(json.dumps(response, ensure_ascii=False), encoding="utf-8")
    raise SystemExit(1 if "error" in response else 0)
