"""Optional molSimplify structure builder, isolated from the RDKit runtime."""
from __future__ import annotations

import importlib.util
import json
import math
import os
from pathlib import Path
import re
import external_process as subprocess
import sys
import tempfile

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from gaussian_core import (InputError, Result, Settings, _trigonal_boron_dihedrals,
                           _validate_settings, _external_preprocessing)


GEOMETRIES = {
    "li": ("直线形", 2), "tpl": ("三角平面", 3),
    "sqp": ("平面四方", 4), "thd": ("四面体", 4),
    "spy": ("四方锥", 5), "tbp": ("三角双锥", 5),
    "oct": ("八面体", 6), "tpr": ("三棱柱", 6),
    "pbp": ("五角双锥", 7), "sqap": ("四方反棱柱", 8),
}

PRESETS = {
    "四氨合铜(II)": ({"metal": "Cu", "oxidation_state": 2, "geometry": "sqp",
                      "ligands": [{"kind": "library", "value": "ammonia", "count": 4}]}, 2, 2),
    "六水合铁(II)": ({"metal": "Fe", "oxidation_state": 2, "geometry": "oct",
                      "ligands": [{"kind": "library", "value": "water", "count": 6}]}, 2, 5),
    "三联吡啶合铁(II)": ({"metal": "Fe", "oxidation_state": 2, "geometry": "oct",
                         "ligands": [{"kind": "library", "value": "bipy", "count": 3}]}, 2, 1),
    "二茂铁（π 配位）": ({"metal": "Fe", "oxidation_state": 2, "geometry": "li",
                          "ligands": [{"kind": "library", "value": "cyclopentadienyl", "count": 2}]}, 0, 1),
}


def validate_spec(spec):
    if not isinstance(spec, dict):
        raise InputError("配合物配置必须是 JSON 对象。")
    metal = spec.get("metal")
    if not isinstance(metal, str) or not re.fullmatch(r"[A-Z][a-z]?", metal):
        raise InputError("金属中心请输入元素符号，例如 Cu、Fe、Ni、Pt。")
    try:
        z = Chem.GetPeriodicTable().GetAtomicNumber(metal)
    except RuntimeError as exc:
        raise InputError("无法识别金属元素符号。") from exc
    if not (21 <= z <= 30 or 39 <= z <= 48 or 57 <= z <= 80 or 89 <= z <= 112):
        raise InputError("当前配合物后端支持单个过渡金属或 f 区金属中心。")
    oxidation = spec.get("oxidation_state")
    if type(oxidation) is not int or not -4 <= oxidation <= 8:
        raise InputError("金属氧化态请输入 -4 到 8 的整数；它不等于配合物总电荷。")
    geometry = spec.get("geometry")
    if geometry not in GEOMETRIES:
        raise InputError("请选择支持的配位构型。")
    ligands = spec.get("ligands")
    if not isinstance(ligands, list) or not 1 <= len(ligands) <= 8:
        raise InputError("请填写 1–8 种配体。")
    cleaned = []
    for ligand in ligands:
        if not isinstance(ligand, dict):
            raise InputError("每种配体必须包含类型、名称或 SMILES、数量。")
        kind, value, count = ligand.get("kind"), ligand.get("value"), ligand.get("count")
        if kind not in ("library", "smiles") or not isinstance(value, str) or not value.strip():
            raise InputError("配体类型请选择 library（配体库）或 smiles，并填写结构。")
        value = value.strip()
        if len(value) > 2000 or any(char.isspace() or ord(char) < 32 for char in value):
            raise InputError("配体名称 / SMILES 不能含空白或换行。")
        if type(count) is not int or not 1 <= count <= 8:
            raise InputError("每种配体的数量请输入 1–8 的整数。")
        item = {"kind": kind, "value": value, "count": count}
        if kind == "library":
            if not re.fullmatch(r"[A-Za-z0-9_+\-]+", value):
                raise InputError("配体库名称请使用库内英文标识，例如 ammonia、water、bipy。")
        else:
            params = Chem.SmilesParserParams()
            params.parseName = False
            mol = Chem.MolFromSmiles(value, params)
            if mol is None or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
                raise InputError("配体必须是 RDKit 可解析的完整、连通 SMILES。")
            if mol.GetNumHeavyAtoms() * count > 200:
                raise InputError("配体规模超过 200 个非氢原子限制。")
            donors = ligand.get("donors")
            if (not isinstance(donors, list) or not 1 <= len(donors) <= 8 or
                    any(type(i) is not int or not 1 <= i <= mol.GetNumAtoms() for i in donors) or
                    len(set(donors)) != len(donors)):
                raise InputError("SMILES 配体必须指定有效且不重复的配位原子编号（从 1 开始，按 SMILES 原子顺序）。")
            if any(mol.GetAtomWithIdx(i - 1).GetAtomicNum() == 1 for i in donors):
                raise InputError("当前自定义配体的配位原子不能是氢。")
            item["donors"] = donors
        cleaned.append(item)
    if sum(ligand["count"] for ligand in cleaned) > 8:
        raise InputError("当前单中心模型最多容纳 8 个配体。")
    return {"metal": metal, "oxidation_state": oxidation, "geometry": geometry, "ligands": cleaned}


def load_spec(value):
    if not isinstance(value, str) or len(value) > 50000:
        raise InputError("配合物配置过长或类型无效。")
    try:
        if value.lstrip().startswith("{"):
            return validate_spec(json.loads(value))
        path = Path(value)
        if path.suffix.lower() != ".json" or not path.is_file() or path.stat().st_size > 50000:
            raise InputError("配合物模式请输入配置 JSON 文件路径，或使用“配置配合物”按钮。")
        return validate_spec(json.loads(path.read_text(encoding="utf-8-sig")))
    except (ValueError, OSError) as exc:
        if isinstance(exc, InputError):
            raise
        raise InputError(f"无法读取配合物配置：{exc}") from exc


def backend_python():
    configured = os.environ.get("MOLSIMPLIFY_PYTHON")
    if configured:
        path = Path(configured)
        if not path.is_file():
            raise InputError("MOLSIMPLIFY_PYTHON 未指向有效的 Python 可执行文件。")
        return str(path.resolve())
    base = Path(__file__).resolve().parent / ".molsimplify"
    path = base / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    if path.is_file():
        return str(path)
    if not getattr(sys, "frozen", False) and importlib.util.find_spec("molSimplify") is not None:
        return sys.executable
    raise InputError("尚未安装 molSimplify 后端。请在「工具与安装」页安装或配置 molSimplify。")


def build_geometry(spec, charge, multiplicity):
    interpreter = backend_python()
    with tempfile.TemporaryDirectory(prefix="gaussian_complex_") as folder:
        folder = Path(folder)
        request, response = folder / "request.json", folder / "response.json"
        request.write_text(json.dumps({"spec": spec, "charge": charge, "multiplicity": multiplicity}), encoding="utf-8")
        environment = os.environ.copy()
        for key in ("PYTHONHOME", "PYTHONPATH", "VIRTUAL_ENV"):
            environment.pop(key, None)
        environment.update(PYTHONUTF8="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", TF_CPP_MIN_LOG_LEVEL="3")
        try:
            completed = subprocess.run(
                [interpreter, str(Path(__file__).with_name("molsimplify_worker.py")), str(request), str(response)],
                cwd=folder, env=environment, capture_output=True, text=True, encoding="utf-8", errors="replace",
                timeout=240, creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        except subprocess.TimeoutExpired as exc:
            raise InputError("molSimplify 构建超过 240 秒；请减少配体复杂度后重试。") from exc
        except OSError as exc:
            raise InputError(f"无法启动 molSimplify：{exc}") from exc
        if not response.is_file():
            raise InputError("molSimplify 未返回结构：" + (completed.stderr or completed.stdout)[-2000:])
        payload = json.loads(response.read_text(encoding="utf-8"))
        if completed.returncode or "error" in payload:
            raise InputError("molSimplify 构建失败：" + payload.get("error", completed.stderr[-1500:]))
        return payload


def generate_complex(value: str, settings: Settings, *, on_stage=None):
    stage = on_stage if on_stage is not None else lambda _: None
    _validate_settings(settings)
    spec = load_spec(value)
    if settings.charge is None or settings.multiplicity is None:
        raise InputError("配合物模式请明确指定总电荷和自旋多重度；金属氧化态不能代替总电荷。")
    stage("委托 molSimplify 组装配合物")
    payload = build_geometry(spec, settings.charge, settings.multiplicity)
    atoms, bonds = payload["atoms"], payload["bonds"]
    if not atoms or len(atoms) > 800 or sum(a["symbol"] != "H" for a in atoms) > 200:
        raise InputError("构建结果超过 200 个非氢原子 / 800 个总原子限制。")
    mol = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))
    for index, atom in enumerate(atoms):
        new_atom = Chem.Atom(atom["symbol"])
        new_atom.SetNoImplicit(True)
        mol.AddAtom(new_atom)
        if not all(math.isfinite(v) for v in atom["xyz"]):
            raise InputError("molSimplify 返回了无效三维坐标。")
        conf.SetAtomPosition(index, atom["xyz"])
    bond_types = {1.0: Chem.BondType.SINGLE, 1.5: Chem.BondType.AROMATIC,
                  2.0: Chem.BondType.DOUBLE, 3.0: Chem.BondType.TRIPLE}
    metal = payload["metal_index"]
    neighbors = []
    for i, j, order in bonds:
        if metal in (i, j):
            donor = j if i == metal else i
            neighbors.append(donor)
            mol.AddBond(donor, metal, Chem.BondType.DATIVE)
        else:
            mol.AddBond(i, j, bond_types.get(float(order), Chem.BondType.SINGLE))
    if len(Chem.GetMolFrags(mol)) != 1:
        raise InputError("molSimplify 返回的拓扑有未连接的配体，已停止生成。")
    mol.AddConformer(conf)
    mol.UpdatePropertyCache(strict=False)
    electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - settings.charge
    if electrons < 1 or settings.multiplicity > electrons + 1 or (electrons - settings.multiplicity + 1) % 2:
        raise InputError(f"配合物电荷 / 多重度不相容：{electrons} 个电子，多重度 {settings.multiplicity}。")
    boron_notes = []
    geometry_method = "molSimplify 配体定向组装（保持构建坐标）"
    if settings.preprocess_method in ("xtb", "tblite", "gfnff"):
        stage("运行配合物几何预处理")
        label, notes = _external_preprocessing(mol, mol.GetConformer().GetId(), settings.preprocess_method,
                                               settings.charge, settings.multiplicity, settings.processors)
        geometry_method = "molSimplify 配体定向组装 + " + label
        boron_notes.extend(notes)
    stage("检查配位几何并组装输入文件")
    constraints = _trigonal_boron_dihedrals(mol, mol.GetConformer(), settings.boron_preprocess,
                                          settings.boron_bond_tolerance, boron_notes)
    job = {"optfreq": "Opt Freq", "opt": "Opt", "sp": "SP", "freq": "Freq"}[settings.job]
    if constraints and settings.job in ("opt", "optfreq"):
        job = "Opt=(ModRedundant)" + (" Freq" if settings.job == "optfreq" else "")
    route = f"#p {settings.method}/{settings.basis} {job} Geom=Connectivity"
    if settings.extra_keywords.strip():
        route += " " + settings.extra_keywords.strip()
    lines = [f"%chk={settings.checkpoint}.chk", f"%mem={settings.memory_gb}GB",
             f"%nprocshared={settings.processors}", route, "", settings.title.strip(), "",
             f"{settings.charge} {settings.multiplicity}"]
    for atom in mol.GetAtoms():
        xyz = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<12} {xyz.x: .8f} {xyz.y: .8f} {xyz.z: .8f}")
    lines.append("")
    connectivity = {i: [] for i in range(len(atoms))}
    for i, j, order in bonds:
        i, j = sorted((i, j))
        # A fractional order records coordination in Gaussian/GaussView's
        # connectivity, without constraining the quantum optimization geometry.
        connectivity[i].append((j, 0.1 if metal in (i, j) else order))
    for i in range(len(atoms)):
        lines.append(str(i + 1) + "".join(f" {j + 1} {order:.1f}" for j, order in sorted(connectivity[i])))
    if constraints and settings.job in ("opt", "optfreq"):
        lines.append("")
        for i, j, k, l, target in constraints:
            coordinate = f"D {i + 1} {j + 1} {k + 1} {l + 1}"
            lines.extend((f"{coordinate} {target:.1f}", f"{coordinate} F"))
    # Formula uses explicit atoms; charge is reported separately, not inferred
    # from RDKit valence assignments to the metal centre.
    formula_mol = Chem.Mol(mol)
    formula = rdMolDescriptors.CalcMolFormula(formula_mol)
    warnings = [
        f"molSimplify {payload['version']}：{spec['metal']}，氧化态 {spec['oxidation_state']:+d}，"
        f"{GEOMETRIES[spec['geometry']][0]}，{len(neighbors)} 条金属–配位原子连接。",
        "按指定配体和构型组装初始结构；未启用神经网络键长 / 自旋预测，不自动确定基态或枚举 cis/trans 等异构体。",
        "金属–配体连接以 0.1 阶写入 Gaussian Connectivity；该数值记录拓扑，不是物理键级或键长冻结约束。",
        "总电荷和多重度由用户指定；请核对质子化态、氧化态、基组 / 赝势及最终优化结构。",
    ]
    warnings.extend(boron_notes)
    result = Result("\n".join(lines) + "\n\n", Chem.MolToSmiles(mol), formula,
                    settings.charge, settings.multiplicity, len(atoms),
                    geometry_method, None, "molSimplify " + payload["version"], warnings)
    result.complex_spec = spec
    return result
