"""SMILES / PubChem name -> reproducible initial geometry -> Gaussian input."""
from __future__ import annotations

import hashlib
import json
import math
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
import os
from pathlib import Path
import shutil
import external_process as subprocess
import tempfile

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors


class InputError(ValueError):
    """An error that can be shown directly to the user."""


ALIASES = {
    "水": "O", "water": "O", "甲烷": "C", "methane": "C",
    "乙烷": "CC", "乙烯": "C=C", "乙炔": "C#C",
    "甲醇": "CO", "methanol": "CO", "乙醇": "CCO", "ethanol": "CCO",
    "丙酮": "CC(=O)C", "acetone": "CC(=O)C",
    "苯": "c1ccccc1", "benzene": "c1ccccc1", "甲苯": "Cc1ccccc1",
    "甲醛": "C=O", "乙醛": "CC=O", "乙酸": "CC(=O)O", "醋酸": "CC(=O)O",
    "乙酸乙酯": "CCOC(C)=O", "乙腈": "CC#N", "二氯甲烷": "ClCCl",
    "三氯甲烷": "ClC(Cl)Cl", "氯仿": "ClC(Cl)Cl", "四氯化碳": "ClC(Cl)(Cl)Cl",
    "四氢呋喃": "C1CCOC1", "二甲基亚砜": "CS(C)=O",
    "二甲基甲酰胺": "CN(C)C=O", "正己烷": "CCCCCC", "环己烷": "C1CCCCC1",
    "苯酚": "Oc1ccccc1", "苯胺": "Nc1ccccc1", "吡啶": "c1ccncc1",
    "氨": "N", "氨气": "N", "ammonia": "N", "氮气": "N#N",
    "二氧化碳": "O=C=O", "氢气": "[H][H]", "铵根": "[NH4+]",
    "乙酸根": "CC(=O)[O-]", "甲基自由基": "[CH3]",
}


@dataclass
class Settings:
    method: str = "B3LYP"
    basis: str = "6-31G(d)"
    job: str = "optfreq"
    memory_gb: int = 4
    processors: int = 4
    charge: int | None = None
    multiplicity: int | None = None
    conformers: int = 20
    extra_keywords: str = "SCF=XQC Int=UltraFine"
    title: str = "Molecule"
    checkpoint: str = "molecule"
    degree: int = 5
    end_group: str = "hydrogen"
    boron_preprocess: str = "strict"
    boron_bond_tolerance: float = 0.15
    # Initial geometry optimizer. External xTB methods are opt-in and use
    # XTB_EXE or PATH; the default keeps the original RDKit-only behavior.
    preprocess_method: str = "auto"
    engine: str = "gaussian"
    cif_block: str = ""
    cif_component: int = 1
    cif_scope: str = "molecule"


@dataclass
class Result:
    text: str
    smiles: str
    formula: str
    charge: int
    multiplicity: int
    atom_count: int
    geometry_method: str
    energy: float | None
    source: str
    warnings: list[str] = field(default_factory=list)
    repeat_unit: str | None = None
    degree: int | None = None
    end_group: str | None = None
    complex_spec: dict | None = None


def _parse_smiles(value: str):
    # Do not leak RDKit parser diagnostics into the desktop console.
    with rdBase.BlockLogs():
        params = Chem.SmilesParserParams()
        params.parseName = False
        return Chem.MolFromSmiles(value, params)


def _pubchem_name(name: str) -> tuple[str, str, list[str]]:
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name, safe="") + "/property/IsomericSMILES/JSON")
    request = urllib.request.Request(url, headers={"User-Agent": "GaussianInputTool/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=15) as response:
            body = response.read(2_000_001)
        if len(body) > 2_000_000:
            raise InputError("名称查询返回内容过大，请改用明确的 SMILES。")
        properties = json.loads(body)["PropertyTable"]["Properties"]
    except urllib.error.HTTPError as exc:
        if exc.code in (400, 404):
            raise InputError("PubChem 未找到该名称。请使用英文名称、IUPAC 名称或 SMILES。") from exc
        raise InputError(f"PubChem 查询失败（HTTP {exc.code}），请稍后重试或改用 SMILES。") from exc
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise InputError("无法连接 PubChem。请检查网络，或使用 SMILES / 内置常见名称。") from exc
    except (KeyError, ValueError, TypeError) as exc:
        raise InputError("PubChem 返回格式异常，请改用 SMILES。") from exc
    if not properties:
        raise InputError("未找到对应分子，请改用 SMILES。")
    if len(properties) != 1:
        raise InputError("该名称对应多个结构，请使用明确的立体化学名称或 SMILES。")
    prop = properties[0]
    smiles = prop.get("SMILES") or prop.get("IsomericSMILES")
    if not isinstance(smiles, str) or not smiles:
        raise InputError("PubChem 未返回含立体化学信息的 SMILES，请直接输入 SMILES。")
    return smiles, f"PubChem CID {prop.get('CID', '未知')}", [
        "名称由 PubChem 解析，请核对返回的 SMILES 和所需异构体、质子化状态。"
    ]


def resolve_molecule(value: str, mode: str = "auto"):
    value = value.strip()
    if not value or len(value) > 4000 or any(c in value for c in "\r\n\x00"):
        raise InputError("请输入单行 SMILES 或分子名称（不超过 4000 字符）。")
    if mode not in ("auto", "smiles", "name"):
        raise InputError("输入类型必须是 auto、smiles 或 name。")
    warnings: list[str] = []
    mol = _parse_smiles(value) if mode != "name" else None
    source = "SMILES"
    if mol is None:
        if mode == "smiles":
            raise InputError("SMILES 无法解析，请检查括号、环编号、芳香性和原子价态。")
        alias = ALIASES.get(value.casefold())
        if alias:
            mol = _parse_smiles(alias)
            source = "内置常见名称"
        else:
            resolved, source, warnings = _pubchem_name(value)
            mol = _parse_smiles(resolved)
    return _check_structure(mol, source, warnings)


def _check_structure(mol, source: str, warnings: list[str]):
    if mol is None or mol.GetNumAtoms() == 0:
        raise InputError("无法得到有效分子结构。")
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise InputError("不支持 *、R 基或未指定元素，请输入完整分子。")
    if len(Chem.GetMolFrags(mol)) > 1:
        warnings.append("输入含点号分隔的多个独立片段；程序会将片段分开放置以避免坐标重叠。SMILES 未定义片段间相对位置或键合关系，输出只能作为片段集合的起始几何，不能据此视作已组装的配合物。")
    if mol.GetNumAtoms() > 200:
        raise InputError("当前小工具限制为 200 个补氢前原子；大分子请使用专用建模软件。")
    stereo = Chem.FindPotentialStereo(mol)
    if any(str(item.specified) == "Unspecified" for item in stereo):
        warnings.append("结构含未指定的立体中心或双键构型；当前三维结构只代表一种可能构型。")
    return mol, source, warnings


def build_polymer(repeat_smiles: str, degree: int, end_group: str = "hydrogen"):
    """Join oriented two-port repeat units using molzip's stereo-aware graph edits.

    Atom maps 1/2, when supplied, define head/tail; otherwise use SMILES parse
    order of the two wildcards. Both ports must be terminal single bonds.
    """
    if type(degree) is not int or not 1 <= degree <= 100:
        raise InputError("聚合度必须是 1–100 的整数，且最终链不能超过 200 个补氢前原子。")
    if end_group not in ("hydrogen", "methyl"):
        raise InputError("封端方式必须是 hydrogen（氢）或 methyl（甲基）。")
    value = repeat_smiles.strip()
    if not value or len(value) > 4000 or any(c in value for c in "\r\n\x00"):
        raise InputError("请输入单行链节 SMILES，以两个 * 标记连接位置。")
    unit = _parse_smiles(value)
    if unit is None or len(Chem.GetMolFrags(unit)) != 1:
        raise InputError("链节必须是可解析的单个连通 SMILES，例如 *CC* 或 [*:1]CCO[*:2]。")
    ports = [a for a in unit.GetAtoms() if a.GetAtomicNum() == 0]
    if len(ports) != 2:
        raise InputError("链节必须恰好含两个连接标记 *；例如聚乙烯用 *CC*，不是乙烯单体 C=C。")
    for atom in ports:
        if (atom.GetDegree() != 1 or atom.GetFormalCharge() or atom.GetIsotope() or
                atom.GetNumExplicitHs() or atom.GetNumRadicalElectrons() or
                atom.GetBonds()[0].GetBondType() != Chem.BondType.SINGLE or
                atom.GetNeighbors()[0].GetAtomicNum() <= 1):
            raise InputError("每个 * 必须中性、无同位素、仅通过一条单键连接到真实非氢原子。")
    maps = [a.GetAtomMapNum() for a in ports]
    if maps == [0, 0]:
        head, tail = [a.GetIdx() for a in ports]
    elif set(maps) == {1, 2}:
        head = next(a.GetIdx() for a in ports if a.GetAtomMapNum() == 1)
        tail = next(a.GetIdx() for a in ports if a.GetAtomMapNum() == 2)
    else:
        raise InputError("连接标记请都用 *，或成对使用 [*:1]（头）和 [*:2]（尾）。")
    estimated_atoms = degree * (unit.GetNumAtoms() - 2) + (2 if end_group == "methyl" else 0)
    if estimated_atoms > 200:
        raise InputError(f"该聚合度将产生约 {estimated_atoms} 个补氢前原子，超过当前 200 原子限制；请降低聚合度。")
    combined = None
    for i in range(degree):
        part = Chem.Mol(unit)
        for atom in part.GetAtoms():
            atom.SetAtomMapNum(0)
        part.GetAtomWithIdx(head).SetAtomMapNum(degree if i == 0 else i)
        part.GetAtomWithIdx(tail).SetAtomMapNum(degree + 1 if i == degree - 1 else i + 1)
        combined = part if combined is None else Chem.CombineMols(combined, part)
    # Zip caps through the same mechanism to preserve neighbor ordering at
    # stereocenters, even if a wildcard itself is one of their ligands.
    cap = "[H]" if end_group == "hydrogen" else "C"
    for label in (degree, degree + 1):
        combined = Chem.CombineMols(combined, _parse_smiles(f"{cap}[*:{label}]"))
    try:
        with rdBase.BlockLogs():
            chain = Chem.molzip(combined)
            Chem.SanitizeMol(chain)
            chain = Chem.RemoveHs(chain)
            Chem.AssignStereochemistry(chain, cleanIt=True, force=True)
    except (ValueError, RuntimeError) as exc:
        raise InputError("链节连接或封端后价态不合法，请检查 * 的位置、端点氢数和链节结构。") from exc
    # Preserve a connected graph and a precise stoichiometry, never export dummies.
    cap_name = "氢" if end_group == "hydrogen" else "甲基"
    source = f"链节 SMILES · 聚合度 {degree} · 两端{cap_name}封端"
    notes = [f"按尾 → 下一链节头的方向重复 {degree} 次，连接处均为单键，两端{cap_name}封端。",
             "这是有限长度的封端分子链，不是周期性无限聚合物；不会模拟聚合反应或自动改变链节内的键。",
             "重复复制输入中的立体信息，不自动生成无规 / 间规序列；长链构象采样需另行充分检查。"]
    return _check_structure(chain, source, notes)


def _validate_settings(s: Settings) -> None:
    if s.preprocess_method not in PREPROCESS_METHODS:
        raise InputError("预处理方法请选择 auto、mmff94s、uff、xtb、tblite、gfnff 或 none。")
    if s.boron_preprocess not in ("strict", "relaxed", "projected", "skip"):
        raise InputError("硼预处理请选择 strict、relaxed、projected 或 skip。")
    if (type(s.boron_bond_tolerance) not in (int, float) or not math.isfinite(s.boron_bond_tolerance)
            or not 0.01 <= s.boron_bond_tolerance <= 1.0):
        raise InputError("放宽保键的容差请输入 0.01–1.00 Å。")
    for label, value in (("内存", s.memory_gb), ("CPU 核数", s.processors), ("构象数", s.conformers)):
        limit = 1024 if label == "内存" else 256 if label == "CPU 核数" else 100
        if type(value) is not int or not 1 <= value <= limit:
            raise InputError(f"{label}必须是 1–{limit} 的整数。")
    for label, token in (("方法", s.method), ("基组", s.basis)):
        if not re.fullmatch(r"[A-Za-z0-9+*(),.\-]+", token):
            raise InputError(f"{label}只能包含英文、数字及 + * ( ) , . -，且不能留空。")
    if s.basis.casefold() in {"gen", "genecp", "chkbasis"}:
        raise InputError("Gen / GenECP / ChkBasis 需要额外基组或检查点；本工具请填写完整内置基组名。")
    if s.method.casefold() in {"oniom", "external", "read", "gen"}:
        raise InputError("该计算方法需要额外输入段，本工具暂不支持。")
    if s.job not in ("optfreq", "opt", "sp", "freq"):
        raise InputError("未知的计算任务。")
    if not re.fullmatch(r"[A-Za-z0-9_\-]{1,80}", s.checkpoint):
        raise InputError("检查点名称请使用 1–80 位英文字母、数字、下划线或连字符（不含扩展名）。")
    if not s.title.strip() or not s.title.isascii() or any(ord(c) < 32 or ord(c) == 127 for c in s.title):
        raise InputError("Gaussian 标题请输入单行可打印英文或数字。")
    if len(s.title) > 200:
        raise InputError("标题不能超过 200 字符。")
    if s.extra_keywords and (not s.extra_keywords.isascii() or
                            any(ord(c) < 32 or ord(c) == 127 for c in s.extra_keywords)):
        raise InputError("附加关键词必须是单行 ASCII 文本。")
    if len(s.extra_keywords) > 1000 or any(c in s.extra_keywords for c in "#%/@!\\"):
        raise InputError("附加关键词过长或包含非法的路由分隔符。")
    # These keywords require extra input sections or conflict with the generated geometry.
    prohibited = (r"\b(?:geom|guess|pseudo|pseudopotential|scrf|gen|genecp|chkbasis|"
                  r"extrabasis|basis|connectivity|counterpoise|oniom|external|"
                  r"modredundant|read|readfc|readopt|readfreeze|readanharm|"
                  r"qst2|qst3|addredundant|cards|include)\b|@")
    if re.search(prohibited, s.extra_keywords, re.I):
        raise InputError("附加关键词含需要检查点、额外数据段或特殊几何的选项，当前版本不支持。")
    if s.charge is not None and type(s.charge) is not int:
        raise InputError("电荷必须是整数。")
    if s.multiplicity is not None and (type(s.multiplicity) is not int or s.multiplicity < 1):
        raise InputError("自旋多重度必须是正整数。")


PREPROCESS_METHODS = {
    "auto": "自动（MMFF94s → UFF）",
    "mmff94s": "MMFF94s（有机分子）",
    "uff": "UFF（元素范围更广）",
    "xtb": "xTB GFN2-xTB（中等/大体系）",
    "tblite": "GFN2-xTB / tblite（严格优化）",
    "gfnff": "xTB GFN-FF（大体系/金属初猜）",
    "none": "仅 ETKDG（跳过力场）",
}

PREPROCESS_SCENARIOS = {
    "auto": "推荐默认：有 MMFF94s 参数时用于常见有机分子，否则尝试 UFF；完全离线。",
    "mmff94s": "适合常见有机分子、构象筛选和小分子；对金属、硼及重元素覆盖有限。",
    "uff": "适合需要更广元素覆盖的结构，包括部分无机/金属体系；参数是通用近似，不等于量化优化。",
    "xtb": "需要安装 xTB。GFN2-xTB 适合中等大小有机、主族和部分金属结构，通常比经典力场更能描述电荷和非键作用。",
    "tblite": "通过 xTB 内置 tblite 运行 GFN2-xTB，采用 tight 几何收敛标准及所设 CPU 核数；适合有机、主族和部分金属体系。是同一量化方法的另一后端，不提高理论级别；需要支持 tblite 的 xTB。",
    "gfnff": "需要安装 xTB。GFN-FF 速度快、适合大分子和聚合物初始几何；不适合作为最终能量或电子结构结论。",
    "none": "只生成 ETKDG 三维坐标，不做力场优化；适合快速检查拓扑或交给 Gaussian/其他量化程序自行优化。",
}


def _user_xtb_executable() -> str:
    # GUI processes already running when xTB was installed may have stale env.
    if os.name == "nt":
        import winreg
        try:
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, "Environment") as key:
                return str(winreg.QueryValueEx(key, "XTB_EXE")[0]).strip()
        except OSError:
            pass
    return ""


def _xtb_executable() -> str:
    configured = os.environ.get("XTB_EXE", "").strip()
    candidates = [configured] if configured else []
    if not configured:
        candidates.append(_user_xtb_executable())
    candidates.extend(("xtb.exe", "xtb"))
    for candidate in candidates:
        if not candidate:
            continue
        path = Path(candidate)
        if path.is_file():
            return str(path.resolve())
        found = shutil.which(candidate)
        if found:
            return found
    raise InputError("已选择 xTB 预处理，但找不到 xtb.exe。请安装 xTB，并将其加入 PATH，或设置环境变量 XTB_EXE。")


def _read_xyz(path: Path):
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        count = int(lines[0].strip())
        rows = [line.split() for line in lines[2:2 + count]]
        if len(rows) != count:
            raise ValueError("XYZ 原子数不完整")
        return [(row[0], tuple(float(v) for v in row[1:4])) for row in rows]
    except (OSError, ValueError, IndexError) as exc:
        raise InputError("xTB 未返回可读取的优化坐标（xtbopt.xyz）。") from exc


def _xtb_optimize(mol, conf_id: int, method: str, charge: int, multiplicity: int,
                  processors: int = 1):
    executable = _xtb_executable()
    with tempfile.TemporaryDirectory(prefix="gaussian_xtb_") as folder_name:
        folder = Path(folder_name)
        xyz = folder / "input.xyz"
        conf = mol.GetConformer(conf_id)
        rows = [str(mol.GetNumAtoms()), "GaussianInputTool xTB preprocessing"]
        for atom in mol.GetAtoms():
            point = conf.GetAtomPosition(atom.GetIdx())
            rows.append(f"{atom.GetSymbol()} {point.x:.10f} {point.y:.10f} {point.z:.10f}")
        xyz.write_text("\n".join(rows) + "\n", encoding="ascii")
        flags = ["--gfnff"] if method == "gfnff" else ["--gfn", "2"]
        flags += ["--opt"]
        if method == "tblite":
            flags += ["tight", "--tblite"]
        flags += ["--chrg", str(charge), "--uhf", str(max(0, multiplicity - 1)),
                  "--parallel", str(processors)]
        environment = {**os.environ, "OMP_NUM_THREADS": str(processors),
                       "MKL_NUM_THREADS": str(processors), "OPENBLAS_NUM_THREADS": str(processors)}
        try:
            completed = subprocess.run([executable, str(xyz), *flags], cwd=folder,
                                       capture_output=True, text=True, encoding="utf-8", errors="replace",
                                       timeout=600, env=environment,
                                       creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        except subprocess.TimeoutExpired as exc:
            raise InputError("xTB 预处理超过 600 秒；可减少体系规模或改用 UFF/MMFF94s。") from exc
        except OSError as exc:
            raise InputError(f"无法启动 xTB：{exc}") from exc
        output = folder / "xtbopt.xyz"
        log = completed.stdout + "\n" + completed.stderr
        if method == "tblite" and re.search(
                r"(?:not compiled with|compiled without|not available|not supported|unknown option|unrecognized).*tblite|"
                r"tblite.*(?:not available|not supported|not enabled)", log, re.I):
            raise InputError("当前 xtb.exe 不支持 tblite 后端，请安装启用 tblite 的 xTB；未自动切换其他方法。")
        if completed.returncode != 0 or not output.is_file():
            detail = (completed.stderr or completed.stdout).strip()[-2000:]
            raise InputError("xTB 预处理失败：" + (detail or "未生成 xtbopt.xyz"))
        if (not (folder / ".xtboptok").is_file() or
                (folder / "NOT_CONVERGED").exists()):
            raise InputError("xTB 几何优化未收敛；已停止导出，避免将未收敛坐标当作优化结果。")
        optimized = _read_xyz(output)
        expected = [atom.GetSymbol() for atom in mol.GetAtoms()]
        if len(optimized) != len(expected) or [row[0] for row in optimized] != expected:
            raise InputError("xTB 返回的原子顺序或数量与输入不一致，已停止导出。")
        for _, point in optimized:
            if len(point) != 3 or not all(math.isfinite(value) for value in point):
                raise InputError("xTB 返回了非有限三维坐标。")
        for index, (_, point) in enumerate(optimized):
            conf.SetAtomPosition(index, point)
        return completed.stdout


def _external_preprocessing(mol, conf_id, method, charge, multiplicity, processors):
    _xtb_optimize(mol, conf_id, method, charge, multiplicity, processors)
    label = {"xtb": "xTB GFN2-xTB", "gfnff": "xTB GFN-FF",
             "tblite": "GFN2-xTB / tblite（tight）"}[method]
    return label, [f"已用 {label} 优化一个初始构象；xTB 坐标仅作 Gaussian 起始结构，不能替代最终量化优化。"]


def _geometry(mol, count: int, preprocess_method="auto", charge=0, multiplicity=1, processors=1, *, on_stage=None):
    stage = on_stage if on_stage is not None else lambda _: None
    if mol.GetNumAtoms() == 1:
        conf = Chem.Conformer(1)
        conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
        conf.Set3D(True)
        mol.AddConformer(conf, assignId=True)
        return 0, "单原子（原点）", None, []
    if mol.GetNumAtoms() == 2:
        # ETKDG handles diatomic bond lengths while force-field coverage often does not.
        count = 1
    if preprocess_method in ("xtb", "tblite", "gfnff"):
        # Only one conformer is optimized by these backends. Avoid generating
        # unused conformers, especially for large rings and polymers.
        count = 1
    stage("生成三维构象")
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261117
    params.numThreads = 1
    params.pruneRmsThresh = 0.25
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=count, params=params))
    if not ids:
        params.useRandomCoords = True
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=count, params=params))
    if not ids:
        raise InputError("三维构象生成失败；请检查结构或使用专用建模软件。")
    notes = []
    if preprocess_method in ("xtb", "tblite", "gfnff"):
        stage(f"运行 {PREPROCESS_METHODS[preprocess_method]} 预处理")
        label, notes = _external_preprocessing(mol, ids[0], preprocess_method, charge, multiplicity, processors)
        return ids[0], f"ETKDGv3 + {label}", None, notes
    if preprocess_method == "none":
        notes.append("已跳过力场预优化，仅使用 ETKDG 三维坐标；请在 Gaussian 中检查并优化键长、键角和非键接触。")
        return ids[0], "ETKDGv3（未进行力场预优化）", None, notes
    use_mmff = preprocess_method in ("auto", "mmff94s") and AllChem.MMFFHasAllMoleculeParams(mol)
    use_uff = preprocess_method in ("auto", "uff") and AllChem.UFFHasAllMoleculeParams(mol)
    if preprocess_method == "mmff94s" and not use_mmff:
        raise InputError("已选择 MMFF94s，但该结构缺少 MMFF94s 参数；可改用 UFF、xTB 或自动模式。")
    if preprocess_method == "uff" and not use_uff:
        raise InputError("已选择 UFF，但该结构缺少 UFF 参数；可改用 xTB 或仅 ETKDG。")
    candidates = []
    if use_mmff or use_uff:
        stage("预优化并筛选构象")
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s") if use_mmff else None
        for cid in ids:
            try:
                ff = (AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid) if use_mmff
                      else AllChem.UFFGetMoleculeForceField(mol, confId=cid))
                if ff is None:
                    continue
                status = ff.Minimize(maxIts=1500)
                energy = float(ff.CalcEnergy())
                if math.isfinite(energy):
                    candidates.append((status != 0, energy, cid))
            except (RuntimeError, ValueError):
                continue
    if candidates:
        unconv, energy, cid = min(candidates)
        method = "MMFF94s" if use_mmff else "UFF"
        if unconv:
            notes.append("力场预优化未收敛，坐标仅用于起始结构；请在 Gaussian 中优化并检查。")
        return cid, f"ETKDGv3 + {method}（{len(ids)} 个构象中选取）", energy, notes
    notes.append("该结构缺少可用力场参数，使用未预优化的 ETKDG 三维坐标；请检查键长并进行量化优化。")
    return ids[0], "ETKDGv3（未进行力场预优化）", None, notes


def _separate_disconnected_fragments(mol, conf):
    """Lay out disconnected SMILES components side by side without changing them internally."""
    fragments = Chem.GetMolFrags(mol, asMols=False)
    if len(fragments) < 2:
        return
    cursor = 0.0
    gap = 3.0
    for fragment in fragments:
        points = np.asarray([conf.GetAtomPosition(i) for i in fragment], dtype=float)
        mins, maxs = points.min(axis=0), points.max(axis=0)
        shift = np.array([cursor - mins[0], -0.5 * (mins[1] + maxs[1]), -0.5 * (mins[2] + maxs[2])])
        for idx in fragment:
            point = np.asarray(conf.GetAtomPosition(idx), dtype=float) + shift
            conf.SetAtomPosition(idx, point)
        cursor += float(maxs[0] - mins[0]) + gap


def _ferrocene_geometry(mol):
    """Build a staggered sandwich starting geometry and explicit Gaussian connectivity."""
    fe = next((atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "Fe"), None)
    rings = [tuple(int(i) for i in ring)
             for ring in Chem.GetSymmSSSR(mol) if len(ring) == 5 and
             all(mol.GetAtomWithIdx(int(i)).GetAtomicNum() == 6 for i in ring)]
    if fe is None or len(rings) != 2:
        raise InputError("二茂铁模板内部结构无效。")
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    conf.SetAtomPosition(fe, (0.0, 0.0, 0.0))
    height, radius, bond_ch = 1.66, 1.20, 1.08
    for ring_index, ring in enumerate(rings):
        z = height if ring_index == 0 else -height
        offset = 0.0 if ring_index == 0 else math.pi / 5.0
        for k, carbon in enumerate(ring):
            angle = offset + 2.0 * math.pi * k / 5.0
            x, y = radius * math.cos(angle), radius * math.sin(angle)
            conf.SetAtomPosition(int(carbon), (x, y, z))
            hydrogens = [n.GetIdx() for n in mol.GetAtomWithIdx(int(carbon)).GetNeighbors()
                         if n.GetAtomicNum() == 1]
            for hydrogen in hydrogens:
                conf.SetAtomPosition(hydrogen, ((radius + bond_ch) * math.cos(angle),
                                                (radius + bond_ch) * math.sin(angle), z))
    cid = mol.AddConformer(conf, assignId=True)
    connectivity = {atom.GetIdx(): {} for atom in mol.GetAtoms()}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = float(bond.GetBondTypeAsDouble())
        connectivity[i][j] = connectivity[j][i] = order
    for ring in rings:
        for carbon in ring:
            connectivity[fe][int(carbon)] = connectivity[int(carbon)][fe] = 1.0
    return cid, connectivity


def _trigonal_boron_dihedrals(mol, conf, preprocess="strict", bond_tolerance=0.15, notes=None):
    """Choose torsions that place each trigonal B substituent in one plane."""
    from rdkit.Chem import rdMolTransforms

    constraints = []
    has_trigonal_boron = False
    for boron_atom in mol.GetAtoms():
        if boron_atom.GetAtomicNum() != 5 or boron_atom.GetDegree() != 3:
            continue
        has_trigonal_boron = True
        boron = boron_atom.GetIdx()
        neighbors = [atom.GetIdx() for atom in boron_atom.GetNeighbors()]
        anchors = []
        for idx in neighbors:
            atom = mol.GetAtomWithIdx(idx)
            outer = [n.GetIdx() for n in atom.GetNeighbors()
                     if n.GetIdx() != boron and n.GetAtomicNum() > 1]
            if outer:
                anchors.append((atom.GetAtomicNum() == 6, idx, outer))
        if not anchors:
            continue
        # Prefer an attached carbon (e.g. the aryl carbon in an aryl boronate)
        # as the torsion reference. Freeze both remaining B substituents.
        anchor_candidates = sorted(anchors, key=lambda item: (not item[0], item[1]))
        _, anchor, outer_candidates = anchor_candidates[0]
        others = sorted(idx for idx in neighbors if idx != anchor)
        # Pick the reference side of the anchor bond for which the required
        # 180-degree torsion is closest to the starting conformation.
        outer = min(outer_candidates, key=lambda idx: abs(
            (rdMolTransforms.GetDihedralDeg(conf, idx, anchor, boron, others[0]) - 180 + 180) % 360 - 180))
        primary = (outer, anchor, boron, others[0])
        secondary = (outer, anchor, boron, others[1])
        constraints.append((*primary, 180.0))
        current = rdMolTransforms.GetDihedralDeg(conf, *secondary)
        target = min((0.0, 180.0), key=lambda angle: abs((current - angle + 180) % 360 - 180))
        constraints.append((*secondary, target))

    if preprocess == "projected" and has_trigonal_boron:
        try:
            from boron_projection import project_boron_geometry
            coords, drift, error, height = project_boron_geometry(
                mol, conf.GetPositions(), constraints, min(0.03, bond_tolerance))
        except ImportError as exc:
            raise InputError("几何投影需要 SciPy；请重新运行 start.bat 安装依赖。") from exc
        except ValueError as exc:
            raise InputError(f"硼几何投影未达标：{exc}") from exc
        for i, point in enumerate(coords):
            conf.SetAtomPosition(i, point)
        if notes is not None:
            notes.append(f"硼几何投影：最大键长变化 {drift:.6f} Å，最大二面角偏差 {error:.4f}°，"
                         f"B 离邻原子平面的最大距离 {height:.6f} Å。键长以整体预处理后的坐标为参考。")
        return constraints
    if not constraints:
        return constraints
    if preprocess == "skip":
        if notes is not None:
            notes.append("已跳过硼平面化预处理，保留构象生成后的坐标；初始二面角可能不共面。优化输入仍设置并冻结主二面角 180° 及辅助共面二面角；未添加键角约束。")
        return constraints

    # Constrained force-field relaxation keeps ring closures and all existing
    # bonds connected while moving the complete structure toward the targets.
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf.GetId())
    if ff is None:
        raise InputError("含三配位硼的结构缺少 UFF 参数，无法在保持键连的情况下生成 180° 二面角起始坐标。")
    ff.Minimize(maxIts=1500, energyTol=1.0e-8, forceTol=1.0e-5)
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf.GetId())
    if ff is None:
        raise InputError("无法建立硼平面化约束所需的 UFF 力场。")
    bonds = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        length = float(np.linalg.norm(np.asarray(conf.GetAtomPosition(i)) -
                                      np.asarray(conf.GetAtomPosition(j))))
        bonds.append((i, j, length))
        window = bond_tolerance if preprocess == "relaxed" else 0.0
        ff.UFFAddDistanceConstraint(i, j, False, max(0.1, length - window), length + window,
                                    1.0e4 if preprocess == "relaxed" else 1.0e5)
    for i, j, k, l, target in constraints:
        ff.UFFAddTorsionConstraint(i, j, k, l, False, target, target, 1.0e5)
    ff.Initialize()
    ff.Minimize(maxIts=6000 if preprocess == "relaxed" else 3000, energyTol=1.0e-9, forceTol=1.0e-6)
    maximum_drift = 0.0
    for i, j, original in bonds:
        current = float(np.linalg.norm(np.asarray(conf.GetAtomPosition(i)) -
                                       np.asarray(conf.GetAtomPosition(j))))
        maximum_drift = max(maximum_drift, abs(current - original))
        limit = bond_tolerance if preprocess == "relaxed" else 0.03
        if abs(current - original) > limit:
            raise InputError(f"硼平面化预处理导致键 {i + 1}-{j + 1} 长度变化 {current - original:+.3f} Å，超过 {limit:.3f} Å；"
                             "可选择“放宽保键”并调整容差，或“跳过平面化预处理”。Gaussian 优化中的共面二面角冻结仍保留。")
    for torsion in constraints:
        i, j, k, l, target = torsion
        current = rdMolTransforms.GetDihedralDeg(conf, i, j, k, l)
        error = abs((current - target + 180) % 360 - 180)
        if error > 0.5:
            raise InputError(f"硼相关二面角 {i + 1}-{j + 1}-{k + 1}-{l + 1} 未能达到 {target:.0f}°（当前 {current:.1f}°）；"
                             "可选择“跳过平面化预处理”，将共面约束交给 Gaussian 优化。")
    if notes is not None:
        notes.append(f"硼平面化已完成（{'放宽保键' if preprocess == 'relaxed' else '严格保键'}）；相对 UFF 初步松弛结构的最大键长变化 {maximum_drift:.4f} Å。"
                     "优化输入冻结主二面角 180° 及辅助共面二面角；未添加键角约束。")
    return constraints


def generate(value: str, mode: str = "auto", settings: Settings | None = None, *, on_stage=None) -> Result:
    """Generate Gaussian text; optionally report real processing stages to the caller."""
    stage = on_stage if on_stage is not None else lambda _: None
    stage("解析结构与计算参数")
    s = settings or Settings()
    _validate_settings(s)
    if s.engine not in ('gaussian', 'orca'):
        raise InputError('未知计算程序。')
    if mode == 'cif' or (mode == 'auto' and value.lower().endswith('.cif')):
        from crystal_input import generate_cif
        result = generate_cif(value, s)
    else:
        result = _generate_molecule(value, mode, s, stage)
    if s.engine == 'orca':
        from orca_support import render
        return render(result, s)
    return result


def _generate_molecule(value, mode, s, stage):
    if mode == "complex":
        from complex_core import generate_complex
        return generate_complex(value, s, on_stage=stage)
    if mode == "polymer":
        mol, source, warnings = build_polymer(value, s.degree, s.end_group)
    elif mode == "ferrocene":
        if value.strip().casefold() not in ("ferrocene", "二茂铁", "[fe+2].[ch-]1cccc1.[ch-]1cccc1"):
            raise InputError("二茂铁模板模式请使用输入值 ferrocene，或选择界面中的“二茂铁模板”。")
        mol = _parse_smiles("[Fe+2].[cH-]1cccc1.[cH-]1cccc1")
        mol, source, warnings = _check_structure(mol, "二茂铁夹心模板", [])
        warnings = [warning for warning in warnings if "多个独立片段" not in warning]
    else:
        mol, source, warnings = resolve_molecule(value, mode)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    formal_charge = Chem.GetFormalCharge(mol)
    mol = Chem.AddHs(mol)
    charge = formal_charge if s.charge is None else s.charge
    electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
    radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    has_metal = any(21 <= a.GetAtomicNum() <= 30 or 39 <= a.GetAtomicNum() <= 48 or
                    57 <= a.GetAtomicNum() <= 80 or a.GetAtomicNum() >= 89 for a in mol.GetAtoms())
    if has_metal and s.multiplicity is None and mode != "ferrocene":
        raise InputError("含过渡金属或 f 区元素，无法可靠推断自旋态；请手动填写自旋多重度。")
    multiplicity = s.multiplicity
    if multiplicity is None:
        multiplicity = (1 if mode == "ferrocene" else
                        (radicals + 1 if radicals else (1 if electrons % 2 == 0 else 2)))
        if radicals > 1:
            warnings.append(f"检测到 {radicals} 个自由基电子，自动采用高自旋多重度 {multiplicity}；请按目标电子态确认。")
    if electrons < 1 or multiplicity > electrons + 1 or (electrons - multiplicity + 1) % 2:
        raise InputError(f"电荷 / 多重度不相容：当前有 {electrons} 个电子，多重度为 {multiplicity}。"
                         "偶数电子配奇数多重度，奇数电子配偶数多重度。")
    if charge != formal_charge:
        warnings.append(f"指定电荷 {charge:+d} 与 SMILES 形式电荷 {formal_charge:+d} 不同；只改变电子数，不改变氢数。")
    warnings.append("自旋多重度需按目标电子态核对；电子数奇偶校验不能判定基态自旋。")
    if has_metal:
        warnings.append("支持以 RDKit 可解析的配位 SMILES 输入金属配合物（推荐使用 -> 表示配位键）；自动生成的配位几何仅为起始猜测，须检查配位数、键长、构型、方法/基组/赝势和自旋态。")
    if s.job in ("sp", "freq"):
        warnings.append("当前任务不进行量化结构优化；单点 / 频率计算前请确认几何满足研究目的。")
    connectivity = None
    if mode == "ferrocene":
        stage("构建二茂铁夹心结构")
        cid, connectivity = _ferrocene_geometry(mol)
        geometry_method, energy = "二茂铁模板（交错式 Cp–Fe–Cp 夹心初始几何）", None
        if s.preprocess_method in ("xtb", "tblite", "gfnff"):
            stage(f"运行 {PREPROCESS_METHODS[s.preprocess_method]} 预处理")
            label, notes = _external_preprocessing(mol, cid, s.preprocess_method, charge, multiplicity, s.processors)
            geometry_method += " + " + label
            warnings.extend(notes)
        warnings.append("已将 Fe 与两个环戊二烯基环的 10 个碳写入 Gaussian Connectivity 连通表，作为优化的分子拓扑；Gaussian 会对整个体系统一进行量子化学优化。初始模板及自旋态仍需按研究目标检查。")
    else:
        cid, geometry_method, energy, notes = _geometry(
            mol, s.conformers, s.preprocess_method, charge, multiplicity, s.processors, on_stage=stage)
        warnings.extend(notes)
        conf = mol.GetConformer(cid)
        _separate_disconnected_fragments(mol, conf)
    conf = mol.GetConformer(cid)
    stage("检查起始几何并组装输入文件")
    boron_dihedrals = _trigonal_boron_dihedrals(mol, conf, s.boron_preprocess, s.boron_bond_tolerance, warnings)
    if boron_dihedrals:
        warnings.append(f"检测到 {len(boron_dihedrals) // 2} 个三配位硼中心；所有原有键写入 Connectivity，保留环闭合拓扑。")
        if connectivity is None:
            connectivity = {atom.GetIdx(): {} for atom in mol.GetAtoms()}
            for bond in mol.GetBonds():
                i, j = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
                connectivity[i][j] = 0.1 if bond.GetBondType() == Chem.BondType.DATIVE else bond.GetBondTypeAsDouble()
    if any(a.GetAtomicNum() == 5 and a.GetDegree() == 4 for a in mol.GetAtoms()):
        warnings.append("检测到四配位硼；未添加平面约束，因为该配位通常为四面体结构。")
    warnings.append("构象筛选仅提供初始几何，不保证全局最低能构象；力场能量不是 Gaussian 能量。")
    if boron_dihedrals and s.job == "optfreq":
        route_job = "Opt=(ModRedundant) Freq"
    elif boron_dihedrals and s.job == "opt":
        route_job = "Opt=(ModRedundant)"
    else:
        route_job = {"optfreq": "Opt Freq", "opt": "Opt", "sp": "SP", "freq": "Freq"}[s.job]
    route = f"#p {s.method}/{s.basis} {route_job}"
    if connectivity is not None:
        route += " Geom=Connectivity"
    if s.extra_keywords.strip():
        route += " " + s.extra_keywords.strip()
    lines = [f"%chk={s.checkpoint}.chk", f"%mem={s.memory_gb}GB", f"%nprocshared={s.processors}",
             route, "", s.title.strip(), "", f"{charge} {multiplicity}"]
    conf = mol.GetConformer(cid)
    for atom in mol.GetAtoms():
        point = conf.GetAtomPosition(atom.GetIdx())
        if not all(math.isfinite(x) for x in (point.x, point.y, point.z)):
            raise InputError("三维坐标包含非有限值，请检查分子结构。")
        symbol = atom.GetSymbol()
        if atom.GetIsotope():
            symbol += f"(Iso={atom.GetIsotope()})"
        lines.append(f"{symbol:<12} {point.x: .8f} {point.y: .8f} {point.z: .8f}")
    if connectivity is not None:
        lines.append("")
        for atom_idx in range(mol.GetNumAtoms()):
            entries = [str(atom_idx + 1)]
            for neighbor_idx, order in sorted(connectivity[atom_idx].items()):
                if neighbor_idx > atom_idx:
                    entries.extend((str(neighbor_idx + 1), f"{order:.1f}"))
            lines.append(" ".join(entries))
    if boron_dihedrals and s.job in ("opt", "optfreq"):
        lines.append("")
        # Gaussian 16 treats a value and F on the same line as a coordinate
        # build operation. Set the value first, then freeze on a separate line.
        for i, j, k, l, target in boron_dihedrals:
            coordinate = f"D {i + 1} {j + 1} {k + 1} {l + 1}"
            lines.extend((f"{coordinate} {target:.1f}", f"{coordinate} F"))
    return Result("\n".join(lines) + "\n\n", canonical, rdMolDescriptors.CalcMolFormula(mol),
                  charge, multiplicity, mol.GetNumAtoms(), geometry_method, energy, source, warnings,
                  value.strip() if mode == "polymer" else None,
                  s.degree if mode == "polymer" else None,
                  s.end_group if mode == "polymer" else None)


def default_checkpoint(value: str, mode: str = "auto", degree: int = 5,
                       end_group: str = "hydrogen") -> str:
    if mode == "polymer":
        key = json.dumps([value.strip(), degree, end_group], ensure_ascii=True)
        return f"poly_n{degree}_" + hashlib.sha256(key.encode("utf-8")).hexdigest()[:10]
    return "mol_" + hashlib.sha256(value.encode("utf-8")).hexdigest()[:10]
