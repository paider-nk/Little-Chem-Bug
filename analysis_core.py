"""Read Gaussian results and export auditable numerical analyses and figures."""
from __future__ import annotations

import csv
import html
import json
import math
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from rdkit import Chem

from gaussian_core import InputError

HARTREE_EV = 27.211386245988
HARTREE_KJMOL = 2625.4996394799
BOHR_ANGSTROM = 0.529177210903
FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[DEde][-+]?\d+)?"


def number(value):
    return float(value.replace("D", "E").replace("d", "e"))


@dataclass
class Analysis:
    path: str
    normal: bool = False
    route: str = ""
    method_basis: str | None = None
    charge: int | None = None
    multiplicity: int | None = None
    energy_hartree: float | None = None
    energy_kind: str = "SCF"
    zpe_corrected_hartree: float | None = None
    enthalpy_hartree: float | None = None
    free_energy_hartree: float | None = None
    temperature_k: float | None = None
    pressure_atm: float | None = None
    solvent: str | None = None
    dipole_debye: float | None = None
    optimization_complete: bool = False
    frequencies_cm1: list[float] = field(default_factory=list)
    ir_intensities: list[float] = field(default_factory=list)
    atoms: list[dict] = field(default_factory=list)
    mulliken: list[dict] = field(default_factory=list)
    orbitals: dict = field(default_factory=dict)
    counterpoise_energy_hartree: float | None = None
    bsse_hartree: float | None = None
    warnings: list[str] = field(default_factory=list)


def orbital_summary(energies, occupied, *, occupations=None):
    if not energies or not 0 <= occupied <= len(energies) or not all(map(math.isfinite, energies)):
        return {}
    is_occupied = ([value > 1e-8 for value in occupations] if occupations is not None
                   else [i < occupied for i in range(len(energies))])
    occupied_indices = [i for i, filled in enumerate(is_occupied) if filled]
    empty_indices = [i for i, filled in enumerate(is_occupied) if not filled]
    homo = max(occupied_indices, key=energies.__getitem__, default=None)
    lumo = min(empty_indices, key=energies.__getitem__, default=None)
    homo_ev = energies[homo] * HARTREE_EV if homo is not None else None
    lumo_ev = energies[lumo] * HARTREE_EV if lumo is not None else None
    return {"occupied": len(occupied_indices), "is_occupied": is_occupied, "energies_hartree": energies,
            "homo_index": homo + 1 if homo is not None else None,
            "lumo_index": lumo + 1 if lumo is not None else None,
            "homo_ev": homo_ev, "lumo_ev": lumo_ev,
            "gap_ev": lumo_ev - homo_ev if homo_ev is not None and lumo_ev is not None else None}


def _gaussian_orbitals(text):
    # Alpha/Beta and occupied/virtual rows must come from the same printed block.
    pattern = re.compile(r"^[ \t]*(Alpha|Beta)[ \t]+(occ|virt)\.[ \t]+eigenvalues[ \t]*--(.*)$", re.M)
    boundaries = list(re.finditer(r"^[ \t]*Normal termination[^\n]*$", text, re.M))
    if boundaries:
        tail = text[boundaries[-1].end():]
        # A later Link1 or concatenated job may have failed before printing any MOs.
        if pattern.search(tail) or re.search(r"Error termination|Entering Link 1|SCF Done:", tail):
            text = tail
        elif len(boundaries) > 1:
            text = text[boundaries[-2].end():boundaries[-1].start()]
    last, block, previous = None, {}, None
    order = {('Alpha', 'occ'): 0, ('Alpha', 'virt'): 1, ('Beta', 'occ'): 2, ('Beta', 'virt'): 3}
    for match in pattern.finditer(text):
        key = (match[1], match[2])
        if previous is not None and (order[key] < order[last] or text[previous:match.start()].strip()):
            block = {}
        block.setdefault(key, []).extend(number(value) for value in re.findall(FLOAT, match[3]))
        last, previous = key, match.end()
    results = {}
    for spin in ('Alpha', 'Beta'):
        occupied = block.get((spin, 'occ'), [])
        virtual = block.get((spin, 'virt'), [])
        if summary := orbital_summary(occupied + virtual, len(occupied)):
            results[spin.lower()] = summary
    return results


def _orca_orbitals(text):
    headers = list(re.finditer(r'^[ \t]*ORBITAL ENERGIES[ \t]*$', text, re.M))
    if not headers:
        return {}
    row_pattern = re.compile(r'[ \t]*(\d+)[ \t]+(' + FLOAT + r')[ \t]+(' + FLOAT + r')[ \t]+' + FLOAT + r'[ \t]*')
    result, rows, spin = {}, [], 'alpha'
    for line in text[headers[-1].end():].splitlines():
        if not line.strip() or not line.strip().strip('-'):
            continue
        if line.strip() == 'SPIN DOWN ORBITALS':
            result[spin] = rows
            rows, spin = [], 'beta'
        elif match := row_pattern.fullmatch(line):
            index = int(match[1])
            if index != len(rows):
                break
            rows.append((number(match[2]), number(match[3])))
        elif rows:
            break
        elif not (line.strip(' -') == '' or line.split()[0] in ('NO', 'SPIN')):
            break
    result[spin] = rows
    summaries = {}
    for spin, rows in result.items():
        if rows:
            occupations, energies = zip(*rows)
            if summary := orbital_summary(list(energies), 0, occupations=occupations):
                summaries[spin] = summary
    return summaries


def parse_orca(path, text):
    a = Analysis(str(path))
    a.normal = bool(re.search(r'\*+ORCA TERMINATED NORMALLY\*+\s*(?:TOTAL RUN TIME:[^\n]*)?\s*\Z', text))
    if not a.normal:
        a.warnings.append('ORCA 未正常结束。')
    routes = re.findall(r'^\s*\|\s*\d+>\s*!\s*(.+)$', text, re.M)
    a.route = ' '.join(routes)
    tokens = a.route.lower().split()
    solvent = [t for t in tokens if t.startswith(('cpcm(', 'smd('))]
    a.solvent = ' '.join(solvent) or None
    ignored = {'opt','sp','freq','numfreq','tightopt'} | set(solvent)
    signature = [t for t in tokens if t not in ignored]
    # Only compare simple, standalone inputs whose theory is fully on the ! line.
    advanced = re.search(r'^\s*\|\s*\d+>\s*(?:%(?!pal\b|maxcore\b|geom\b)|\$)', text, re.M|re.I)
    if len(signature) >= 2 and not advanced:
        a.method_basis = 'orca:' + ' '.join(signature)
    fields = {'energy_hartree':r'FINAL SINGLE POINT ENERGY\s+', 'enthalpy_hartree':r'Total Enthalpy\s+\.+\s*',
              'free_energy_hartree':r'Final Gibbs free energy\s+\.+\s*', 'temperature_k':r'Temperature\s+\.+\s*',
              'pressure_atm':r'Pressure\s+\.+\s*','dipole_debye':r'Magnitude \(Debye\)\s*:\s*'}
    for field, prefix in fields.items():
        values = re.findall(prefix + '('+FLOAT+')', text)
        if values:
            setattr(a,field,number(values[-1]))
    a.energy_kind = 'ORCA'
    a.optimization_complete = 'OPTIMIZATION HAS CONVERGED' in text
    for field, key in [('charge','Total Charge\\s+Charge'), ('multiplicity','Multiplicity\\s+Mult')]:
        matches = re.findall(key + r'\s+\.+\s*(-?\d+)',text)
        if matches:
            setattr(a,field,int(matches[-1]))
    def last_table(title):
        matches = re.findall(title+r'\s*\n-+\n(.*?)(?:\n[ \t]*\n|\Z)', text, re.S)
        return matches[-1].splitlines() if matches else []
    for row in last_table(r'CARTESIAN COORDINATES \(ANGSTROEM\)'):
        tokens = row.split()
        if len(tokens)==4:
            symbol = tokens[0].capitalize()
            a.atoms.append(dict(index=len(a.atoms)+1, atomic_number=Chem.GetPeriodicTable().GetAtomicNumber(symbol), symbol=symbol, xyz=list(map(number,tokens[1:]))))
    for row in last_table(r'MULLIKEN ATOMIC CHARGES(?: AND SPIN POPULATIONS)?'):
        cells = row.replace(':',' ').split()
        if cells and cells[0].isdigit() and len(cells) in (3,4):
            a.mulliken.append(dict(index=int(cells[0])+1,symbol=cells[1],charge=number(cells[2]),spin=number(cells[3]) if len(cells)==4 else None))
    a.orbitals = _orca_orbitals(text)
    modes = re.findall(r'^\s*(\d+):\s+('+FLOAT+r')\s+cm\*\*-1',text,re.M)
    first = re.findall(r'first frequency considered to be a vibration is (\d+)',text)
    vib = [(int(i),number(v)) for i,v in modes if not first or int(i)>=int(first[-1]) or number(v)<0]
    a.frequencies_cm1 = [v for i,v in vib]
    intensities = dict(re.findall(r'^\s*(\d+):\s+'+FLOAT+r'\s+'+FLOAT+r'\s+('+FLOAT+r')\s+'+FLOAT+r'\s+\(',text,re.M))
    if vib and all(str(i) in intensities for i,v in vib):
        a.ir_intensities = [number(intensities[str(i)]) for i,v in vib]
    if any(v<0 for i,v in vib):
        a.warnings.append('振动分析含虚频，请核对最终结构。')
    return a


def parse_log(path: str | Path) -> Analysis:
    path = Path(path).resolve()
    text = path.read_text(encoding="utf-8", errors="replace")
    if 'ORCA TERMINATED NORMALLY' in text or ('O   R   C   A' in text and 'Program Version' in text) or 'FINAL SINGLE POINT ENERGY' in text:
        return parse_orca(path,text)
    a = Analysis(str(path))
    a.normal = bool(re.search(r"Normal termination[^\n]*\s*\Z", text))
    if not a.normal:
        a.warnings.append("日志未以 Normal termination 结束；当前结果可能不完整，不能作为完成的计算。")
    routes = re.findall(r"^\s*(#[^\n]*(?:\n(?!\s*-{5,})[^\n]*)*)\n\s*-{5,}", text, re.M)
    if routes:
        a.route = " ".join(routes[-1].split())
    match = re.search(r"\b([A-Za-z0-9+*().-]+/[A-Za-z0-9+*(),.\-]+)", a.route)
    if match:
        a.method_basis = match[1].casefold()
    solvents = re.findall(r"\bsolvent\s*=\s*([A-Za-z0-9-]+)", a.route, re.I)
    a.solvent = solvents[-1].casefold() if solvents else None
    charges = re.findall(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", text)
    if charges:
        a.charge, a.multiplicity = map(int, charges[-1])
    energies = list(re.finditer(r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(" + FLOAT + ")", text))
    if energies:
        a.energy_hartree = number(energies[-1][2])
        a.energy_kind = energies[-1][1]
        tail = text[energies[-1].end():]
        for kind, expression in (("MP2", r"EUMP2\s*=\s*(" + FLOAT + ")"),
                                 ("CCSD(T)", r"CCSD\(T\)\s*=\s*(" + FLOAT + ")")):
            correlated = re.findall(expression, tail)
            if correlated:
                a.energy_hartree, a.energy_kind = number(correlated[-1]), kind
    fields = {"zpe_corrected_hartree": "Sum of electronic and zero-point Energies",
              "enthalpy_hartree": "Sum of electronic and thermal Enthalpies",
              "free_energy_hartree": "Sum of electronic and thermal Free Energies",
              "counterpoise_energy_hartree": "Counterpoise corrected energy",
              "bsse_hartree": "BSSE energy"}
    for attribute, label in fields.items():
        matches = re.findall(re.escape(label) + r"\s*=\s*(" + FLOAT + ")", text, re.I)
        if matches:
            setattr(a, attribute, number(matches[-1]))
    conditions = re.findall(r"Temperature\s+(" + FLOAT + r")\s+Kelvin\.\s+Pressure\s+(" + FLOAT + ")", text)
    if conditions:
        a.temperature_k, a.pressure_atm = map(number, conditions[-1])
    dipoles = re.findall(r"Dipole moment[^\n]*\n[^\n]*Tot=\s*(" + FLOAT + ")", text)
    if dipoles:
        a.dipole_debye = number(dipoles[-1])
    a.optimization_complete = "Optimization completed" in text
    # Keep the final complete vibrational analysis if a log contains multiple links.
    freq_start = text.rfind("Harmonic frequencies")
    freq_text = text[freq_start:] if freq_start >= 0 else text
    for line in freq_text.splitlines():
        if re.match(r"\s*Frequencies\s+--", line):
            a.frequencies_cm1.extend(number(v) for v in re.findall(FLOAT, line.split("--", 1)[1]))
        elif re.match(r"\s*IR Inten\s+--", line):
            a.ir_intensities.extend(number(v) for v in re.findall(FLOAT, line.split("--", 1)[1]))
    if any(f < 0 for f in a.frequencies_cm1):
        a.warnings.append("存在虚频；需判断是否为目标过渡态或未达到极小点。")
    orientations = re.findall(r"(?:Standard|Input) orientation:\s*\n\s*-+\s*\n.*?\n.*?\n\s*-+\s*\n(.*?)\n\s*-+", text, re.S)
    if orientations:
        periodic = Chem.GetPeriodicTable()
        for line in orientations[-1].splitlines():
            tokens = line.split()
            if len(tokens) == 6:
                z = int(tokens[1])
                a.atoms.append({"index": int(tokens[0]), "atomic_number": z,
                                "symbol": periodic.GetElementSymbol(z) if z > 0 else "Bq",
                                "xyz": [number(v) for v in tokens[3:]]})
    charge_blocks = re.findall(r"Mulliken (?:atomic )?charges(?: and spin densities)?:\s*\n(.*?)(?:Sum of Mulliken|Electronic spatial extent)", text, re.S)
    if charge_blocks:
        for line in charge_blocks[-1].splitlines():
            row = re.match(r"\s*(\d+)\s+([A-Za-z]+)\s+(" + FLOAT + r")(?:\s+(" + FLOAT + "))?", line)
            if row:
                a.mulliken.append({"index": int(row[1]), "symbol": row[2], "charge": number(row[3]),
                                   "spin": number(row[4]) if row[4] else None})
    a.orbitals = _gaussian_orbitals(text)
    if a.orbitals:
        a.warnings.append("轨道能隙是当前计算方法的 HOMO–LUMO 能隙，不直接等于光学带隙或实验电离能。")
    return a


def parse_fchk(path):
    """Read selected fields without loading density matrices / MO coefficients."""
    wanted = {"Number of atoms", "Charge", "Multiplicity", "Number of alpha electrons",
              "Number of beta electrons", "Atomic numbers", "Current cartesian coordinates",
              "Alpha Orbital Energies", "Beta Orbital Energies", "Total Energy"}
    lines = Path(path).read_text(encoding="ascii", errors="replace").splitlines()
    data = {}
    i = 2
    while i < len(lines):
        line = lines[i]
        i += 1
        label = line[:43].strip()
        if label not in wanted:
            continue
        field_text = line[43:].strip()
        array = re.search(r"N=\s*(\d+)", field_text)
        if array:
            count = int(array[1])
            values = []
            while len(values) < count and i < len(lines):
                values.extend(number(v) for v in lines[i].split())
                i += 1
            if len(values) != count:
                raise InputError(f"fchk 数据不完整：{label}")
            data[label] = values
        else:
            data[label] = number(field_text.split()[-1])
    if "Atomic numbers" not in data or "Current cartesian coordinates" not in data:
        raise InputError("fchk 中缺少原子或坐标数据。")
    data["orbitals"] = {}
    for spin in ("Alpha", "Beta"):
        energies = data.get(spin + " Orbital Energies", [])
        occ = int(data.get("Number of " + spin.lower() + " electrons", 0))
        if summary := orbital_summary(energies, occ):
            data["orbitals"][spin.lower()] = summary
    return data


def _plotting():
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none"})
    return plt


def write_orbital_levels(orbitals, output):
    """Write a shared energy-axis diagram and the complete orbital table."""
    plt = _plotting()
    fig, axes = plt.subplots(figsize=(7.2, 4.25), layout="constrained")
    rows, spin_labels = [], []
    for lane, (spin, data) in enumerate(orbitals.items()):
        energies = data["energies_hartree"]
        homo, lumo = data["homo_index"], data["lumo_index"]
        occupied = data["is_occupied"]
        near_homo = sorted((i for i, filled in enumerate(occupied, 1) if filled),
                           key=lambda i: energies[i - 1], reverse=True)[:3]
        near_lumo = sorted((i for i, filled in enumerate(occupied, 1) if not filled),
                           key=lambda i: energies[i - 1])[:3]
        for i, energy in enumerate(energies, 1):
            rows.append((spin, i, occupied[i - 1], energy, energy * HARTREE_EV))
        nearby = [i for i in set(near_homo + near_lumo) if i not in (homo, lumo)]
        axes.hlines([energies[i - 1] * HARTREE_EV for i in nearby], lane - .14, lane + .14,
                    color="#a6b4c3", linewidth=1.25)
        for label, i, offset, align, color in (("HOMO", homo, -.18, 'right', '#1e578c'),
                                               ("LUMO", lumo, .18, 'left', '#a23957')):
            if i is None:
                continue
            y = energies[i - 1] * HARTREE_EV
            axes.hlines(y, lane - .14, lane + .14, color=color, linewidth=3)
            axes.text(lane + offset, y, f"{label}\n{y:.4f}", color=color, weight="bold",
                      ha=align, va="center", fontsize=10)
        spin_name = {'alpha': 'α spin', 'beta': 'β spin', 'shared': 'Shared orbitals'}[spin]
        gap = data['gap_ev']
        spin_labels.append(spin_name + (f"\nΔE = {gap:.4f} eV" if gap is not None else "\nΔE unavailable"))
    axes.set(xticks=range(len(orbitals)), xticklabels=spin_labels,
             xlim=(-.54, len(orbitals) - .46), ylabel="Orbital energy (eV)", title="Frontier molecular orbitals")
    axes.margins(y=.25)
    axes.grid(axis="y", color="#e1e7ee", linewidth=.7)
    axes.set_axisbelow(True)
    axes.tick_params(axis='x', length=0, pad=8)
    axes.spines[["top", "right"]].set_visible(False)
    try:
        for ext in ("svg", "png"):
            fig.savefig(output / ("orbital_levels." + ext), dpi=240)
    finally:
        plt.close(fig)
    with (output / "orbitals.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["spin", "orbital_1based", "occupied", "energy_hartree", "energy_ev"])
        writer.writerows(rows)


def write_report(a: Analysis, output: str | Path):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "analysis.json").write_text(json.dumps(asdict(a), ensure_ascii=False, indent=2), encoding="utf-8")
    plots = []
    plt = _plotting()
    if a.atoms:
        xyz = f"{len(a.atoms)}\nFinal geometry; Angstrom\n"
        xyz += "\n".join(f"{r['symbol']} " + " ".join(f"{v:.10f}" for v in r["xyz"]) for r in a.atoms) + "\n"
        (output / "final_structure.xyz").write_text(xyz, encoding="ascii")
    if a.orbitals:
        write_orbital_levels(a.orbitals, output)
        plots.append(("HOMO / LUMO 能级", "orbital_levels.svg"))
    if a.mulliken:
        with (output / "mulliken_charges.csv").open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["index", "symbol", "charge", "spin"])
            writer.writeheader()
            writer.writerows(a.mulliken)
        if len(a.mulliken) <= 150:
            fig, ax = plt.subplots(figsize=(max(6, len(a.mulliken) * .2), 4))
            charges = [r["charge"] for r in a.mulliken]
            ax.bar(range(len(charges)), charges, color=["#b33b55" if q > 0 else "#2464a8" for q in charges])
            ax.set_xticks(range(len(charges)), [str(r["index"]) + r["symbol"] for r in a.mulliken], rotation=90)
            ax.set_ylabel("Mulliken charge / e")
            ax.axhline(0, color="#666666", lw=.7)
            fig.tight_layout()
            fig.savefig(output / "atomic_charges.svg")
            fig.savefig(output / "atomic_charges.png", dpi=180)
            plt.close(fig)
            plots.append(("Mulliken 原子电荷（依赖基组）", "atomic_charges.svg"))
    if a.frequencies_cm1:
        with (output / "frequencies.csv").open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["frequency_cm-1", "IR_intensity_km_mol"])
            for i, value in enumerate(a.frequencies_cm1):
                writer.writerow([value, a.ir_intensities[i] if i < len(a.ir_intensities) else ""])
        if len(a.ir_intensities) == len(a.frequencies_cm1):
            x = np.linspace(0, max(4000, max(a.frequencies_cm1) + 100), 3000)
            y = np.zeros_like(x)
            for f, intensity in zip(a.frequencies_cm1, a.ir_intensities):
                if f > 0:
                    y += intensity * 10**2 / ((x - f)**2 + 10**2)
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(x, y, color="#2464a8")
            ax.invert_xaxis()
            ax.set(xlabel="Wavenumber / cm$^{-1}$", ylabel="Broadened IR intensity / a.u.", title="Unscaled harmonic IR; Lorentzian HWHM = 10 cm$^{-1}$")
            fig.tight_layout()
            fig.savefig(output / "ir_spectrum.svg")
            fig.savefig(output / "ir_spectrum.png", dpi=180)
            plt.close(fig)
            plots.append(("谐振 IR 光谱（未乘频率校正因子）", "ir_spectrum.svg"))
    values = [("来源", a.path), ("正常结束", a.normal), ("计算设置", a.route),
              ("电荷 / 多重度", f"{a.charge} / {a.multiplicity}"), ("电子能 / Hartree", a.energy_hartree),
              ("能量类型", a.energy_kind), ("含零点能 / Hartree", a.zpe_corrected_hartree),
              ("焓 / Hartree", a.enthalpy_hartree), ("自由能 / Hartree", a.free_energy_hartree),
              ("温度 / K", a.temperature_k), ("压力 / atm", a.pressure_atm),
              ("偶极矩 / Debye", a.dipole_debye), ("优化完成标记", a.optimization_complete),
              ("虚频数量", sum(v < 0 for v in a.frequencies_cm1))]
    for spin, orbitals in a.orbitals.items():
        for key in ("homo_ev", "lumo_ev", "gap_ev"):
            values.append((spin + " " + key, orbitals[key]))
    document = "<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>量化分析报告</title>"
    document += "<style>body{font:16px system-ui;max-width:1040px;margin:40px auto;padding:0 24px;color:#193047}table{border-collapse:collapse;width:100%}td{padding:9px;border-bottom:1px solid #ddd;overflow-wrap:anywhere}td:first-child{width:210px}img{max-width:100%}.note{background:#fff3d9;padding:12px;white-space:pre-wrap}</style>"
    document += "<h1>量化结果分析</h1><table>" + "".join("<tr><td>" + html.escape(k) + "</td><td>" + html.escape(str(v) if v is not None else "未找到") + "</td></tr>" for k, v in values) + "</table>"
    document += "<p class='note'>" + html.escape("\n".join(a.warnings)) + "</p>"
    for label, image in plots:
        document += f"<h2>{html.escape(label)}</h2><img src='{image}' alt='{html.escape(label)}'>"
    document += "<p>完整数值见 analysis.json 和 CSV；空间密度及轨道需从同一任务的 chk/fchk 生成。</p></html>"
    (output / "report.html").write_text(document, encoding="utf-8")
    return output / "report.html"


@dataclass
class Cube:
    origin: np.ndarray
    axes: np.ndarray
    atoms: np.ndarray
    values: np.ndarray
    comments: tuple[str, str] = ("Gaussian workflow", "Scalar grid; atomic units")


def read_cube(path) -> Cube:
    with Path(path).open(encoding="ascii", errors="strict") as stream:
        comments = (stream.readline().rstrip(), stream.readline().rstrip())
        row = stream.readline().split()
        natoms = int(row[0])
        if natoms < 0 or len(row) > 4 and int(row[4]) != 1:
            raise InputError("只支持单标量 cube；请每个轨道单独导出，不使用多轨道 cube。")
        origin = np.array([number(v) for v in row[1:4]])
        counts, axes = [], []
        for _ in range(3):
            row = stream.readline().split()
            counts.append(int(row[0]))
            axes.append([number(v) for v in row[1:4]])
        if any(n <= 0 for n in counts) or math.prod(counts) > 30_000_000:
            raise InputError("cube 必须使用 Bohr 单位的正网格数，且最多 3000 万网格点。")
        atoms = np.array([[number(v) for v in stream.readline().split()[:5]] for _ in range(natoms)]).reshape((-1, 5))
        # Fortran output may use D exponents.
        values = np.fromstring(stream.read().replace("D", "E").replace("d", "e"), sep=" ")
        if len(values) != math.prod(counts) or not np.isfinite(values).all():
            raise InputError("cube 数值数量或有效性不符合网格头。")
        return Cube(origin, np.array(axes), atoms, values.reshape(counts), comments)


def write_cube(cube: Cube, path):
    with Path(path).open("w", encoding="ascii", newline="\n") as stream:
        stream.write("\n".join(cube.comments) + "\n")
        stream.write(f"{len(cube.atoms):5d}" + "".join(f"{v:13.6f}" for v in cube.origin) + "\n")
        for n, axis in zip(cube.values.shape, cube.axes):
            stream.write(f"{n:5d}" + "".join(f"{v:13.6f}" for v in axis) + "\n")
        for atom in cube.atoms:
            stream.write(f"{int(atom[0]):5d}" + "".join(f"{v:13.6f}" for v in atom[1:]) + "\n")
        flat = cube.values.ravel()
        for start in range(0, len(flat), 6):
            stream.write(" ".join(f"{v: .6E}" for v in flat[start:start+6]) + "\n")


def matching_grid(left: Cube, right: Cube):
    if left.values.shape != right.values.shape or not np.allclose(left.origin, right.origin, atol=2e-6, rtol=0) or not np.allclose(left.axes, right.axes, atol=2e-6, rtol=0):
        raise InputError("cube 网格原点、步长或点数不一致。需用共同网格、同一坐标系重新生成，不能直接相减。")


def difference_density(ab_path, a_path, b_path, output):
    ab, a, b = (read_cube(p) for p in (ab_path, a_path, b_path))
    matching_grid(ab, a)
    matching_grid(ab, b)
    # Fragment geometries must partition the complex at the same fixed positions.
    real = lambda c: sorted((int(r[0]), *np.round(r[2:5], 5)) for r in c.atoms if r[0] > 0)
    if real(ab) != sorted(real(a) + real(b)):
        raise InputError("片段 A / B 的真实原子及坐标必须恰好组成 AB；不接受单独优化或平移后的片段密度。")
    delta = Cube(ab.origin, ab.axes, ab.atoms, ab.values - a.values - b.values,
                 ("Delta rho = rho(AB) - rho(A) - rho(B)", "Positive: electron accumulation; e/bohr^3"))
    write_cube(delta, output)
    voxel = abs(float(np.linalg.det(ab.axes)))
    result = {"definition": "rho(AB) - rho(A) - rho(B)", "integrated_delta_electrons": float(delta.values.sum() * voxel),
              "positive_electrons": float(delta.values.clip(0).sum() * voxel),
              "negative_electrons": float(delta.values.clip(max=0).sum() * voxel),
              "sources": [str(Path(p).resolve()) for p in (ab_path, a_path, b_path)],
              "note": "电子密度差，非电荷密度；正值为电子积累。计算方法/基组需由使用者确认一致，有限网格积分存在截断误差。"}
    Path(output).with_suffix(".json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def energy_difference(ab: Analysis, fragments: list[Analysis], *, free_energy=False):
    all_results = [ab, *fragments]
    attribute = "free_energy_hartree" if free_energy else "energy_hartree"
    if not all(r.normal and getattr(r, attribute) is not None for r in all_results):
        raise InputError("必须提供正常结束且含所选能量的全部日志。自由能需要频率热化学数据。")
    if not ab.method_basis or any(r.method_basis != ab.method_basis or r.energy_kind != ab.energy_kind or r.solvent != ab.solvent for r in fragments):
        raise InputError("方法 / 基组、能量类型或溶剂不一致（或无法识别），不能直接比较。")
    if not ab.atoms or sorted(a["atomic_number"] for a in ab.atoms) != sorted(a["atomic_number"] for r in fragments for a in r.atoms):
        raise InputError("AB 与各片段的元素组成不守恒。")
    if ab.charge is None or any(r.charge is None for r in fragments) or ab.charge != sum(r.charge for r in fragments):
        raise InputError("AB 电荷必须等于片段电荷之和。")
    if free_energy:
        if any(not r.frequencies_cm1 or any(v < 0 for v in r.frequencies_cm1) for r in all_results):
            raise InputError("自由能比较要求所有体系有频率数据且无虚频。")
        if any(r.temperature_k != ab.temperature_k or r.pressure_atm != ab.pressure_atm for r in fragments):
            raise InputError("各日志热化学温度 / 压力不一致。")
    delta = getattr(ab, attribute) - sum(getattr(r, attribute) for r in fragments)
    return {"definition": "G(AB)-sum(G(fragment))" if free_energy else "E(AB)-sum(E(fragment))",
            "hartree": delta, "kj_mol": delta * HARTREE_KJMOL, "kcal_mol": delta * HARTREE_KJMOL / 4.184,
            "note": "负值表示该能量定义下结合有利。未校正 BSSE；孤立片段给结合能，冻结复合物几何片段给相互作用能。热自由能沿用日志标准态，未自动做 1 atm→1 M 修正。"}


def solvation_difference(gas: Analysis, solution: Analysis, *, free_energy=False):
    attribute = "free_energy_hartree" if free_energy else "energy_hartree"
    if not gas.normal or not solution.normal or getattr(gas, attribute) is None or getattr(solution, attribute) is None:
        raise InputError("气相与溶液日志必须正常结束并含所需能量。")
    if gas.solvent or not solution.solvent:
        raise InputError("第一份日志应为气相，第二份应含可识别的溶剂设置。")
    if not gas.method_basis or gas.method_basis != solution.method_basis or gas.energy_kind != solution.energy_kind or gas.charge != solution.charge or gas.multiplicity != solution.multiplicity:
        raise InputError("气相 / 溶液的方法、基组、能量类型、电荷和多重度必须一致。")
    if sorted(a["atomic_number"] for a in gas.atoms) != sorted(a["atomic_number"] for a in solution.atoms) or not gas.atoms:
        raise InputError("两体系元素组成不一致或缺少坐标。")
    if free_energy and (not gas.frequencies_cm1 or not solution.frequencies_cm1 or
                        any(v < 0 for r in (gas, solution) for v in r.frequencies_cm1) or
                        gas.temperature_k != solution.temperature_k or gas.pressure_atm != solution.pressure_atm):
        raise InputError("自由能差需要两份无虚频、温压一致的频率热化学结果。")
    delta = getattr(solution, attribute) - getattr(gas, attribute)
    return {"definition": ("G" if free_energy else "E") + "(solution)-" + ("G" if free_energy else "E") + "(gas)",
            "hartree": delta, "kj_mol": delta * HARTREE_KJMOL, "kcal_mol": delta * HARTREE_KJMOL / 4.184,
            "solvent": solution.solvent,
            "note": "此值为日志所给总自由能/电子能之差；是否为垂直或弛豫差取决于输入几何。未自动加入浓度标准态修正，电子能差不等于完整实验溶剂化自由能。"}
