"""Read small-molecule CIF blocks, expand symmetry and unwrap periodic bonds."""
from pathlib import Path
import re

import gemmi
import numpy as np
from ase import Atoms
from ase.geometry import cellpar_to_cell
from ase.neighborlist import natural_cutoffs, neighbor_list

from gaussian_core import InputError, Result


def blocks(path):
    try:
        doc = gemmi.cif.read_file(str(Path(path)))
    except (OSError, RuntimeError, ValueError) as exc:
        raise InputError(f"无法读取 CIF：{exc}") from exc
    result = {b.name: b for b in doc if any(b.find_values('_atom_site_' + k + '_x') for k in ('fract', 'Cartn'))}
    if not result:
        raise InputError('未找到小分子 CIF 原子坐标表；蛋白质 mmCIF 不在此导入器范围内。')
    return result


def read_cell(block):
    def column(tag):
        return list(block.find_values(tag))

    def numeric(raw):
        value = gemmi.cif.as_number(raw or '?')
        if not np.isfinite(value):
            raise InputError(f'CIF 数值缺失或无效：{raw}')
        return value

    params = [numeric(block.find_value('_cell_' + key)) for key in
              ('length_a', 'length_b', 'length_c', 'angle_alpha', 'angle_beta', 'angle_gamma')]
    try:
        if min(params[:3]) <= 0 or any(not 0 < a < 180 for a in params[3:]):
            raise ValueError()
        cell = cellpar_to_cell(params)
        if np.linalg.det(cell) < 1e-6:
            raise ValueError()
    except (ValueError, AssertionError) as exc:
        raise InputError('CIF 晶胞长度或角度无效。') from exc
    ops = column('_space_group_symop_operation_xyz') or column('_symmetry_equiv_pos_as_xyz')
    if ops:
        operations = [gemmi.Op(gemmi.cif.as_string(op)) for op in ops]
    else:
        structure = gemmi.make_small_structure_from_block(block)
        if not structure.spacegroup:
            raise InputError('CIF 缺少可识别的空间群或对称操作，无法展开晶胞。')
        operations = list(structure.spacegroup.operations())
    labels = [gemmi.cif.as_string(x) for x in column('_atom_site_label')]
    species = [gemmi.cif.as_string(x) for x in column('_atom_site_type_symbol')] or labels
    fractional = bool(column('_atom_site_fract_x'))
    prefix = '_atom_site_fract_' if fractional else '_atom_site_Cartn_'
    axes = [column(prefix + axis) for axis in 'xyz']
    if not species or any(len(axis) != len(species) for axis in axes):
        raise InputError('CIF 元素或坐标列不完整。')
    occupancy = column('_atom_site_occupancy')
    disorder = column('_atom_site_disorder_group')
    if (occupancy and (len(occupancy) != len(species) or any(abs(numeric(x)-1) > 1e-6 for x in occupancy))) or any(x not in ('.','?','0') for x in disorder):
        raise InputError('存在无序或部分占位，请先选定一个有序构型再导入。')
    positions, symbols, expanded_labels = [], [], []
    for i, raw in enumerate(species):
        match = re.fullmatch(r'([A-Za-z]{1,2})(?:[\d+\-].*)?', raw)
        element = gemmi.Element(match[1]) if match else None
        if not element or not element.atomic_number or element.name == 'D':
            raise InputError(f'未支持的元素或同位素标记：{raw}')
        frac = np.array([numeric(axis[i]) for axis in axes]); site = []
        if not fractional:
            frac = np.linalg.solve(cell.T, frac)
        for op in operations:
            pos = np.array(op.apply_to_xyz(frac)) % 1
            if site:
                delta = pos - np.asarray(site)
                if np.min(np.linalg.norm((delta - np.rint(delta)) @ cell, axis=1)) < 0.01:
                    continue  # symmetry-equivalent copy of the same special position
            site.append(pos)
            positions.append(pos); symbols.append(element.name)
            expanded_labels.append(labels[i] if labels else str(i+1))
    atoms = Atoms(symbols, scaled_positions=positions, cell=cell, pbc=True)
    if len(neighbor_list('d', atoms, 0.1)):
        raise InputError('不同原子位点重叠，请检查 CIF 无序、占位及对称操作。')
    return atoms, expanded_labels


def components(atoms):
    indices, neighbors, offsets = neighbor_list('ijS', atoms, natural_cutoffs(atoms, mult=1.18))
    graph = [[] for _ in atoms]
    for i, j, shift in zip(indices, neighbors, offsets):
        graph[i].append((j, shift))
    found = []
    seen = set()
    for i in range(len(atoms)):
        if i in seen:
            continue
        images = {i: np.zeros(3, dtype=int)}; stack = [i]; periodic = False
        while stack:
            a = stack.pop()
            for b, shift in graph[a]:
                offset = images[a] + shift
                if b not in images:
                    images[b] = offset; stack.append(b)
                elif not np.array_equal(images[b], offset):
                    periodic = True
        ids = sorted(images); seen.update(ids)
        fragment = atoms[ids]
        fragment.positions += np.asarray([images[j] for j in ids]) @ atoms.cell.array
        fragment.pbc = False
        found.append((fragment, periodic, ids))
    return found


def describe(parts, labels):
    return [f'{i}. {f.get_chemical_formula()} · {len(f)} 原子 · {labels[ids[0]]}' + (' · 周期网络' if net else '')
            for i, (f, net, ids) in enumerate(parts, 1)]


def generate_cif(value, s):
    options = blocks(value)
    key = s.cif_block or next(iter(options))
    if key not in options or (not s.cif_block and len(options) > 1):
        raise InputError('请选择 CIF 数据块：' + '、'.join(options))
    atoms, labels = read_cell(options[key])
    periodic = s.cif_scope == 'periodic'
    if s.cif_scope not in ('molecule', 'periodic'):
        raise InputError('未知 CIF 模型。')
    if periodic and (s.engine != 'gaussian' or s.job != 'sp'):
        raise InputError('周期晶胞请选择 Gaussian 单点能，ORCA 可使用有限分子。')
    notes = ['保留 CIF 坐标和已有氢。未补氢、未作力场预优化或硼平面化；请核对实验结构的氢、总电荷及自旋。']
    if not periodic:
        parts = components(atoms)
        if not 1 <= s.cif_component <= len(parts):
            raise InputError('请选择有效的片段编号：' + '；'.join(describe(parts, labels)[:30]))
        atoms, net, _ = parts[s.cif_component-1]
        if net:
            raise InputError('所选片段沿晶格无限连接，请选择周期晶胞，或导入已截取并封端的簇。')
        notes.append(f'按共价半径和跨胞距离提取第 {s.cif_component}/{len(parts)} 个连通片段；金属配位及弱连接请核对。')
    charge = 0 if s.charge is None else s.charge
    if s.charge is None:
        notes.append('电荷按中性 0 处理；提取离子时请手动指定。')
    if s.multiplicity is None and any(21 <= z <= 30 or 39 <= z <= 48 or 57 <= z <= 80 or z >= 89 for z in atoms.numbers):
        raise InputError('含金属的 CIF 请手动填写自旋多重度。')
    electrons = int(sum(atoms.numbers)) - charge
    mult = s.multiplicity or (1 + electrons % 2)
    if electrons < 1 or mult > electrons + 1 or (electrons - mult + 1) % 2:
        raise InputError('CIF 模型的电荷与自旋多重度不相容。')
    job = {'optfreq':'Opt Freq','opt':'Opt','freq':'Freq','sp':'SP'}[s.job]
    lines = [f'%chk={s.checkpoint}.chk', f'%mem={s.memory_gb}GB', f'%nprocshared={s.processors}',
             f'#p {s.method}/{s.basis} {job} {s.extra_keywords}'.rstrip(), '', s.title, '', f'{charge} {mult}']
    coords = list(zip(atoms.get_chemical_symbols(), atoms.positions))
    if periodic:
        coords.extend(('Tv', vector) for vector in atoms.cell.array)
    lines.extend(f'{symbol:<4} {x: .8f} {y: .8f} {z: .8f}' for symbol,(x,y,z) in coords)
    return Result('\n'.join(lines)+'\n\n', '', atoms.get_chemical_formula(), charge, mult, len(atoms),
                  'CIF 周期晶胞' if periodic else 'CIF 展开对称并还原跨胞分子', None, f'{Path(value).name} / {key}', notes)
def choose_crystal(parent, variables):
    from tkinter import ttk, filedialog, messagebox, StringVar, Toplevel
    path = filedialog.askopenfilename(parent=parent, title='导入小分子 CIF', filetypes=[('晶体文件', '*.cif')])
    if not path:
        return
    try:
        data = blocks(path)
    except Exception as exc:
        messagebox.showerror('CIF 导入失败', str(exc), parent=parent)
        return
    window = Toplevel(parent)
    window.title('选择 CIF 结构'); window.transient(parent); window.grab_set()
    panel = ttk.Frame(window, padding=16); panel.pack(fill='both', expand=True)
    block = StringVar(window, next(iter(data))); model = StringVar(window); note = StringVar(window)
    ttk.Label(panel,text='数据块').grid(row=0,column=0,sticky='w')
    selector = ttk.Combobox(panel,textvariable=block,values=list(data),state='readonly',width=44)
    selector.grid(row=0,column=1,pady=5)
    ttk.Label(panel,text='计算模型').grid(row=1,column=0,sticky='w')
    models = ttk.Combobox(panel,textvariable=model,state='readonly',width=44); models.grid(row=1,column=1,pady=5)
    ttk.Label(panel,textvariable=note,wraplength=370).grid(row=2,columnspan=2,pady=6,sticky='w')
    def update(*_):
        try:
            atoms,labels = read_cell(data[block.get()]); parts = components(atoms)
            items = [label for label, (_,net,_) in zip(describe(parts,labels),parts) if not net]
            models.configure(values=[*items,'完整周期晶胞（Gaussian 单点）']); models.current(0)
            note.set(f'完整晶胞 {atoms.get_chemical_formula()}，共 {len(atoms)} 个原子。请核对氢、电荷和多重度。'); ok.configure(state='normal')
        except Exception as exc:
            model.set(''); models.configure(values=()); note.set(str(exc)); ok.configure(state='disabled')
    def accept():
        periodic = model.get().startswith('完整周期')
        variables['cif_scope'].set('periodic' if periodic else 'molecule')
        variables['cif_component'].set(1 if periodic else int(model.get().split('.')[0]))
        variables['cif_block'].set(block.get()); variables['input'].set(path); variables['mode'].set('晶体 CIF')
        if periodic:
            variables['engine'].set('Gaussian'); variables['job'].set('单点能')
        window.destroy()
    ok = ttk.Button(panel,text='使用此结构',command=accept); ok.grid(row=3,column=1,sticky='e',pady=5)
    selector.bind('<<ComboboxSelected>>',update); update()
