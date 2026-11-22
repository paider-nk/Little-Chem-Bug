"""Export a frontier orbital analysis from a calculation result."""
from pathlib import Path
import csv
import re
import shutil
import tempfile

from analysis_core import number, orbital_summary, parse_fchk, parse_log, write_orbital_levels
from gaussian_core import InputError
from gaussian_workflow import to_fchk
from orca_support import gbw_to_molden

SPIN_LABELS = {"alpha": "α", "beta": "β", "shared": "单组轨道"}


def read_molden_orbitals(path):
    records = {"alpha": [], "beta": []}
    fields = {}
    in_mo = False
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            line = line.strip()
            if line.startswith('['):
                if in_mo:
                    break
                in_mo = line.casefold() == '[mo]'
            if not in_mo:
                continue
            match = re.fullmatch(r'(Ene|Spin|Occup)\s*=\s*(\S+)', line, re.I)
            if not match:
                continue
            key = match[1].lower()
            if key in fields:
                raise InputError('Molden 轨道记录缺少能量、自旋或占据数。')
            fields[key] = match[2]
            if len(fields) == 3:
                spin = fields['spin'].lower()
                if spin not in records:
                    raise InputError('Molden 文件含无法识别的自旋标记。')
                records[spin].append((number(fields['ene']), number(fields['occup'])))
                fields.clear()
    if fields:
        raise InputError('Molden 轨道记录不完整。')
    return {spin: orbital_summary([energy for energy, _ in rows], 0,
                                  occupations=[occupation for _, occupation in rows])
            for spin, rows in records.items() if rows}


def generate_frontier_levels(source, output=None, *, gaussian=None, orca=None, on_stage=lambda _: None):
    source = Path(source).resolve()
    if not source.is_file():
        raise InputError(f'找不到计算结果：{source}')
    suffix = '.molden' if source.name.lower().endswith('.molden.input') else source.suffix.lower()
    on_stage('读取轨道能量与占据数')
    warnings = []
    if suffix in ('.log', '.out'):
        analysis = parse_log(source)
        orbitals = analysis.orbitals
        if not analysis.normal:
            warnings.append('日志未正常结束，以下是最后输出的轨道数据。')
    elif suffix in ('.fchk', '.fch'):
        orbitals = parse_fchk(source)['orbitals']
    elif suffix == '.molden':
        orbitals = read_molden_orbitals(source)
    elif suffix in ('.chk', '.gbw'):
        with tempfile.TemporaryDirectory(prefix='frontier-') as temporary:
            staged = Path(temporary) / ('wave' + suffix)
            shutil.copyfile(source, staged)
            if suffix == '.chk':
                orbitals = parse_fchk(to_fchk(staged, gaussian=gaussian))['orbitals']
            else:
                orbitals = read_molden_orbitals(gbw_to_molden(staged, orca=orca))
    else:
        raise InputError('请选择 Gaussian .log/.chk/.fchk 或 ORCA .out/.gbw/Molden 文件。')
    orbitals = {spin: data for spin, data in orbitals.items() if data}
    if not orbitals:
        raise InputError('未找到轨道能量。请选择含轨道表的日志，或使用同一任务的波函数文件。')
    if set(orbitals) == {'alpha'}:
        orbitals = {'shared': orbitals['alpha']}
    if any(data['homo_index'] is None or data['lumo_index'] is None for data in orbitals.values()):
        warnings.append('部分通道缺少占据或空轨道，对应能级和能隙留空。')
    on_stage('绘制轨道能级图并导出数值')
    parent = Path(output).resolve() if output else source.parent
    parent.mkdir(parents=True, exist_ok=True)
    folder = Path(tempfile.mkdtemp(prefix=f'{source.stem[:60]}_HOMO-LUMO_', dir=parent))
    write_orbital_levels(orbitals, folder)
    rows = [{"spin": spin, **{key: data[key] for key in ('homo_index', 'lumo_index', 'homo_ev', 'lumo_ev', 'gap_ev')}}
            for spin, data in orbitals.items()]
    summary = folder / 'frontier_summary.csv'
    with summary.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['source', 'status', *rows[0]])
        writer.writeheader()
        for row in rows:
            writer.writerow({'source': str(source), 'status': '；'.join(warnings), **row})
    return {'folder': folder, 'png': folder / 'orbital_levels.png', 'summary': summary,
            'rows': rows, 'warning': ' '.join(warnings)}
