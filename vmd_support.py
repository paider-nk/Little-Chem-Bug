"""Prepare a result and start VMD with a verified molecular representation."""
from pathlib import Path
import os
import math
import shutil
import external_process as subprocess
import tempfile
import time

from analysis_core import parse_log
from gaussian_core import InputError


def locate_vmd(explicit=None):
    paths = [explicit] if explicit else [os.environ.get('VMD_EXE'), shutil.which('vmd')]
    for candidate in filter(None, paths):
        path = Path(candidate)
        if path.is_dir():
            path = path / 'vmd.exe'
        if path.is_file():
            return str(path.resolve())
    raise InputError('找不到 VMD，请在可视化页选择 vmd.exe。')


def write_xyz_from_log(source, output=None):
    source = Path(source).resolve()
    a = parse_log(source)
    atoms = [atom for atom in a.atoms if atom['atomic_number'] > 0]
    if not atoms:
        raise InputError('日志中没有可显示的原子坐标。')
    target = Path(output).resolve() if output else source.with_suffix('.final.xyz')
    if target.suffix.lower() != '.xyz' or target == source:
        raise InputError('结构导出的目标文件必须是单独的 .xyz 文件。')
    lines = [str(len(atoms)), 'Final geometry (Angstrom)' if a.normal else 'Last recorded geometry (incomplete job)']
    for atom in atoms:
        x, y, z = atom['xyz']
        lines.append(f"{atom['symbol']:<3} {x: .10f} {y: .10f} {z: .10f}")
    with tempfile.NamedTemporaryFile(dir=target.parent, delete=False, mode='w', encoding='ascii', newline='\n') as tmp:
        tmp.write('\n'.join(lines)+'\n')
    try:
        os.replace(tmp.name, target)
    finally:
        Path(tmp.name).unlink(missing_ok=True)
    return target


def prepare_vmd_input(source, output=None):
    path = Path(source).resolve()
    if not path.is_file():
        raise InputError(f'文件不存在：{path}')
    suffix = path.suffix.lower()
    if suffix in ('.log', '.out'):
        return write_xyz_from_log(path, output), 'xyz'
    if suffix == '.gbw':
        from orca_support import gbw_to_molden
        return gbw_to_molden(path), 'molden'
    formats = {'.xyz':'xyz', '.pdb':'pdb', '.cube':'cube', '.cub':'cube', '.molden':'molden'}
    if path.name.lower().endswith('.molden.input'):
        return path, 'molden'
    if suffix in formats:
        return path, formats[suffix]
    raise InputError('请选择 .log、.out、.xyz、.pdb、.cube、.cub、.gbw 或 Molden 波函数。')


def vmd_environment(executable):
    env = os.environ.copy()
    bindir = Path(executable).parent
    for root in (bindir, bindir.parent):
        if (root/'scripts/vmd/vmdinit.tcl').is_file():
            env.update(VMDDIR=str(root), TCL_LIBRARY=str(root/'scripts/tcl8.6'), TK_LIBRARY=str(root/'scripts/tk8.6'))
            break
    env['PATH'] = str(bindir) + os.pathsep + env.get('PATH','')
    return env


def open_vmd(source, *, vmd=None, output=None, isovalue=0.01):
    executable = locate_vmd(vmd)
    path, kind = prepare_vmd_input(source, output)
    if kind == 'cube' and not (math.isfinite(isovalue) and isovalue > 0):
        raise InputError('等值面阈值必须是有限正数。')
    env = vmd_environment(executable)
    script = Path(__file__).resolve().parent/'assets/workflow_vmd.tcl'
    with tempfile.TemporaryDirectory(prefix='vmd-result-') as temporary:
        workspace = Path(temporary)
        status = workspace/'loaded.txt'
        env.update(WORKFLOW_VMD_FILE=path.as_posix().encode().hex(), WORKFLOW_VMD_TYPE=kind,
                   WORKFLOW_VMD_ISO=str(isovalue), WORKFLOW_VMD_STATUS=status.as_posix().encode().hex())
        # Local ASCII script name avoids VMD's legacy Windows command-line path encoding.
        proc = subprocess.Popen([executable, '-e', script.name], cwd=script.parent, env=env)
        deadline = time.monotonic()+120
        while not status.exists():
            if proc.poll() is not None or time.monotonic() > deadline:
                raise InputError('VMD 未确认加载结果。请检查 VMD 窗口和脚本、插件安装目录。')
            time.sleep(0.1)
        code, detail = status.read_text(encoding='utf-8').strip().split(' ',1)
        if code != '0':
            raise InputError('VMD 无法显示此文件：'+detail)
    return {'opened': str(path), 'pid': proc.pid, 'atoms': int(detail.split()[0])}
