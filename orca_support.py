"""ORCA inputs, execution and wavefunction conversion."""
from pathlib import Path
import os
import re
import shutil
import external_process as subprocess
import tempfile

from gaussian_core import InputError


def locate_orca(explicit=None):
    candidate = explicit or os.environ.get('ORCA_EXE') or shutil.which('orca')
    if candidate:
        path = Path(candidate)
        if path.is_dir():
            path = path / 'orca.exe'
        if path.is_file():
            return str(path.resolve())
    raise InputError('未找到 ORCA，请选择 orca.exe 或设置 ORCA_EXE。')


def resources(cores, memory):
    return f'%pal nprocs {cores} end\n%maxcore {max(1, int(memory * 1024 * 0.75 / cores))}'


def render(result, s):
    sections = result.text.strip().split('\n\n')
    geom = sections[2].splitlines()[1:]
    if any('(' in line.split()[0] for line in geom):
        raise InputError('ORCA 导出尚不支持同位素指定。')
    method = {'pbe1pbe':'PBE0','m062x':'M06-2X'}.get(s.method.lower(), s.method)
    if s.method.lower() == 'wb97xd':
        raise InputError('wB97XD 两端命名不同，请直接填写已确认的 ORCA 泛函名。')
    basis = {'def2svp':'def2-SVP','def2tzvp':'def2-TZVP','6-31g(d)':'6-31G*'}.get(s.basis.lower(), s.basis)
    extra = s.extra_keywords.strip()
    if extra == 'SCF=XQC Int=UltraFine':
        extra = 'TightSCF'
    if not re.fullmatch(r'[A-Za-z0-9 ()+*.,/\-]*', extra):
        raise InputError('请使用 ORCA 的 ! 行关键词，例如 TightSCF CPCM(Water)。')
    job = {'opt':'Opt','freq':'Freq','optfreq':'Opt Freq','sp':'SP'}[s.job]
    lines = [f'# {s.title}',f'! {method} {basis} {job} {extra}'.rstrip(),resources(s.processors,s.memory_gb)]
    constraints = []
    for tail in sections[3:]:
        for line in tail.splitlines():
            m = re.fullmatch(r'D (\d+) (\d+) (\d+) (\d+) ([-\d.]+)', line)
            if m:
                ids = ' '.join(str(int(x)-1) for x in m.groups()[:4])
                constraints.append('    { D '+ids+' '+m[5]+' C }')
    if constraints:
        lines += ['%geom','  Constraints',*constraints,'  end','end']
    lines += [f'* xyz {result.charge} {result.multiplicity}', *geom, '*']
    if 'Geom=Connectivity' in sections[0]:
        result.warnings.append('ORCA 使用全部原子的 XYZ 坐标；Gaussian 连通表和配位键级仅用于 Gaussian 导出。')
    result.text = '\n'.join(lines)+'\n'
    return result


def run_one(executable, source, root, cores, memory):
    # A fresh folder prevents reuse of an unrelated GBW through ORCA AutoStart.
    folder = Path(tempfile.mkdtemp(prefix=source.stem+'-', dir=root))
    inp, out = folder/'job.inp', folder/'job.out'
    text = source.read_text(encoding='utf-8-sig')
    if (not re.search(r'^\s*\*\s+xyz\s+[-+]?\d+\s+\d+\s*$', text, re.M|re.I) or
        re.search(r'(?im)^\s*(?:\*\s+xyzfile|%\s*(?:moinp|base|include|coords)|\$new_job|\$compound)|\bmoread\b', text)):
        raise InputError('批处理需要内嵌 * xyz 坐标的独立 ORCA 输入文件。')
    text = re.sub(r'(?ims)^\s*%pal\b.*?\bend\b', '', text)
    text = re.sub(r'(?im)^\s*%maxcore\s+\d+\s*$', '', text)
    text = re.sub(r'(?im)^([ \t]*!.*)$', lambda m: re.sub(r'\bPAL\d+\b','',m[0],flags=re.I), text)
    inp.write_text(resources(cores,memory)+'\n'+text, encoding='utf-8')
    with out.open('wb') as stream:
        done = subprocess.run([executable, inp.name], cwd=folder, stdout=stream, stderr=subprocess.STDOUT)
    normal = done.returncode == 0 and 'ORCA TERMINATED NORMALLY' in out.read_text(errors='replace')[-3000:]
    return dict(input=str(source), log=str(out), return_code=done.returncode, normal_termination=normal, status='normal' if normal else 'incomplete')


def gbw_to_molden(path, *, orca=None):
    utility = Path(locate_orca(orca)).with_name('orca_2mkl.exe')
    if not utility.is_file():
        raise InputError('ORCA 目录中未找到 orca_2mkl.exe。')
    with tempfile.TemporaryDirectory(dir=path.parent) as tmp:
        scratch = Path(tmp)
        shutil.copyfile(path,scratch/'wave.gbw')
        done = subprocess.run([str(utility),'wave','-molden'],cwd=scratch,capture_output=True)
        converted = scratch/'wave.molden.input'
        if done.returncode or not converted.is_file():
            raise InputError('无法将 ORCA GBW 转换为 Molden 波函数。')
        target = path.with_suffix('.molden.input')
        os.replace(converted,target)
        return target
