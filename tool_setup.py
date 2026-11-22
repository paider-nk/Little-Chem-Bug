"""Per-user settings and Windows tool installation for the desktop release."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import subprocess
import tempfile
from urllib.request import Request, urlopen
from zipfile import ZipFile

from external_process import run

# Environment variables are the public interface shared with the CLI backends.
TOOLS = {
    'GAUSSIAN_EXE': ('Gaussian 16', 'g16.exe', 'https://gaussian.com/g16w/', '运行 Gaussian 计算'),
    'FORMCHK_EXE': ('Gaussian formchk', 'formchk.exe', '', '转换 Gaussian chk，通常和 g16.exe 在同一目录'),
    'ORCA_EXE': ('ORCA', 'orca.exe', 'https://www.faccts.de/orca/', '运行 ORCA 计算'),
    'MULTIWFN_EXE': ('Multiwfn', 'Multiwfn.exe', 'https://www.umsyar.com/multiwfn/download.html', '波函数、静电势和密度分析'),
    'VMD_EXE': ('VMD', 'vmd.exe', 'https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD', '查看分子和渲染静电势表面'),
    'XTB_EXE': ('xTB', 'xtb.exe', 'https://github.com/grimme-lab/xtb/releases/tag/v6.8.0', 'GFN-xTB、tblite 和 GFN-FF 预处理'),
    'MOLSIMPLIFY_PYTHON': ('molSimplify', 'python.exe', 'https://pypi.org/project/molSimplify/2.0.0/#files', '按配位拓扑和供体构建配合物'),
}
XTB_URL = 'https://github.com/grimme-lab/xtb/releases/download/v6.8.0/xtb-6.8.0-windows-x86_64-openblas.zip'
XTB_HASH = 'abd7f74d60908925a7e88a030f83d448253ff4e087451f9b13aeda755eb7b98d'
UV_URL = ('https://files.pythonhosted.org/packages/06/a8/d72229dd90d1e5a3c8368d51a70219018d579380945e67c8dcffbe8e53c0/'
          'uv-0.9.28-py3-none-win_amd64.whl')
UV_HASH = '203ab59710c0c1b3c5ecc684f9cfc9264340a69c8706aaa8aea75415779f0d74'
MOLSIMPLIFY_WHEEL = 'molsimplify-2.0.0-py3-none-any.whl'
MOLSIMPLIFY_URL = ('https://files.pythonhosted.org/packages/77/4e/acdd1152abf168be1ddd244b7fcbf7c35992a87419c7f87fe10ddc5852ca/'
                  + MOLSIMPLIFY_WHEEL)
MOLSIMPLIFY_HASH = 'b2d5e63efae66ef5ee7b03168d54daaabc8e52eee6cd9137fb4430683e1e8d8d'
PACKAGE_ROOT = Path(__file__).resolve().parent
PROBE_MOLSIMPLIFY = 'from molSimplify.Scripts.generator import startgen_pythonic; from openbabel import openbabel'


def data_dir():
    override = os.environ.get('CHEMWORKBENCH_HOME')
    return Path(override).resolve() if override else Path(os.environ.get('LOCALAPPDATA') or Path.home()) / 'ChemWorkbench'


def read_settings():
    try:
        raw = json.loads((data_dir() / 'tools.json').read_text(encoding='utf-8'))
        return {key: raw[key] for key in TOOLS if isinstance(raw.get(key), str)}
    except (FileNotFoundError, json.JSONDecodeError, AttributeError):
        return {}


class ToolSettings:
    def __init__(self):
        self.inherited = {key: os.environ.get(key, '') for key in TOOLS}
        self.paths = read_settings()
        self.apply()

    def apply(self):
        for key in TOOLS:
            path = self.paths.get(key, self.inherited[key])
            if path:
                os.environ[key] = path
            else:
                os.environ.pop(key, None)

    def save(self, key, path):
        updated = {**self.paths, key: str(path)}
        folder = data_dir()
        folder.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=folder, delete=False) as stream:
            tmp = Path(stream.name)
            json.dump(updated, stream, ensure_ascii=False, indent=2)
        try:
            tmp.replace(folder / 'tools.json')
        finally:
            tmp.unlink(missing_ok=True)
        self.paths = updated
        self.apply()


def download(url, target, digest, report):
    report('正在下载并校验官方软件包…')
    system_dir = Path(os.environ.get('SystemRoot', r'C:\Windows')) / 'System32'
    curl = system_dir / 'curl.exe'
    if os.name == 'nt' and curl.is_file():
        checked_run([curl, '--fail', '--location', '--silent', '--show-error',
                     '--proto', '=https', '--proto-redir', '=https', '--connect-timeout', '30',
                     '--max-time', '300', '--output', target, url], timeout=310)
    else:
        request = Request(url, headers={'User-Agent': 'ChemWorkbench/1.0'})
        with urlopen(request, timeout=45) as source, target.open('wb') as output:
            shutil.copyfileobj(source, output, length=1024 * 1024)
    with target.open('rb') as stream:
        checksum = hashlib.file_digest(stream, 'sha256')
    if checksum.hexdigest() != digest:
        raise ValueError('下载文件的 SHA-256 不匹配。请重试。')


def unpack(archive, target):
    """Extract a vendor ZIP after validating Windows paths and entry types."""
    target.mkdir(parents=True, exist_ok=True)
    with ZipFile(archive) as source:
        for item in source.infolist():
            path = PurePosixPath(item.filename.replace('\\', '/'))
            if (not path.parts or path.is_absolute() or '..' in path.parts
                    or any(os.path.isreserved(part) if hasattr(os.path, 'isreserved')
                           else ':' in part or part.endswith((' ', '.')) for part in path.parts)
                    or stat.S_ISLNK(item.external_attr >> 16)):
                raise ValueError(f'压缩包包含不支持的路径：{item.filename}')
        source.extractall(target)


def find_in_folder(folder, key):
    name = TOOLS[key][1]
    folder = Path(folder).resolve()
    candidates = (folder / name, folder / 'bin' / name, folder / 'Scripts' / name)
    for path in candidates:
        if path.is_file():
            return path
    matches = sorted((p for p in folder.rglob(name) if p.is_file()), key=lambda p: (len(p.parts), str(p).casefold()))
    if matches:
        return matches[0]
    raise ValueError(f'所选目录中找不到 {name}。')


def detect_tool(key):
    name = TOOLS[key][1]
    path = os.environ.get(key, '')
    if path and Path(path).is_file():
        return Path(path).resolve()
    if key == 'FORMCHK_EXE':
        gaussian = detect_tool('GAUSSIAN_EXE')
        if gaussian and gaussian.with_name(name).is_file():
            return gaussian.with_name(name)
    if key == 'MOLSIMPLIFY_PYTHON':
        path = PACKAGE_ROOT / '.molsimplify' / 'Scripts' / name
        return path if path.is_file() else None
    found = shutil.which(name)
    if found:
        return Path(found).resolve()
    return None


def tool_run(args, *, timeout=120, env=None):
    return run([str(arg) for arg in args], capture_output=True, text=True,
               encoding='utf-8', errors='replace', timeout=timeout, env=env,
               creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)


def checked_run(args, *, timeout=120, env=None):
    try:
        result = tool_run(args, timeout=timeout, env=env)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f'当前步骤运行超过 {timeout} 秒。若正在下载，请检查网络后重试。') from exc
    if result.returncode:
        raise RuntimeError((result.stderr.strip() or result.stdout.strip() or f'退出代码 {result.returncode}')[-2500:])
    return result.stdout


def install_zip(key, archive=None, report=lambda _: None):
    folder = data_dir() / 'tools'
    folder.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=folder, prefix='download-') as stage:
        stage = Path(stage)
        if archive is None:
            if key != 'XTB_EXE':
                raise ValueError('请先下载该工具的官方 Windows ZIP。')
            archive = stage / 'xtb.zip'
            download(XTB_URL, archive, XTB_HASH, report)
        report('正在解压并检查程序…')
        extracted = stage / 'package'
        unpack(archive, extracted)
        exe = find_in_folder(extracted, key)
        if key == 'XTB_EXE':
            checked_run([exe, '--version'])
        relative = exe.relative_to(extracted)
        with tempfile.TemporaryDirectory(dir=folder, prefix='reserve-') as reservation:
            destination = folder / f'{exe.stem}-{Path(reservation).name[8:]}'
            extracted.rename(destination)
        return destination / relative


def install_molsimplify(report=lambda _: None, *, local_wheel=None):
    if local_wheel is None:
        existing = detect_tool('MOLSIMPLIFY_PYTHON')
        if existing:
            report('正在检查已配置的 molSimplify 环境…')
            try:
                checked_run([existing, '-I', '-c', PROBE_MOLSIMPLIFY], timeout=30)
            except (OSError, RuntimeError):
                report('已有环境未通过检查，正在准备新环境…')
            else:
                return existing
    if local_wheel is not None:
        with Path(local_wheel).open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != MOLSIMPLIFY_HASH:
                raise ValueError(f'请选择官方发布的 {MOLSIMPLIFY_WHEEL}。')
    root = data_dir() / 'tools' / 'molsimplify'
    root.mkdir(parents=True, exist_ok=True)
    # A venv cannot be moved after creation. Give each install a permanent path.
    folder = Path(tempfile.mkdtemp(dir=root, prefix='env-'))
    env = {**os.environ, 'UV_PYTHON_INSTALL_DIR': str(root / 'python'),
           'UV_CACHE_DIR': str(root / 'cache'), 'UV_NATIVE_TLS': 'true', 'PYTHONUTF8': '1'}
    for key in ('PYTHONHOME', 'PYTHONPATH', 'VIRTUAL_ENV'):
        env.pop(key, None)
    try:
        with tempfile.TemporaryDirectory(dir=root, prefix='uv-') as helper:
            helper = Path(helper)
            archive = PACKAGE_ROOT / 'uv-bootstrap.whl'
            if not archive.is_file():
                archive = helper / 'uv.whl'
                download(UV_URL, archive, UV_HASH, report)
            unpack(archive, helper / 'package')
            uv = next((helper / 'package').rglob('uv.exe'))
            report('正在查找本机可用的 Python 3.12…')
            found = tool_run([uv, 'python', 'find', '--no-python-downloads', '3.12'], timeout=30, env=env)
            if found.returncode:
                report('未找到 Python 3.12，正在联网下载独立运行环境…')
                checked_run([uv, 'python', 'install', '3.12'], timeout=600, env=env)
            report('正在创建 molSimplify 的独立 Python 3.12 环境…')
            checked_run([uv, 'venv', '--python', '3.12', '--no-python-downloads', folder], timeout=120, env=env)
            python = folder / 'Scripts' / 'python.exe'
            report('正在安装 molSimplify 和 Open Babel…')
            cmd = [uv, 'pip', 'install', '--python', python]
            checked_run([*cmd, '-r', PACKAGE_ROOT / 'requirements-molsimplify.txt'], timeout=720, env=env)
            wheel = helper / MOLSIMPLIFY_WHEEL
            if local_wheel is None:
                download(MOLSIMPLIFY_URL, wheel, MOLSIMPLIFY_HASH, report)
            else:
                shutil.copyfile(local_wheel, wheel)
            report('正在安装并验证 molSimplify…')
            checked_run([*cmd, '--no-deps', '--link-mode', 'copy', wheel], timeout=900, env=env)
            checked_run([python, '-c', PROBE_MOLSIMPLIFY], env=env)
        return python
    except Exception:
        shutil.rmtree(folder, ignore_errors=True)
        raise
