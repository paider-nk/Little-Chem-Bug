"""Build a Windows x64 portable release from a development environment."""
from pathlib import Path
import argparse
import hashlib
import os
import shutil
import struct
import subprocess
import sys
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parent
VERSION = '1.0.1'


def build(output, build_root, uv_wheel=None):
    if os.name != 'nt' or struct.calcsize('P') != 8:
        raise RuntimeError('Build with a Windows x64 Python environment.')
    from tool_setup import UV_URL, UV_HASH, download
    from workflow_gui import configure_tk_paths
    configure_tk_paths()
    build_root.mkdir(parents=True, exist_ok=True)
    bootstrap = build_root / 'uv-bootstrap.whl'
    if uv_wheel is not None:
        shutil.copyfile(uv_wheel, bootstrap)
    elif not bootstrap.is_file():
        download(UV_URL, bootstrap, UV_HASH, print)
    with bootstrap.open('rb') as stream:
        if hashlib.file_digest(stream, 'sha256').hexdigest() != UV_HASH:
            raise RuntimeError('Unexpected hash for the Windows uv wheel.')
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    cmd = [sys.executable, '-m', 'PyInstaller', '--clean', '--noconfirm', '--windowed', '--onedir',
           '--name', 'ChemWorkbench', '--distpath', str(build_root / 'dist'),
           '--workpath', str(build_root / 'work'), '--specpath', str(build_root / 'spec'),
           '--collect-data', 'rdkit']
    for dependency in ('rdkit', 'numpy', 'scipy', 'ase', 'gemmi', 'matplotlib', 'Pillow', 'pyinstaller'):
        cmd.extend(('--copy-metadata', dependency))
    for filename in ('assets', 'molsimplify_worker.py', 'requirements-molsimplify.txt'):
        dest = filename if filename == 'assets' else '.'
        cmd.extend(('--add-data', f'{ROOT / filename}:{dest}'))
    cmd.extend(('--add-data', f'{bootstrap}:.'))
    subprocess.run([*cmd, str(ROOT / 'chem_workbench.py')], cwd=ROOT, env=env, check=True)
    folder = build_root / 'dist' / 'ChemWorkbench'
    internal = folder / '_internal'
    if not (internal / '_tkinter.pyd').is_file():
        raise RuntimeError('The current Python installation did not provide a working Tcl/Tk.')
    shutil.copy2(ROOT / '使用说明.md', folder / '使用说明.md')
    shutil.copytree(ROOT / 'examples', folder / 'examples')
    license_file = Path(sys.base_prefix) / 'LICENSE.txt'
    if license_file.is_file():
        shutil.copy2(license_file, internal / 'PYTHON-LICENSE.txt')
    output.mkdir(parents=True, exist_ok=True)
    archive = output / f'ChemWorkbench-{VERSION}-Windows-x64.zip'
    partial = archive.with_suffix('.zip.tmp')
    with ZipFile(partial, 'w', compression=ZIP_DEFLATED, compresslevel=5) as zf:
        for file in sorted(folder.rglob('*')):
            if file.is_file():
                zf.write(file, file.relative_to(folder.parent).as_posix())
    partial.replace(archive)
    print(archive)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT.parent)
    parser.add_argument('--build-root', type=Path, default=ROOT / 'build-release')
    parser.add_argument('--uv-wheel', type=Path, help='Path to the official uv 0.9.28 Windows x64 wheel')
    args = parser.parse_args()
    build(args.output.resolve(), args.build_root.resolve(), args.uv_wheel)
