"""Map wavefunction ESP onto an electron density surface and export VMD's render."""
from pathlib import Path
import math
import os
import shutil
import external_process as subprocess
import tempfile
import time

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from analysis_core import read_cube
from gaussian_core import InputError
from gaussian_workflow import locate_multiwfn, to_fchk
from orca_support import gbw_to_molden
from vmd_support import locate_vmd, vmd_environment


def export_field(executable, wave, folder, field, grid, name):
    log = folder / (name + '.log')
    with log.open('wb') as stream:
        result = subprocess.run([executable, wave.name], cwd=folder,
                                env=dict(os.environ, Multiwfnpath=str(Path(executable).parent)),
                                input=f'5\n{field}\n{grid}\n2\n0\nq\n'.encode('ascii'),
                                stdout=stream, stderr=subprocess.STDOUT,
                                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
    candidates = list(folder.glob('*' + name + '.cub'))
    if result.returncode or len(candidates) != 1:
        raise InputError(f'Multiwfn 未能导出 {name}，详见 {log}')
    candidates[0].replace(folder / (name + '.cub'))


def scene_parameters(folder, iso):
    density, esp = (read_cube(folder / name) for name in ('density.cub', 'totesp.cub'))
    if (density.values.shape != esp.values.shape or density.atoms.shape != esp.atoms.shape or
            any(not np.allclose(a, b, atol=1e-6, rtol=0) for a, b in
                ((density.origin, esp.origin), (density.axes, esp.axes), (density.atoms, esp.atoms)))):
        raise InputError('电子密度与 ESP 的网格或原子坐标不一致。')
    rho = density.values
    if not rho.min() < iso < rho.max():
        raise InputError('此密度等值面不存在。请调整密度等值或重新计算波函数。')
    for axis in range(3):
        if np.take(rho, [0, -1], axis=axis).max() >= iso:
            raise InputError('等值面触及网格边界，请提高密度等值或增大 Multiwfn 网格外延。')
    xyz = density.atoms[:, 2:5]
    # View along the shortest principal axis so planar molecules appear face-on.
    rotation = np.eye(4)
    if len(xyz) > 1:
        _, axes = np.linalg.eigh(np.cov(xyz.T))
        rotation[:2, :3] = axes[:, ::-1].T[:2]
        rotation[2, :3] = np.cross(rotation[0, :3], rotation[1, :3])
    return ' '.join('{' + ' '.join(f'{v:.9g}' for v in row) + '}' for row in rotation)


def run_scene(executable, folder, render):
    env = vmd_environment(executable)
    status = folder / ('render.status' if render else 'viewer.status')
    env['ESP_STATUS'] = status.name
    args = [executable]
    if render:
        env['ESP_RENDER'] = '1'
        args += ['-dispdev', 'text', '-size', '1600', '1000']
    else:
        env.pop('ESP_RENDER', None)
    with (folder / (status.stem + '.log')).open('wb') as log:
        proc = subprocess.Popen([*args, '-e', 'scene.tcl'], cwd=folder, env=env,
                                stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL if render else None,
                                creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0) if render else 0)
        deadline = time.monotonic() + (900 if render else 120)
        while not status.exists() and proc.poll() is None:
            if time.monotonic() >= deadline:
                proc.kill(); proc.wait()
                break
            time.sleep(0.1)
    result = status.read_text(errors='replace').strip() if status.exists() else 'VMD 未完成加载或渲染'
    if not result.startswith('0 '):
        raise InputError(f'{result}。详见 {folder / (status.stem + ".log")}')
    if render:
        proc.wait(timeout=30)


def save_images(folder, limit, iso, legend):
    with Image.open(folder / 'surface.tga') as raw:
        image = raw.convert('RGB')
    if legend:
        palette = np.loadtxt(folder / 'palette.txt')
        width = image.width // 3
        indices = np.linspace(0, len(palette)-1, width).astype(int)
        colors = np.clip(palette[indices]*255, 0, 255).astype('uint8')
        bar = Image.fromarray(colors[None, :, :]).resize((width, 15))
        footer = Image.new('RGB', (image.width, image.height + 94), 'white')
        footer.paste(image)
        x = (image.width - width)//2
        footer.paste(bar, (x, image.height + 12))
        draw = ImageDraw.Draw(footer)
        from matplotlib import font_manager
        font = ImageFont.truetype(font_manager.findfont('DejaVu Sans'), 18)
        for dx, text in [(0, f'-{limit:g}'), (width//2, '0'), (width, f'+{limit:g}')]:
            draw.text((x+dx, image.height+31), text, font=font, fill='#232535', anchor='mt')
        draw.text((image.width//2, image.height+60), f'ESP (a.u.)     |     Density = {iso:g} e/bohr³',
                  font=font, fill='#232535', anchor='mt')
        image = footer
    image.save(folder / 'ESP.png', dpi=(300, 300))
    image.save(folder / 'ESP.tif', compression='tiff_lzw', dpi=(300, 300))


def generate_esp(source, output=None, *, iso=0.001, limit=0.05, quality=2, legend=True,
                 interactive=True, vmd=None, gaussian=None, orca=None, on_stage=lambda _: None):
    source = Path(source).resolve()
    suffix = '.molden.input' if source.name.lower().endswith('.molden.input') else source.suffix.lower()
    if not source.is_file() or suffix not in ('.chk', '.fchk', '.fch', '.gbw', '.molden', '.molden.input'):
        raise InputError('请选择 Gaussian 检查点或 ORCA GBW / Molden 波函数文件。')
    if not all(math.isfinite(x) and x > 0 for x in (iso, limit)) or quality not in (1, 2, 3):
        raise InputError('密度等值和色标半宽必须为有限正数，请选择有效的网格精度。')
    executable, multi = locate_vmd(vmd), locate_multiwfn()
    parent = Path(output).resolve() if output else source.parent
    parent.mkdir(parents=True, exist_ok=True)
    folder = Path(tempfile.mkdtemp(prefix=source.stem + '_ESP_', dir=parent))
    on_stage('准备波函数')
    # ASCII local filenames also work with legacy Fortran readers on Windows.
    wave = folder / ('wave' + suffix)
    shutil.copyfile(source, wave)
    if suffix == '.chk':
        wave = to_fchk(wave, gaussian=gaussian)
    elif suffix == '.gbw':
        wave = gbw_to_molden(wave, orca=orca)
    on_stage('计算电子密度网格')
    export_field(multi, wave, folder, 1, str(quality), 'density')
    on_stage('在同一网格计算总静电势')
    export_field(multi, wave, folder, 12, '8\ndensity.cub', 'totesp')
    on_stage('校验两份网格并用 VMD 渲染')
    rotation = scene_parameters(folder, iso)
    template = (Path(__file__).parent / 'assets/esp_surface.tcl').read_text('ascii')
    scene = template.replace('@ISO@', str(iso)).replace('@LIMIT@', str(limit)).replace('@ROTATION@', rotation)
    (folder / 'scene.tcl').write_text(scene, encoding='ascii')
    run_scene(executable, folder, True)
    save_images(folder, limit, iso, legend)
    result = {'folder': str(folder), 'png': str(folder / 'ESP.png'), 'tiff': str(folder / 'ESP.tif')}
    if interactive:
        on_stage('打开可旋转的 VMD 表面')
        try:
            run_scene(executable, folder, False)
        except (InputError, OSError) as exc:
            result['viewer_warning'] = str(exc)
    return result

