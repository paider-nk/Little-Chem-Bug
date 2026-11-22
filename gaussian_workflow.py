"""Run Gaussian jobs and create practical analyses from Gaussian outputs."""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import external_process as subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from analysis_core import (difference_density, energy_difference, parse_log,
                          solvation_difference, write_report)
from gaussian_core import InputError


def locate_gaussian(explicit=None):
    candidates = [explicit, os.environ.get("GAUSSIAN_EXE")]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    for name in ("g16.exe", "g09.exe", "g16"):
        found = shutil.which(name)
        if found:
            return found
    raise InputError("找不到 Gaussian。请用 --gaussian 指定 g16.exe，或设置 GAUSSIAN_EXE。")


def locate_formchk(explicit=None, gaussian=None):
    candidates = [explicit, os.environ.get("FORMCHK_EXE")]
    if gaussian:
        candidates.append(str(Path(gaussian).resolve().parent / "formchk.exe"))
    if not gaussian:
        try:
            candidates.append(str(Path(locate_gaussian()).with_name("formchk.exe")))
        except InputError:
            pass
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    found = shutil.which("formchk.exe") or shutil.which("formchk")
    if found:
        return found
    raise InputError("找不到 formchk.exe。请用 --formchk 指定路径或设置 FORMCHK_EXE。")


def locate_multiwfn(explicit=None):
    candidates = [explicit, os.environ.get("MULTIWFN_EXE")]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    found = shutil.which("Multiwfn.exe") or shutil.which("Multiwfn")
    if found:
        return found
    raise InputError("找不到 Multiwfn。请用 --multiwfn 指定 Multiwfn.exe 或设置 MULTIWFN_EXE。")


def _prepared_input(source: Path, target: Path, cores: int, memory: int):
    text = source.read_text(encoding="ascii", errors="strict")
    lines = text.splitlines()
    lines = [line for line in lines if not re.match(r"^\s*%(?:mem|nprocshared|chk)\s*=", line, re.I)]
    name = re.sub(r"[^A-Za-z0-9_.-]", "_", source.stem)[:60] or "job"
    lines = [f"%chk={name}.chk", f"%mem={memory}GB", f"%nprocshared={cores}", *lines]
    # Gaussian requires a blank record after the geometry (and after later
    # optional sections). Preserve it when inserting resource directives.
    target.write_text("\n".join(lines).rstrip() + "\n\n", encoding="ascii", newline="\n")


def _run_one(executable, source: Path, root: Path, cores, memory):
    folder = root / source.stem
    folder.mkdir(parents=True, exist_ok=True)
    gjf, log = folder / (source.stem + ".gjf"), folder / (source.stem + ".log")
    _prepared_input(source, gjf, cores, memory)
    env = os.environ.copy()
    env["GAUSS_EXEDIR"] = str(Path(executable).parent)
    env["GAUSS_SCRDIR"] = str(folder)
    try:
        with log.open("wb") as output:
            done = subprocess.run([executable, gjf.name, log.name], cwd=folder, env=env,
                                   stdout=output, stderr=subprocess.STDOUT, check=False)
    except OSError as exc:
        return {"input": str(source), "log": str(log), "status": "launch_error", "error": str(exc)}
    normal = False
    if log.exists():
        normal = "Normal termination" in log.read_text(encoding="utf-8", errors="replace")
    return {"input": str(source), "log": str(log), "return_code": done.returncode,
            "status": "normal" if normal else "incomplete", "normal_termination": normal}


def run_jobs(inputs, output, *, gaussian=None, orca=None, cores=4, memory=4, slots=1, on_progress=None):
    """Run a batch and optionally report (settled, total, abnormal) job counts."""
    paths = [Path(p).resolve() for p in inputs]
    if not paths or any(not p.is_file() or p.suffix.lower() not in (".gjf", ".com", ".inp") for p in paths):
        raise InputError("请提供至少一个存在的 .gjf / .com / .inp 输入文件。")
    if cores < 1 or memory < 1 or slots < 1:
        raise InputError("核心数、内存和并行槽位均须为正整数。")
    if cores * slots > (os.cpu_count() or cores):
        raise InputError(f"总并行核心数 {cores * slots} 超过本机逻辑处理器数 {os.cpu_count()}。")
    from orca_support import locate_orca, run_one as run_orca
    jobs = [(run_orca, locate_orca(orca), p) if p.suffix.lower() == ".inp" else
            (_run_one, locate_gaussian(gaussian), p) for p in paths]
    root = Path(output).resolve()
    root.mkdir(parents=True, exist_ok=True)
    results = []
    abnormal = 0
    if on_progress is not None:
        on_progress(0, len(paths), 0)
    with ThreadPoolExecutor(max_workers=min(slots, len(paths))) as pool:
        futures = [pool.submit(runner, exe, p, root, cores, memory) for runner, exe, p in jobs]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            abnormal += result["status"] != "normal"
            if on_progress is not None:
                on_progress(len(results), len(paths), abnormal)
    results.sort(key=lambda x: x["input"])
    (root / "run_report.json").write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    return results


def to_fchk(checkpoint, output=None, *, formchk=None, gaussian=None):
    source = Path(checkpoint).resolve()
    if not source.is_file() or source.suffix.lower() != ".chk":
        raise InputError("请提供已存在的 Gaussian .chk 检查点文件。")
    target = Path(output).resolve() if output else source.with_suffix(".fchk")
    exe = locate_formchk(formchk, gaussian)
    result = subprocess.run([exe, str(source), str(target)], capture_output=True, text=True, check=False)
    if result.returncode != 0 or not target.is_file():
        raise InputError("formchk 转换失败：" + (result.stderr or result.stdout or f"退出码 {result.returncode}"))
    return target


def launch_multiwfn(source, *, multiwfn=None, formchk=None, gaussian=None):
    path = Path(source).resolve()
    if not path.is_file():
        raise InputError(f"文件不存在：{path}")
    if path.suffix.lower() == ".gbw":
        from orca_support import gbw_to_molden
        path = gbw_to_molden(path)
    if path.suffix.lower() == ".chk":
        path = to_fchk(path, formchk=formchk, gaussian=gaussian)
    exe = locate_multiwfn(multiwfn)
    env = os.environ.copy()
    env["Multiwfnpath"] = str(Path(exe).parent)
    return subprocess.run([exe, str(path)], cwd=path.parent, env=env, check=False).returncode


def run_multiwfn_analysis(source, analysis, output=None, *, multiwfn=None, formchk=None, gaussian=None, on_stage=None):
    """Run a supported non-interactive Multiwfn analysis and save its transcript."""
    path = Path(source).resolve()
    if not path.is_file() or path.suffix.lower() not in (".chk", ".fchk", ".fch", ".gbw", ".molden", ".input"):
        raise InputError("请选择已存在的 Gaussian 检查点或 ORCA .gbw / Molden 文件。")
    stage = on_stage if on_stage is not None else lambda _: None
    if path.suffix.lower() == ".gbw":
        from orca_support import gbw_to_molden
        path = gbw_to_molden(path)
    if path.suffix.lower() == ".chk":
        stage("转换 .chk 为 .fchk")
        path = to_fchk(path, formchk=formchk, gaussian=gaussian)
    scripts = {
        "esp_surface": ("12\n0\n-1\n-1\nq\n", "Summary of surface analysis"),
        "esp_extrema": ("5\n12\n4\n0.2\n0\nq\n", "The minimum is"),
    }
    if analysis not in scripts:
        raise InputError("不支持的 Multiwfn 自动分析类型。")
    commands, marker = scripts[analysis]
    exe = locate_multiwfn(multiwfn)
    env = os.environ.copy()
    env["Multiwfnpath"] = str(Path(exe).parent)
    stage("运行 Multiwfn ESP 分析")
    result = subprocess.run([exe, str(path)], input=commands, text=True, encoding="utf-8",
                            errors="replace", cwd=path.parent, env=env,
                            capture_output=True, check=False, timeout=600)
    stage("保存 Multiwfn 完整输出")
    transcript = result.stdout + ("\n[stderr]\n" + result.stderr if result.stderr else "")
    target = Path(output).resolve() if output else path.with_name(path.stem + f"_{analysis}.txt")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(transcript, encoding="utf-8")
    if result.returncode != 0:
        raise InputError(f"Multiwfn 运行失败（退出码 {result.returncode}），完整输出：{target}")
    if marker not in transcript:
        raise InputError(f"Multiwfn 未生成预期的分析结果，请检查输入波函数及完整输出：{target}")
    summary_start = transcript.find(marker)
    summary_text = transcript[summary_start:]
    terminator = "Surface analysis finished!" if analysis == "esp_surface" else "Summing up all value"
    end = summary_text.find(terminator)
    if end >= 0:
        summary_text = summary_text[:end + len(terminator)]
    summary = summary_text.strip()
    return {"analysis": analysis, "wavefunction": str(path), "output": str(target),
            "summary": summary.strip()}


def cli(argv=None):
    parser = argparse.ArgumentParser(description="Gaussian 运行与 Multiwfn 辅助后处理工具")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="运行一个或多个 .gjf/.com 文件")
    run.add_argument("inputs", nargs="+")
    run.add_argument("-o", "--output", required=True, help="日志及检查点输出目录")
    run.add_argument("--gaussian", help="Gaussian 可执行文件；默认读取 GAUSSIAN_EXE 或 PATH")
    run.add_argument("--orca", help="ORCA 可执行文件或安装目录")
    run.add_argument("--cores", type=int, default=4)
    run.add_argument("--memory", type=int, default=4, help="每个作业使用的 GB")
    run.add_argument("--slots", type=int, default=1, help="并行作业数；总核心=cores×slots")

    analyze = sub.add_parser("analyze", help="分析 Gaussian 日志并生成 HTML/CSV/图")
    analyze.add_argument("log")
    analyze.add_argument("-o", "--output", help="报告目录，默认与日志相邻的 analysis 文件夹")

    energy = sub.add_parser("binding-energy", help="结合能：E(AB)-sum(E(fragment))")
    energy.add_argument("logs", nargs="+", help="依次为复合物日志、片段 A、片段 B……")
    energy.add_argument("--free-energy", action="store_true")
    energy.add_argument("-o", "--output", help="可选 JSON 输出文件")

    solv = sub.add_parser("solvation", help="溶剂化能/自由能差：solution - gas")
    solv.add_argument("gas_log")
    solv.add_argument("solution_log")
    solv.add_argument("--free-energy", action="store_true")
    solv.add_argument("-o", "--output")

    density = sub.add_parser("density-difference", help="同网格电子密度差：rho(AB)-rho(A)-rho(B)")
    density.add_argument("complex_cube")
    density.add_argument("fragment_a_cube")
    density.add_argument("fragment_b_cube")
    density.add_argument("-o", "--output", required=True)

    convert = sub.add_parser("formchk", help="Gaussian .chk 转换为 Multiwfn 可读的 .fchk")
    convert.add_argument("checkpoint")
    convert.add_argument("-o", "--output")
    convert.add_argument("--formchk", help="formchk.exe 路径")
    convert.add_argument("--gaussian", help="用于从 Gaussian 目录定位 formchk.exe")

    multi = sub.add_parser("multiwfn", help="启动 Multiwfn；输入 .chk 时自动先转换为 .fchk")
    multi.add_argument("wavefunction", help="Gaussian .fchk/.fch 或 .chk 检查点")
    multi.add_argument("--multiwfn", help="Multiwfn.exe 路径")
    multi.add_argument("--formchk", help="formchk.exe 路径")
    multi.add_argument("--gaussian", help="用于从 Gaussian 目录定位 formchk.exe")

    vmd = sub.add_parser("vmd", help="用 VMD 打开计算结果或结构文件")
    vmd.add_argument("source", help=".log/.out/.gbw/.xyz/.pdb/.cube/.molden 等文件")
    vmd.add_argument("--vmd", help=r"vmd.exe 路径；默认读取 VMD_EXE 或 PATH")
    vmd.add_argument("-o", "--output", help="日志提取的最终结构 .xyz 输出路径")

    esp = sub.add_parser("esp-map", help="生成 VMD 彩色 ESP 表面图并导出 PNG/TIFF")
    esp.add_argument("wavefunction", help="Gaussian .chk/.fchk 或 ORCA .gbw/.molden")
    esp.add_argument("-o", "--output", help="输出目录，默认在波函数旁新建目录")
    esp.add_argument("--iso", type=float, default=0.001, help="电子密度等值，默认 0.001 e/bohr^3")
    esp.add_argument("--limit", type=float, default=0.05, help="ESP 色标半宽，默认 0.05 a.u.")
    esp.add_argument("--quality", type=int, choices=(1, 2, 3), default=2, help="Multiwfn 网格精度：1 低、2 中、3 高")
    esp.add_argument("--no-legend", action="store_true", help="不在导出图片下方添加色标")
    esp.add_argument("--no-vmd", action="store_true", help="只导出图片，不打开 VMD 交互窗口")
    esp.add_argument("--vmd", help=r"vmd.exe 路径；默认读取 VMD_EXE 或 PATH")
    esp.add_argument("--orca", help=r"orca.exe 路径；默认读取 ORCA_EXE 或 PATH")

    mw_analyze = sub.add_parser("multiwfn-analyze", help="自动运行支持的 Multiwfn ESP 分析")
    mw_analyze.add_argument("wavefunction", help="Gaussian .chk/.fchk/.fch 波函数")
    mw_analyze.add_argument("analysis", choices=("esp_surface", "esp_extrema"))
    mw_analyze.add_argument("-o", "--output", help="完整终端输出保存路径")
    mw_analyze.add_argument("--multiwfn", help="Multiwfn.exe 路径")
    mw_analyze.add_argument("--formchk", help="formchk.exe 路径")
    mw_analyze.add_argument("--gaussian", help="用于从 Gaussian 目录定位 formchk.exe")

    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            result = run_jobs(args.inputs, args.output, gaussian=args.gaussian, orca=args.orca, cores=args.cores,
                              memory=args.memory, slots=args.slots)
            print(json.dumps(result, ensure_ascii=False, indent=2))
            return 0 if all(r["status"] == "normal" for r in result) else 2
        if args.command == "analyze":
            a = parse_log(args.log)
            report = write_report(a, args.output or (Path(args.log).resolve().parent / "analysis"))
            print(f"报告：{report}")
            return 0 if a.normal else 2
        if args.command in ("binding-energy", "solvation"):
            if args.command == "binding-energy":
                if len(args.logs) < 3:
                    raise InputError("结合能至少需要复合物和两个片段的日志。")
                value = energy_difference(parse_log(args.logs[0]), [parse_log(p) for p in args.logs[1:]], free_energy=args.free_energy)
            else:
                value = solvation_difference(parse_log(args.gas_log), parse_log(args.solution_log), free_energy=args.free_energy)
            rendered = json.dumps(value, ensure_ascii=False, indent=2)
            if args.output:
                Path(args.output).write_text(rendered, encoding="utf-8")
            print(rendered)
            return 0
        if args.command == "density-difference":
            value = difference_density(args.complex_cube, args.fragment_a_cube, args.fragment_b_cube, args.output)
            print(json.dumps(value, ensure_ascii=False, indent=2))
            return 0
        if args.command == "formchk":
            target = to_fchk(args.checkpoint, args.output, formchk=args.formchk, gaussian=args.gaussian)
            print(f"Multiwfn 波函数文件：{target}")
            return 0
        if args.command == "vmd":
            from vmd_support import open_vmd
            print(json.dumps(open_vmd(args.source, vmd=args.vmd, output=args.output), ensure_ascii=False, indent=2))
            return 0
        if args.command == "multiwfn":
            return launch_multiwfn(args.wavefunction, multiwfn=args.multiwfn,
                                   formchk=args.formchk, gaussian=args.gaussian)
        if args.command == "esp-map":
            from esp_map import generate_esp
            value = generate_esp(args.wavefunction, args.output, iso=args.iso, limit=args.limit,
                                 quality=args.quality, legend=not args.no_legend,
                                 interactive=not args.no_vmd, vmd=args.vmd, orca=args.orca,
                                 on_stage=lambda message: print(message, flush=True))
            print(json.dumps(value, ensure_ascii=False, indent=2))
            return 0
        if args.command == "multiwfn-analyze":
            print(json.dumps(run_multiwfn_analysis(args.wavefunction, args.analysis, args.output,
                          multiwfn=args.multiwfn, formchk=args.formchk, gaussian=args.gaussian),
                          ensure_ascii=False, indent=2))
            return 0
    except (InputError, OSError, ValueError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(cli())
