"""Desktop interface and command-line entry point."""
from __future__ import annotations

import argparse
import json
import os
import queue
import sys
import threading
from dataclasses import asdict
from pathlib import Path

from gaussian_core import (InputError, Settings, default_checkpoint, generate,
                           PREPROCESS_METHODS, PREPROCESS_SCENARIOS)


def optional_integer(text: str, label: str):
    if text.strip().casefold() in ("", "auto", "自动"):
        return None
    try:
        return int(text)
    except ValueError as exc:
        raise InputError(f"{label}请输入整数，或保留“自动”。") from exc


def run_gui(parent=None, on_send=None):
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    from mascot_ui import MascotDock

    # A moved Windows Python installation can inherit stale global Tcl/Tk paths.
    # Repair only this process; do not change the machine's environment variables.
    for env_name, version, filename in (("TCL_LIBRARY", tk.TclVersion, "init.tcl"),
                                        ("TK_LIBRARY", tk.TkVersion, "tk.tcl")):
        existing = os.environ.get(env_name)
        if existing and (Path(existing) / filename).is_file():
            continue
        folder = ("tcl" if env_name == "TCL_LIBRARY" else "tk") + str(version)
        for parent in (Path(sys.base_prefix) / "tcl", Path(sys.base_prefix) / "Library" / "lib"):
            candidate = parent / folder
            if (candidate / filename).is_file():
                os.environ[env_name] = str(candidate)
                break

    root = parent if parent is not None else tk.Tk()
    if parent is None:
        root.title("Gaussian Input · 分子输入文件生成器")
        root.geometry("1060x850")
        root.minsize(850, 720)
        root.configure(bg="#f3f6fb")
    style = ttk.Style()
    style.theme_use("clam")
    style.configure("TFrame", background="#f3f6fb")
    style.configure("TLabel", background="#f3f6fb", font=("Microsoft YaHei UI", 10))
    style.configure("TButton", font=("Microsoft YaHei UI", 10), padding=(12, 7))
    style.configure("TEntry", padding=5)
    style.configure("TCombobox", padding=5)
    style.configure("Title.TLabel", font=("Microsoft YaHei UI", 21, "bold"), foreground="#162e51")
    style.configure("Muted.TLabel", foreground="#52647c")
    style.configure("Accent.TButton", foreground="white", background="#245ddd")
    style.map("Accent.TButton", background=[("active", "#164bb8"), ("disabled", "#8da5d6")])

    shell = ttk.Frame(root, padding=24)
    shell.pack(fill="both", expand=True)
    mascot = MascotDock(shell, "输入 SMILES 或分子名称，再点击「生成输入文件」。")
    mascot.pack(side="bottom", fill="x", pady=(10, 0))
    if parent is not None:
        mascot._toggle()
    main = ttk.Frame(shell)
    main.pack(fill="both", expand=True)
    main.columnconfigure(0, weight=1)
    main.rowconfigure(6, weight=1)
    ttk.Label(main, text="分子与晶体 → Gaussian / ORCA", style="Title.TLabel").grid(row=0, column=0, sticky="w")
    ttk.Label(main, text="SMILES / 名称 / 聚合物 / 配合物 / 晶体 CIF", style="Muted.TLabel").grid(
        row=1, column=0, sticky="w", pady=(4, 17))

    variables = {key: tk.StringVar(value=value) for key, value in {
        "input": "CCO", "mode": "自动识别", "method": "B3LYP", "basis": "6-31G(d)",
        "job": "优化 + 频率", "charge": "自动", "multiplicity": "自动", "memory": "4",
        "processors": "4", "conformers": "20", "extra": "SCF=XQC Int=UltraFine",
        "title": "Molecule", "checkpoint": "auto", "engine": "Gaussian",
        "cif_scope": "molecule", "cif_block": "", "cif_component": "1",
        "degree": "5", "end_group": "氢封端", "complex_spec": "",
        "boron_preprocess": "严格保键", "boron_bond_tolerance": "0.15",
        "preprocess_method": "自动（MMFF94s → UFF）",
    }.items()}

    input_frame = ttk.Frame(main)
    input_frame.grid(row=2, column=0, sticky="ew")
    input_frame.columnconfigure(0, weight=1)
    ttk.Label(input_frame, text="RDKit 标准 SMILES / 分子名称 / 链节 SMILES（两个 * 为连接点）").grid(
        row=0, column=0, sticky="w", pady=(0, 5))
    input_entry = ttk.Entry(input_frame, textvariable=variables["input"], font=("Consolas", 12))
    input_entry.grid(row=1, column=0, sticky="ew", padx=(0, 12))
    ttk.Combobox(input_frame, textvariable=variables["mode"], state="readonly", width=12,
                 values=("自动识别", "SMILES", "分子名称", "聚合物链节", "二茂铁模板", "配合物 molSimplify", "晶体 CIF")).grid(row=1, column=1)
    examples = ttk.Frame(input_frame)
    examples.grid(row=2, column=0, columnspan=2, sticky="w", pady=(7, 10))
    ttk.Label(examples, text="示例：", style="Muted.TLabel").pack(side="left")

    def use_example(value, polymer=False):
        variables["input"].set(value)
        variables["mode"].set("聚合物链节" if polymer else "自动识别")
        if polymer:
            variables["degree"].set("5")
            variables["end_group"].set("氢封端")

    def use_ferrocene_template():
        variables["input"].set("ferrocene")
        variables["mode"].set("二茂铁模板")
        variables["charge"].set("自动")
        variables["multiplicity"].set("自动")

    for label, value in (("乙醇", "CCO"), ("苯", "苯"), ("乙酸根", "CC(=O)[O-]"), ("咖啡因（联网）", "caffeine")):
        ttk.Button(examples, text=label, command=lambda v=value: use_example(v)).pack(side="left", padx=3)
    ttk.Button(examples, text="聚乙烯", command=lambda: use_example("*CC*", True)).pack(side="left", padx=3)
    ttk.Button(examples, text="聚醚链", command=lambda: use_example("[*:1]CCO[*:2]", True)).pack(side="left", padx=3)
    ttk.Button(examples, text="二茂铁（夹心模板）", command=use_ferrocene_template).pack(side="left", padx=3)

    def pick_cif():
        from crystal_input import choose_crystal
        choose_crystal(root, variables)
    ttk.Button(examples, text="导入 CIF…", command=pick_cif).pack(side="left", padx=3)
    polymer_options = ttk.Frame(input_frame)
    polymer_options.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(0, 5))
    ttk.Label(polymer_options, text="聚合度 n").pack(side="left", padx=(0, 8))
    degree_entry = ttk.Spinbox(polymer_options, textvariable=variables["degree"], from_=1, to=100, width=6,
                               state="disabled")
    degree_entry.pack(side="left", padx=(0, 16))
    ttk.Label(polymer_options, text="两端封端").pack(side="left", padx=(0, 8))
    cap_combo = ttk.Combobox(polymer_options, textvariable=variables["end_group"],
                            values=("氢封端", "甲基封端"), width=10, state="disabled")
    cap_combo.pack(side="left", padx=(0, 16))
    ttk.Label(polymer_options, text="尾 → 头单键连接；最终 ≤200 个补氢前原子", style="Muted.TLabel").pack(side="left")
    def configure_complex():
        from complex_gui import open_complex_editor
        open_complex_editor(root, variables)
    ttk.Button(polymer_options, text="配置配合物（molSimplify）", command=configure_complex).pack(side="right")

    def update_polymer_controls(*_):
        enabled = variables["mode"].get() == "聚合物链节"
        degree_entry.configure(state="normal" if enabled else "disabled")
        cap_combo.configure(state="readonly" if enabled else "disabled")
        input_entry.configure(state="readonly" if variables["mode"].get() == "配合物 molSimplify" else "normal")

    variables["mode"].trace_add("write", update_polymer_controls)

    options = ttk.Frame(main)
    options.grid(row=3, column=0, sticky="ew", pady=(2, 12))
    for col in range(8):
        options.columnconfigure(col, weight=1 if col % 2 else 0)

    def field(label, key, row, col, choices=None, width=13):
        ttk.Label(options, text=label).grid(row=row, column=col * 2, sticky="w", padx=(0 if col == 0 else 12, 6), pady=6)
        if choices:
            widget = ttk.Combobox(options, textvariable=variables[key], values=choices, width=width,
                                  state="readonly" if key in ("job", "boron_preprocess", "preprocess_method") else "normal")
        else:
            widget = ttk.Entry(options, textvariable=variables[key], width=width)
        widget.grid(row=row, column=col * 2 + 1, sticky="ew", pady=6)

    field("方法", "method", 0, 0, ("B3LYP", "PBE1PBE", "M062X", "wB97XD", "HF"))
    field("基组", "basis", 0, 1, ("6-31G(d)", "6-31+G(d,p)", "6-311++G(d,p)", "Def2SVP", "Def2TZVP"))
    field("任务", "job", 0, 2, ("优化 + 频率", "几何优化", "单点能", "频率"))
    field("构象数", "conformers", 0, 3, width=7)
    field("电荷", "charge", 1, 0)
    field("多重度", "multiplicity", 1, 1)
    field("内存 / GB", "memory", 1, 2)
    field("CPU 核数", "processors", 1, 3, width=7)
    field("标题（英文）", "title", 2, 0)
    field("检查点名称", "checkpoint", 2, 1)
    ttk.Label(options, text="附加关键词").grid(row=2, column=4, sticky="w", padx=(12, 6))
    ttk.Entry(options, textvariable=variables["extra"]).grid(row=2, column=5, columnspan=3, sticky="ew")
    boron_modes = {"严格保键": "strict", "放宽保键": "relaxed",
                   "几何投影（保键共面）": "projected", "跳过平面化预处理": "skip"}
    field("硼预处理", "boron_preprocess", 3, 0, tuple(boron_modes), width=17)
    field("放宽容差 / Å", "boron_bond_tolerance", 3, 1)
    preprocess_choices = tuple(PREPROCESS_METHODS.values())
    preprocess_codes = {label: code for code, label in PREPROCESS_METHODS.items()}
    field("整体预处理", "preprocess_method", 4, 0, preprocess_choices, width=25)
    def update_engine(*_):
        current = variables['extra'].get()
        if current in ('TightSCF', 'SCF=XQC Int=UltraFine'):
            variables['extra'].set('TightSCF' if variables['engine'].get() == 'ORCA' else 'SCF=XQC Int=UltraFine')
    variables['engine'].trace_add('write', update_engine)
    preprocess_help = tk.StringVar()
    ttk.Label(options, textvariable=preprocess_help, style="Muted.TLabel", wraplength=650,
              justify="left").grid(row=4, column=2, columnspan=6, sticky="w", padx=12)

    def update_preprocess_help(*_):
        preprocess_help.set(PREPROCESS_SCENARIOS.get(preprocess_codes.get(variables["preprocess_method"].get()), ""))
    variables["preprocess_method"].trace_add("write", update_preprocess_help)
    update_preprocess_help()

    toolbar = ttk.Frame(main)
    toolbar.grid(row=4, column=0, sticky="ew", pady=(0, 10))

    info = tk.Text(main, height=4, wrap="word", font=("Microsoft YaHei UI", 10), relief="flat",
                   bg="#e8eef9", fg="#244064", padx=12, pady=9, state="disabled")
    info.grid(row=5, column=0, sticky="ew", pady=(0, 10))
    notebook = ttk.Notebook(main)
    notebook.grid(row=6, column=0, sticky="nsew")
    preview_frame = ttk.Frame(notebook)
    notes_frame = ttk.Frame(notebook)
    notebook.add(preview_frame, text=" 计算输入预览 ")
    notebook.add(notes_frame, text=" 结构与生成说明 ")

    def text_area(parent, font):
        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)
        widget = tk.Text(parent, wrap="none", font=font, bg="#ffffff", fg="#172b49", width=35, height=6,
                         relief="flat", padx=12, pady=10, undo=True)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=widget.yview)
        xscroll = ttk.Scrollbar(parent, orient="horizontal", command=widget.xview)
        widget.configure(yscrollcommand=scrollbar.set, xscrollcommand=xscroll.set)
        widget.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        xscroll.grid(row=1, column=0, sticky="ew")
        return widget

    preview = text_area(preview_frame, ("Consolas", 11))
    preview.configure(state="disabled")
    notes = text_area(notes_frame, ("Microsoft YaHei UI", 10))
    notes.configure(wrap="word", state="disabled")

    messages = queue.Queue()
    state = {"result": None, "settings": None, "input": None, "mode": None, "snapshot": None, "busy": False}

    def set_text(widget, content):
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", content)
        widget.configure(state="disabled")

    def refresh_result():
        # Never permit exporting a result that belongs to older form values.
        current = {k: v.get() for k, v in variables.items()}
        valid = state["result"] is not None and state["snapshot"] == current and not state["busy"]
        save_button.configure(state="normal" if valid else "disabled")
        copy_button.configure(state="normal" if valid else "disabled")
        if send_button is not None:
            send_button.configure(state="normal" if valid else "disabled")
        if state["busy"] or state["result"] is None:
            return
        if not valid and mascot.state in ("idle", "done"):
            mascot.set_state("idle", "参数已修改", "重新生成后即可保存当前参数对应的输入。")
        elif valid and mascot.state == "idle":
            result = state["result"]
            mascot.set_state("done", "输入文件已生成", f"{result.formula} · 可以保存当前计算输入。")

    def read_settings():
        try:
            is_polymer = variables["mode"].get() == "聚合物链节"
            try:
                degree = int(variables["degree"].get()) if is_polymer else 5
            except ValueError as exc:
                raise InputError("聚合度请输入 1–100 的整数。") from exc
            end_group = "methyl" if variables["end_group"].get() == "甲基封端" else "hydrogen"
            return Settings(
                method=variables["method"].get().strip(), basis=variables["basis"].get().strip(),
                job={"优化 + 频率": "optfreq", "几何优化": "opt", "单点能": "sp", "频率": "freq"}[variables["job"].get()],
                memory_gb=int(variables["memory"].get()), processors=int(variables["processors"].get()),
                charge=optional_integer(variables["charge"].get(), "电荷"),
                multiplicity=optional_integer(variables["multiplicity"].get(), "多重度"),
                conformers=int(variables["conformers"].get()), extra_keywords=variables["extra"].get().strip(),
                title=variables["title"].get().strip(),
                checkpoint=(default_checkpoint(variables["input"].get(), "polymer" if is_polymer else "auto",
                                               degree, end_group)
                            if variables["checkpoint"].get().strip().casefold() in ("", "auto", "自动")
                            else variables["checkpoint"].get().strip()),
                degree=degree, end_group=end_group, engine=variables["engine"].get().lower(),
                cif_block=variables["cif_block"].get(), cif_component=int(variables["cif_component"].get()), cif_scope=variables["cif_scope"].get(),
                boron_preprocess=boron_modes[variables["boron_preprocess"].get()],
                boron_bond_tolerance=float(variables["boron_bond_tolerance"].get()),
                preprocess_method=preprocess_codes[variables["preprocess_method"].get()],
            )
        except ValueError as exc:
            if isinstance(exc, InputError):
                raise
            raise InputError("内存、CPU 核数和构象数请输入整数；放宽容差请输入数值。") from exc

    def on_generate():
        if state["busy"]:
            return
        try:
            settings = read_settings()
        except InputError as exc:
            mascot.set_state("error", "输入参数有误", str(exc))
            messagebox.showerror("参数有误", str(exc))
            return
        value = variables["input"].get()
        mode = {"自动识别": "auto", "SMILES": "smiles", "分子名称": "name", "聚合物链节": "polymer",
                "二茂铁模板": "ferrocene", "配合物 molSimplify": "complex", "晶体 CIF": "cif"}[variables["mode"].get()]
        if mode == "complex":
            value = variables["complex_spec"].get()
            if not value:
                mascot.set_state("error", "缺少配合物配方", "请先配置金属、配体和配位构型。")
                messagebox.showerror("缺少配方", "请先点击“配置配合物（molSimplify）”。")
                return
        snapshot = {k: v.get() for k, v in variables.items()}
        state["busy"] = True
        generate_button.configure(state="disabled")
        refresh_result()
        mascot.set_state("running", "开始生成输入文件", "构象与预处理耗时因分子而异；进度条不代表百分比。")

        def worker():
            try:
                result = generate(value, mode, settings,
                                  on_stage=lambda stage: messages.put(("stage", stage)))
                messages.put(("ok", (result, settings, value, mode, snapshot)))
            except Exception as exc:
                messages.put(("error", str(exc) if isinstance(exc, InputError) else f"生成失败：{exc}"))

        threading.Thread(target=worker, daemon=True).start()

    def poll():
        try:
            while True:
                kind, payload = messages.get_nowait()
                if kind != "stage":
                    break
                mascot.set_state("running", payload, "处理时间取决于分子大小和所选方法。")
        except queue.Empty:
            root.after(100, poll)
            return
        state["busy"] = False
        generate_button.configure(state="normal")
        if kind == "error":
            state["result"] = None
            set_text(preview, "")
            set_text(info, "生成失败，请按提示修改输入或参数。")
            set_text(notes, "")
            mascot.set_state("error", "未能生成输入文件", "请查看错误详情，检查输入结构与计算参数后重试。")
            messagebox.showerror("无法生成", payload)
        else:
            result, settings, value, mode, snapshot = payload
            state.update(result=result, settings=settings, input=value, mode=mode, snapshot=snapshot)
            set_text(preview, result.text)
            set_text(info, f"分子式：{result.formula}     原子数（含氢）：{result.atom_count}     "
                          f"电荷：{result.charge:+d}     自旋多重度：{result.multiplicity}\n"
                          f"SMILES：{result.smiles or '晶体坐标输入'}\n来源：{result.source}   |   {result.geometry_method}")
            energy = f"\n力场相对筛选用能量：{result.energy:.6f} kcal/mol\n" if result.energy is not None else "\n"
            polymer_info = (f"链节 SMILES：{result.repeat_unit}\n聚合度：{result.degree}\n"
                            f"封端：{'氢' if result.end_group == 'hydrogen' else '甲基'}\n\n"
                            if result.degree is not None else "")
            set_text(notes, polymer_info + f"完整分子链 SMILES\n{result.smiles}\n\n三维结构生成\n{result.geometry_method}\n" + energy +
                     "\n" + "\n\n".join(f"• {note}" for note in result.warnings) +
                     "\n\n默认计算：B3LYP/6-31G(d)，Opt Freq。方法与基组是否适合目标体系需自行选择。"
                     "\n保存输入后，可在工作流界面提交 Gaussian / ORCA 计算；最终优化结果应检查收敛与虚频。")
            mascot.set_state("done", "输入文件已生成",
                             f"{result.formula} · {result.atom_count} 个原子；请核对预览后保存。")
            notebook.select(preview_frame)
        refresh_result()
        root.after(100, poll)

    def on_save(send=False):
        refresh_result()
        if str(save_button["state"]) == "disabled":
            return
        is_orca = state['settings'].engine == 'orca'
        ext = '.inp' if is_orca else '.gjf'
        path = filedialog.asksaveasfilename(
            parent=root, title="保存计算输入文件", initialfile=state["settings"].checkpoint + ext,
            defaultextension=ext, filetypes=[("ORCA input", "*.inp")] if is_orca else [("Gaussian input", "*.gjf"), ("Gaussian input", "*.com")])
        if not path:
            return
        try:
            Path(path).write_text(state["result"].text, encoding="ascii", newline="\n")
        except OSError as exc:
            mascot.set_state("error", "无法写入文件", "请检查目标目录和写入权限，再试一次。")
            messagebox.showerror("保存失败", str(exc))
            return
        mascot.set_state("done", "输入文件已保存", str(path))
        if send and on_send is not None:
            on_send(path)

    def on_copy():
        refresh_result()
        if str(copy_button["state"]) == "disabled":
            return
        root.clipboard_clear()
        root.clipboard_append(state["result"].text)
        mascot.set_state("done", "已复制输入内容", "可以粘贴到文本编辑器或其他工具。")

    ttk.Combobox(toolbar, textvariable=variables["engine"], values=("Gaussian", "ORCA"), state="readonly", width=9).grid(row=0,column=3,padx=8)
    generate_button = ttk.Button(toolbar, text="生成输入文件", style="Accent.TButton", command=on_generate)
    generate_button.grid(row=0, column=0, padx=(0, 8))
    save_button = ttk.Button(toolbar, text="保存输入文件", command=on_save, state="disabled")
    save_button.grid(row=0, column=1, padx=(0, 8))
    copy_button = ttk.Button(toolbar, text="复制内容", command=on_copy, state="disabled")
    copy_button.grid(row=0, column=2)
    send_button = None
    if on_send is not None:
        send_button = ttk.Button(toolbar, text="保存并送入计算", command=lambda: on_save(True), state="disabled")
        send_button.grid(row=0, column=4, padx=(6, 0))
    set_text(info, "输入 RDKit 支持的标准 SMILES 可离线生成，包含合法的点号分隔多片段写法；独立片段仅分开放置，不推断其结合/配位关系。内置常见中文名称；其他名称通过 PubChem 联网查询。\n"
             "聚合物：选择“聚合物链节”，例如 *CC*，再填聚合度与封端。电荷 / 多重度指整条分子链。")
    for variable in variables.values():
        variable.trace_add("write", lambda *_: refresh_result())
    root.bind("<Control-Return>", lambda _: on_generate() if not state["busy"] else None)
    root.after(100, poll)
    if parent is None:
        root.mainloop()


def cli(argv=None):
    parser = argparse.ArgumentParser(description="根据 SMILES、名称或聚合物链节生成 Gaussian .gjf / .com 输入。无参数启动图形界面。")
    parser.add_argument("molecule", nargs="?", help="SMILES、分子名称或含两个 * 的链节 SMILES")
    parser.add_argument("--gui", action="store_true", help="启动图形界面")
    parser.add_argument("--mode", choices=["auto", "smiles", "name", "polymer", "ferrocene", "complex", "cif"], default="auto",
                        help="complex 模式：molecule 为配合物配置 JSON 文件，须指定总电荷及多重度")
    parser.add_argument("--engine", choices=("gaussian","orca"), help="默认按输出后缀选择，未指定时为 Gaussian")
    parser.add_argument("--cif-block", default="", help="CIF 数据块名称")
    parser.add_argument("--cif-component", type=int, default=1, help="展开晶胞后的片段编号，从 1 开始")
    parser.add_argument("--cif-scope", choices=("molecule","periodic"), default="molecule")
    parser.add_argument("--degree", type=int, help="聚合度 1–100；仅 --mode polymer，默认 5")
    parser.add_argument("--end-group", choices=["hydrogen", "methyl"], help="两端封端；仅聚合物模式，默认 hydrogen")
    parser.add_argument("--boron-preprocess", choices=["strict", "relaxed", "projected", "skip"], default="strict",
                        help="硼平面化预处理：严格保键、放宽保键、几何投影、跳过；优化保留共面二面角冻结")
    parser.add_argument("--boron-bond-tolerance", type=float, default=0.15,
                        help="允许的最大键长变化，Å；适用于 relaxed 和 projected（后者上限 0.03）")
    parser.add_argument("--preprocess-method", choices=list(PREPROCESS_METHODS), default="auto",
                        help="整体三维预处理：auto、mmff94s、uff、xtb(GFN2)、tblite(GFN2/tight)、gfnff 或 none(仅ETKDG)")
    parser.add_argument("-o", "--output", type=Path, help="输出 .gjf 或 .com 文件")
    parser.add_argument("--method", default="B3LYP")
    parser.add_argument("--basis", default="6-31G(d)")
    parser.add_argument("--job", choices=["optfreq", "opt", "sp", "freq"], default="optfreq")
    parser.add_argument("--charge", type=int, default=None)
    parser.add_argument("--multiplicity", type=int, default=None)
    parser.add_argument("--memory", type=int, default=4, help="内存，GB")
    parser.add_argument("--processors", type=int, default=4)
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--extra", default="SCF=XQC Int=UltraFine")
    parser.add_argument("--title", default="Molecule")
    parser.add_argument("--checkpoint", help="不含扩展名的 ASCII 检查点名称")
    parser.add_argument("--metadata", type=Path, help="额外保存结构、参数与生成说明为 JSON")
    parser.add_argument("--force", action="store_true", help="允许覆盖已有输出文件")
    args = parser.parse_args(argv)
    if args.gui or args.molecule is None:
        run_gui()
        return 0
    checkpoint = args.checkpoint or default_checkpoint(args.molecule, args.mode,
                                                       args.degree if args.degree is not None else 5,
                                                       args.end_group or "hydrogen")
    engine = args.engine or ("orca" if args.output and args.output.suffix.lower() == ".inp" else "gaussian")
    output = args.output or Path(checkpoint + (".inp" if engine == "orca" else ".gjf"))
    try:
        if args.mode not in ("polymer",) and (args.degree is not None or args.end_group is not None):
            raise InputError("--degree 和 --end-group 需要同时指定 --mode polymer。")
        if output.suffix.lower() not in ((".inp",) if engine == "orca" else (".gjf", ".com")):
            raise InputError("Gaussian 使用 .gjf/.com，ORCA 使用 .inp。")
        if args.metadata and args.metadata.resolve() == Path(args.molecule).resolve() and (args.mode == "cif" or args.molecule.lower().endswith(".cif")):
            raise InputError("JSON 不可覆盖源 CIF。")
        if args.metadata and output.resolve() == args.metadata.resolve():
            raise InputError("输入文件和 JSON 元数据不能使用同一路径。")
        for path in (output, args.metadata):
            if path and path.exists() and not args.force:
                raise InputError(f"文件已存在：{path}。使用 --force 可覆盖。")
        settings = Settings(args.method, args.basis, args.job, args.memory, args.processors,
                            args.charge, args.multiplicity, args.conformers, args.extra,
                            args.title, checkpoint, args.degree if args.degree is not None else 5,
                            args.end_group or "hydrogen", args.boron_preprocess, args.boron_bond_tolerance,
                            args.preprocess_method, engine, args.cif_block, args.cif_component, args.cif_scope)
        result = generate(args.molecule, args.mode, settings)
        output.write_text(result.text, encoding="ascii", newline="\n")
        if args.metadata:
            metadata = asdict(result)
            metadata.pop("text")
            metadata["settings"] = asdict(settings)
            metadata["input"] = args.molecule
            metadata["mode"] = args.mode
            args.metadata.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"已生成：{output.resolve()}")
        print(f"{result.formula} | charge={result.charge} | multiplicity={result.multiplicity} | {result.smiles}")
        for note in result.warnings:
            print(f"提示：{note}")
        return 0
    except (InputError, OSError, RuntimeError) as exc:
        print(f"错误：{exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    # Legacy Windows consoles cannot encode every unit/symbol in diagnostics.
    # Keep their encoding, but never fail after successfully writing a file.
    for stream in (sys.stdout, sys.stderr):
        if stream is not None and hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    raise SystemExit(cli())
