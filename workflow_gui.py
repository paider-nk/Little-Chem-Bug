"""Tkinter interface for Gaussian execution and result post-processing."""
from __future__ import annotations

import json
import os
import queue
import sys
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from PIL import Image, ImageTk

from frontier_levels import SPIN_LABELS, generate_frontier_levels
from analysis_core import difference_density, energy_difference, parse_log, solvation_difference, write_report
from gaussian_workflow import run_jobs, run_multiwfn_analysis
from esp_map import generate_esp
from vmd_support import open_vmd
from mascot_ui import MascotDock
from tool_setup import data_dir, detect_tool


def configure_tk_paths():
    """Repair stale Tcl/Tk environment variables after Python was moved."""
    for env_name, version, filename, folder in (
        ("TCL_LIBRARY", tk.TclVersion, "init.tcl", "tcl"),
        ("TK_LIBRARY", tk.TkVersion, "tk.tcl", "tk"),
    ):
        current = os.environ.get(env_name)
        if current and (Path(current) / filename).is_file():
            continue
        candidates = []
        for prefix in (Path(sys.base_prefix), Path(sys.prefix), Path(sys.executable).parent):
            candidates.extend((prefix / "tcl" / f"{folder}{version}",
                               prefix / "Lib" / "tkinter" / folder))
        found = next((p for p in candidates if (p / filename).is_file()), None)
        if found:
            os.environ[env_name] = str(found)


class WorkflowApp:
    def __init__(self, root, parent=None):
        self.root = root
        if parent is None:
            root.title("Gaussian / ORCA 计算与后处理工作流")
            root.geometry("900x720")
            root.minsize(760, 620)
        self.events = queue.Queue()
        self.busy = False

        style = ttk.Style()
        if parent is None:
            try:
                style.theme_use("vista")
            except tk.TclError:
                pass
        style.configure("Title.TLabel", font=("Microsoft YaHei UI", 18, "bold"))
        style.configure("Hint.TLabel", foreground="#53657a")

        outer = ttk.Frame(parent if parent is not None else root, padding=18)
        outer.pack(fill="both", expand=True)
        ttk.Label(outer, text="Gaussian / ORCA 计算与后处理", style="Title.TLabel").pack(anchor="w")
        ttk.Label(outer, text="运行 Gaussian 16 或 ORCA、生成日志报告、计算能量差，并通过 Multiwfn 处理波函数。",
                  style="Hint.TLabel").pack(anchor="w", pady=(3, 12))
        self.mascot = MascotDock(outer, "选择 Gaussian 计算或后处理，小化区会在这里显示进度。")
        self.mascot.pack(side="bottom", fill="x", pady=(10, 0))
        self.tabs = ttk.Notebook(outer)
        self.tabs.pack(fill="both", expand=True)
        self._run_tab()
        self._analysis_tab()
        self._frontier_tab()
        self._energy_tab()
        self._density_tab()
        self._multiwfn_tab()
        self._esp_map_tab()
        self._vmd_tab()

        self.root.after(100, self._poll)

    def use_input(self, path):
        self.run_inputs.set(str(Path(path).resolve()))
        self.tabs.select(0)

    def update_tool_path(self, key):
        variables = {"GAUSSIAN_EXE": self.gaussian, "ORCA_EXE": self.orca, "VMD_EXE": self.vmd_exe}
        if key in variables:
            variables[key].set(str(detect_tool(key) or ""))

    def _tab(self, name):
        frame = ttk.Frame(self.tabs, padding=16)
        self.tabs.add(frame, text=name)
        return frame

    def _row(self, parent, row, label, variable, button_text=None, command=None):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", pady=6)
        entry = ttk.Entry(parent, textvariable=variable)
        entry.grid(row=row, column=1, sticky="ew", padx=8, pady=6)
        if button_text:
            ttk.Button(parent, text=button_text, command=command).grid(row=row, column=2, pady=6)
        parent.columnconfigure(1, weight=1)
        return entry

    def _pick_files(self, variable, types, multiple=False):
        paths = filedialog.askopenfilenames(filetypes=types) if multiple else filedialog.askopenfilename(filetypes=types)
        if paths:
            variable.set(";".join(paths) if multiple else paths)

    def _pick_dir(self, variable):
        path = filedialog.askdirectory()
        if path:
            variable.set(path)

    def _run_tab(self):
        frame = self._tab("运行计算")
        self.run_inputs = tk.StringVar()
        self.run_output = tk.StringVar(value=str(data_dir() / "Runs"))
        self.gaussian = tk.StringVar(value=str(detect_tool("GAUSSIAN_EXE") or ""))
        self.cores = tk.IntVar(value=4)
        self.memory = tk.IntVar(value=4)
        self.slots = tk.IntVar(value=1)
        self._row(frame, 0, "输入文件（可多选）", self.run_inputs, "选择…",
                  lambda: self._pick_files(self.run_inputs, [("Gaussian / ORCA", "*.gjf *.com *.inp")], True))
        self._row(frame, 1, "结果目录", self.run_output, "选择…", lambda: self._pick_dir(self.run_output))
        self._row(frame, 2, "Gaussian 程序", self.gaussian, "浏览…",
                  lambda: self._pick_files(self.gaussian, [("Gaussian executable", "g16*.exe"), ("All files", "*.*")]))
        self.orca = tk.StringVar(value=str(detect_tool("ORCA_EXE") or ""))
        self._row(frame, 6, "ORCA 程序", self.orca, "浏览…", lambda: self._pick_files(self.orca, [("ORCA", "orca.exe")]))
        resources = ttk.Frame(frame)
        resources.grid(row=3, column=1, sticky="w", padx=8, pady=6)
        for label, var in (("每任务核心", self.cores), ("内存 GB", self.memory), ("并行槽位", self.slots)):
            ttk.Label(resources, text=label).pack(side="left", padx=(0, 4))
            ttk.Spinbox(resources, from_=1, to=256, width=5, textvariable=var).pack(side="left", padx=(0, 14))
        ttk.Label(frame, text="输入文件中的方法、基组、电荷、多重度和坐标会保留。总核心数 = 每任务核心 × 并行槽位。",
                  style="Hint.TLabel", wraplength=680).grid(row=4, column=1, sticky="w", padx=8)
        ttk.Button(frame, text="开始计算", command=self._start_run).grid(row=5, column=1, sticky="w", padx=8, pady=18)

    def _analysis_tab(self):
        frame = self._tab("日志报告")
        self.log_path, self.report_dir = tk.StringVar(), tk.StringVar()
        self._row(frame, 0, "计算日志", self.log_path, "选择…",
                  lambda: self._pick_files(self.log_path, [("Gaussian / ORCA", "*.log *.out"), ("All files", "*.*")]))
        self._row(frame, 1, "报告目录（可选）", self.report_dir, "选择…", lambda: self._pick_dir(self.report_dir))
        ttk.Label(frame, text="输出 HTML、JSON、CSV、能级/电荷/IR 图和最终结构 XYZ（数据存在时）。",
                  style="Hint.TLabel", wraplength=680).grid(row=2, column=1, sticky="w", padx=8, pady=8)
        ttk.Button(frame, text="生成分析报告", command=self._analyze).grid(row=3, column=1, sticky="w", padx=8, pady=12)

    def _frontier_tab(self):
        frame = self._tab("HOMO / LUMO")
        self.frontier_source = tk.StringVar()
        ttk.Label(frame, text="计算结果").grid(row=0, column=0, sticky="w", pady=6)
        self.frontier_recent = ttk.Combobox(frame, textvariable=self.frontier_source)
        self.frontier_recent.grid(row=0, column=1, sticky="ew", padx=8, pady=6)
        result_types = [("Gaussian / ORCA 结果", "*.log *.out *.chk *.fchk *.fch *.gbw *.molden *.molden.input")]
        ttk.Button(frame, text="选择…", command=lambda: self._pick_files(self.frontier_source, result_types)).grid(row=0, column=2)
        self.frontier_output = tk.StringVar()
        self._row(frame, 1, "保存位置（可选）", self.frontier_output, "选择…", lambda: self._pick_dir(self.frontier_output))
        actions = ttk.Frame(frame)
        actions.grid(row=2, column=0, columnspan=3, sticky="w", pady=(4, 6))
        ttk.Button(actions, text="提取能级并绘图", command=self._generate_frontier).pack(side="left")
        self.frontier_open = ttk.Button(actions, text="打开结果目录", state="disabled",
                                       command=lambda: os.startfile(self.frontier_folder))
        self.frontier_open.pack(side="left", padx=8)
        self.frontier_full_image = ttk.Button(actions, text="打开能级图", state="disabled",
                                             command=lambda: os.startfile(self.frontier_folder / 'orbital_levels.png'))
        self.frontier_full_image.pack(side="left")
        results = ttk.Frame(frame)
        results.grid(row=3, column=0, columnspan=3, sticky="nsew")
        results.columnconfigure(1, weight=1)
        results.rowconfigure(1, weight=1)
        self.frontier_table = ttk.Treeview(results, columns=("spin", "homo", "lumo", "gap"), show="headings", height=2)
        for key, label in zip(("spin", "homo", "lumo", "gap"), ("自旋通道", "HOMO (eV)", "LUMO (eV)", "能隙 (eV)")):
            self.frontier_table.heading(key, text=label)
            self.frontier_table.column(key, width=75, anchor="center")
        self.frontier_table.grid(row=0, column=0, sticky="ew")
        self.frontier_note = ttk.Label(results, style="Hint.TLabel", wraplength=280, justify="left")
        self.frontier_note.grid(row=1, column=0, sticky="nw", pady=(6, 0))
        self.frontier_preview = tk.Canvas(results, width=1, height=1, background="white", highlightthickness=0)
        self.frontier_image_id = self.frontier_preview.create_image(0, 0)
        self.frontier_preview.grid(row=0, column=1, rowspan=2, sticky="nsew", padx=(10, 0))
        frame.rowconfigure(3, weight=1)
        self.frontier_image = None
        self.frontier_preview.bind("<Configure>", self._resize_frontier)

    def _resize_frontier(self, event=None):
        if self.frontier_image is None:
            return
        width = self.frontier_preview.winfo_width() if event is None else event.width
        height = self.frontier_preview.winfo_height() if event is None else event.height
        image = self.frontier_image.copy()
        image.thumbnail((max(1, width-4), max(1, height-4)), Image.Resampling.LANCZOS)
        self.frontier_photo = ImageTk.PhotoImage(image)
        self.frontier_preview.coords(self.frontier_image_id, width / 2, height / 2)
        self.frontier_preview.itemconfigure(self.frontier_image_id, image=self.frontier_photo)

    def _generate_frontier(self):
        if self.busy:
            return
        source = self.frontier_source.get().strip()
        if not source:
            messagebox.showwarning("需要计算结果", "请选择计算日志、检查点或波函数文件。")
            return
        output = self.frontier_output.get().strip() or None
        gaussian, orca = self.gaussian.get().strip() or None, self.orca.get().strip() or None
        self.frontier_table.delete(*self.frontier_table.get_children())
        self.frontier_open.configure(state="disabled")
        self.frontier_full_image.configure(state="disabled")
        self.frontier_note.configure(text="")
        self.frontier_preview.itemconfigure(self.frontier_image_id, image="")
        self.frontier_image = None
        self._submit("HOMO / LUMO", lambda stage: generate_frontier_levels(source, output, gaussian=gaussian, orca=orca,
                                                                           on_stage=stage))

    def _show_frontier(self, result):
        self.frontier_folder = result['folder']
        self.frontier_open.configure(state="normal")
        self.frontier_full_image.configure(state="normal")
        for row in result['rows']:
            values = [f"{row[key]:.4f}" if row[key] is not None else "未输出" for key in ('homo_ev', 'lumo_ev', 'gap_ev')]
            self.frontier_table.insert("", "end", values=[SPIN_LABELS[row['spin']], *values])
        warning = result['warning']
        self.frontier_note.configure(text=warning)
        with Image.open(result['png']) as image:
            self.frontier_image = image.copy()
        self._resize_frontier()
        self.mascot.set_state("warning" if warning else "done", "轨道数值和能级图已保存", "可打开原图或结果目录查看导出文件。")

    def _energy_tab(self):
        frame = self._tab("结合能 / 溶剂化")
        self.binding_logs, self.gas_log, self.solution_log = tk.StringVar(), tk.StringVar(), tk.StringVar()
        self.energy_free = tk.BooleanVar(value=False)
        self._row(frame, 0, "结合能日志（按 AB、片段 A、片段 B 顺序）", self.binding_logs, "多选日志…",
                  lambda: self._pick_files(self.binding_logs, [("Gaussian / ORCA", "*.log *.out")], True))
        self._row(frame, 1, "气相日志", self.gas_log, "选择…",
                  lambda: self._pick_files(self.gas_log, [("Gaussian / ORCA", "*.log *.out")]))
        self._row(frame, 2, "溶液日志", self.solution_log, "选择…",
                  lambda: self._pick_files(self.solution_log, [("Gaussian / ORCA", "*.log *.out")]))
        ttk.Checkbutton(frame, text="使用 Gibbs 自由能（需要兼容的频率热化学结果）", variable=self.energy_free).grid(
            row=3, column=1, sticky="w", padx=8, pady=8)
        ttk.Label(frame, text="分别计算结合能和溶剂化差；结合能对话框中选择顺序必须是复合物 AB、片段 A、片段 B……。",
                  style="Hint.TLabel", wraplength=680).grid(row=4, column=1, sticky="w", padx=8, pady=8)
        buttons = ttk.Frame(frame)
        buttons.grid(row=5, column=1, sticky="w", padx=8, pady=10)
        ttk.Button(buttons, text="计算结合能", command=self._binding).pack(side="left", padx=(0, 10))
        ttk.Button(buttons, text="计算溶剂化差", command=self._solvation).pack(side="left")

    def _density_tab(self):
        frame = self._tab("电子密度差")
        self.cube_ab, self.cube_a, self.cube_b, self.cube_out = (tk.StringVar() for _ in range(4))
        for row, label, variable in ((0, "AB 密度 cube", self.cube_ab), (1, "片段 A 冻结几何 cube", self.cube_a),
                                     (2, "片段 B 冻结几何 cube", self.cube_b)):
            self._row(frame, row, label, variable, "选择…",
                      lambda v=variable: self._pick_files(v, [("Cube file", "*.cube *.cub")]))
        self._row(frame, 3, "差分 cube 输出", self.cube_out, "保存为…", self._save_cube)
        ttk.Label(frame, text="要求三份 cube 使用相同网格；两个片段的原子和坐标必须恰好拼成 AB。结果为 ρ(AB)-ρ(A)-ρ(B)，可在 Multiwfn 中可视化。",
                  style="Hint.TLabel", wraplength=680).grid(row=4, column=1, sticky="w", padx=8, pady=8)
        ttk.Button(frame, text="计算密度差", command=self._density).grid(row=5, column=1, sticky="w", padx=8, pady=12)

    def _save_cube(self):
        path = filedialog.asksaveasfilename(defaultextension=".cube", filetypes=[("Cube file", "*.cube")])
        if path:
            self.cube_out.set(path)

    def _multiwfn_tab(self):
        frame = self._tab("Multiwfn / ESP")
        self.wavefunction = tk.StringVar()
        self._row(frame, 0, "波函数文件", self.wavefunction, "选择…",
                  lambda: self._pick_files(self.wavefunction,
                                           [("Wavefunction", "*.chk *.fchk *.fch *.gbw *.molden *.molden.input"), ("All files", "*.*")]))
        ttk.Label(frame, text="自动转换 Gaussian .chk 和 ORCA .gbw。下面的分析会自动操作 Multiwfn，并把完整输出保存为文本供查看。",
                  style="Hint.TLabel", wraplength=680).grid(row=1, column=1, sticky="w", padx=8, pady=8)
        self.esp_output = tk.StringVar()
        self._row(frame, 2, "结果文本（可选）", self.esp_output, "保存为…", self._save_text)
        ttk.Label(frame, text="ESP 表面分析给出范德华面极值、平均值、方差及正负区域统计。网格极值包含核附近网格点，最大值可能受核处势能奇点影响，应谨慎解释。",
                  style="Hint.TLabel", wraplength=680).grid(row=3, column=1, sticky="w", padx=8, pady=8)
        buttons = ttk.Frame(frame)
        buttons.grid(row=4, column=1, sticky="w", padx=8, pady=12)
        ttk.Button(buttons, text="ESP 表面定量分析", command=lambda: self._multiwfn("esp_surface")).pack(side="left", padx=(0, 10))
        ttk.Button(buttons, text="ESP 网格极值（含核附近点）", command=lambda: self._multiwfn("esp_extrema")).pack(side="left")

    def _esp_map_tab(self):
        frame = self._tab("ESP 表面图")
        self._row(frame, 0, "波函数文件", self.wavefunction, "选择…",
                  lambda: self._pick_files(self.wavefunction, [("Gaussian / ORCA 波函数", "*.chk *.fchk *.fch *.gbw *.molden *.molden.input")]))
        self.map_folder = tk.StringVar()
        self._row(frame, 1, "保存位置（可选）", self.map_folder, "选择…", lambda: self._pick_dir(self.map_folder))
        self.map_iso = tk.StringVar(value="0.001")
        self.map_limit = tk.StringVar(value="0.05")
        self._row(frame, 2, "密度等值（e/bohr³）", self.map_iso)
        self._row(frame, 3, "色标范围 ±（a.u.）", self.map_limit)
        ttk.Label(frame, text="网格精度").grid(row=4, column=0, sticky="w")
        self.map_grid = ttk.Combobox(frame, values=("低", "中", "高"), state="readonly", width=9)
        self.map_grid.current(1)
        self.map_grid.grid(row=4, column=1, sticky="w", padx=8, pady=6)
        options = ttk.Frame(frame)
        options.grid(row=5, column=1, sticky="w", padx=8)
        self.map_legend = tk.BooleanVar(value=True)
        self.map_viewer = tk.BooleanVar(value=True)
        ttk.Checkbutton(options, text="图片附色标", variable=self.map_legend).pack(side="left")
        ttk.Checkbutton(options, text="生成后在 VMD 中旋转查看", variable=self.map_viewer).pack(side="left", padx=8)
        ttk.Label(frame, text="红色为负势，青蓝色为正势。默认在波函数旁新建目录，保存 PNG、TIFF 和两份 Cube。",
                  wraplength=540, style="Hint.TLabel").grid(row=6, column=1, sticky="w", padx=8, pady=8)
        actions = ttk.Frame(frame)
        actions.grid(row=7, column=1, sticky="w", padx=8)
        ttk.Button(actions, text="生成静电势表面图", command=self._generate_map).pack(side="left")
        self.map_open = ttk.Button(actions, text="打开图片目录", state="disabled", command=lambda: os.startfile(self.last_map_folder))
        self.map_open.pack(side="left", padx=8)

    def _generate_map(self):
        source = self.wavefunction.get().strip()
        if not source:
            messagebox.showwarning("需要波函数", "请选择 Gaussian 或 ORCA 波函数。")
            return
        try:
            iso, limit = float(self.map_iso.get()), float(self.map_limit.get())
        except ValueError:
            messagebox.showwarning("绘图参数", "密度等值和色标范围需填写正数。")
            return
        options = dict(iso=iso, limit=limit, quality=self.map_grid.current()+1,
                       legend=self.map_legend.get(), interactive=self.map_viewer.get(),
                       gaussian=self.gaussian.get(), orca=self.orca.get(), vmd=self.vmd_exe.get() or None)
        output = self.map_folder.get().strip() or None
        self._submit("ESP 表面图", lambda stage: generate_esp(source, output, on_stage=stage, **options))

    def _vmd_tab(self):
        frame = self._tab("VMD 可视化")
        self.vmd_source = tk.StringVar()
        self.vmd_exe = tk.StringVar(value=str(detect_tool("VMD_EXE") or ""))
        self.vmd_auto = tk.BooleanVar(value=True)
        self.vmd_recent = ttk.Combobox(frame, textvariable=self.vmd_source, state="readonly")
        ttk.Label(frame, text="本批已完成的计算").grid(row=5, column=0, sticky="w")
        self.vmd_recent.grid(row=5, column=1, sticky="ew", padx=8)
        self._row(frame, 0, "结果文件", self.vmd_source, "选择…",
                  lambda: self._pick_files(self.vmd_source, [("计算结果", "*.log *.out *.gbw *.xyz *.pdb *.cube *.cub *.molden *.molden.input"), ("All files", "*.*")]))
        self._row(frame, 1, "VMD 程序", self.vmd_exe, "浏览…",
                  lambda: self._pick_files(self.vmd_exe, [("VMD", "vmd.exe"), ("All files", "*.*")]))
        ttk.Checkbutton(frame, text="单个任务正常结束后自动用 VMD 打开最终结构", variable=self.vmd_auto).grid(
            row=2, column=1, sticky="w", padx=8, pady=8)
        ttk.Label(frame, text=".log/.out 会先提取最终坐标为 .final.xyz；.gbw 会先转为 Molden。Cube 默认显示 +0.01（蓝）和 −0.01（红）等值面。",
                  style="Hint.TLabel", wraplength=680).grid(row=3, column=1, sticky="w", padx=8, pady=8)
        ttk.Button(frame, text="用 VMD 打开", command=self._open_vmd).grid(row=4, column=1, sticky="w", padx=8, pady=12)

    def _open_vmd(self):
        source = self.vmd_source.get().strip()
        if not source:
            messagebox.showwarning("需要结果文件", "请选择计算日志、波函数或结构文件。")
            return
        executable = self.vmd_exe.get().strip() or None
        self._submit("VMD 可视化", lambda stage: open_vmd(source, vmd=executable))

    def _save_text(self):
        path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text file", "*.txt")])
        if path:
            self.esp_output.set(path)

    def _submit(self, label, function, *, batch=False):
        if self.busy:
            return
        self.busy = True
        self.mascot.set_state("running", label, "检查输入与运行环境…")

        def report(stage):
            self.events.put(("stage", stage))

        def report_counts(settled, total, abnormal):
            self.events.put(("batch", (settled, total, abnormal)))

        def worker():
            try:
                result = function(report_counts if batch else report)
                self.events.put(("result", (label, result, batch)))
            except Exception as exc:
                self.events.put(("error", (label, str(exc))))
        threading.Thread(target=worker, daemon=True).start()

    def _poll(self):
        try:
            while True:
                event, payload = self.events.get_nowait()
                if event == "stage":
                    self.mascot.set_state("running", payload, "请等待当前步骤完成。")
                elif event == "batch":
                    settled, total, abnormal = payload
                    detail = f"已结束 {settled} / {total} 项 · 正常 {settled - abnormal} · 异常 {abnormal}"
                    if settled == total:
                        detail += " · 正在写入批次报告"
                    self.mascot.set_state("running", "计算批量计算", detail, counts=(settled, total))
                elif event == "error":
                    self.busy = False
                    label, reason = payload
                    self.mascot.set_state("error", f"{label}失败", "请查看错误详情，检查输入或外部程序。")
                    messagebox.showerror(label, reason, parent=self.root)
                else:
                    self.busy = False
                    label, result, batch = payload
                    self._show_result(label, result, batch)
        except queue.Empty:
            pass
        self.root.after(100, self._poll)

    def _show_result(self, label, result, batch):
        if batch:
            total = len(result)
            abnormal = sum(job["status"] != "normal" for job in result)
            summary = f"已结束 {total} / {total} 项 · 正常 {total - abnormal} · 异常 {abnormal} · 批次报告已写入"
            self.mascot.set_state("warning" if abnormal else "done", "计算批次已结束",
                                  summary, counts=(total, total))
            logs = [job['log'] for job in result if job.get('status') == 'normal']
            self.frontier_recent.configure(values=logs)
            if logs:
                self.frontier_source.set(logs[0])
                self.log_path.set(logs[0])
                log = Path(logs[0])
                checkpoints = sorted(log.parent.glob('*.gbw' if log.suffix.casefold() == '.out' else '*.chk'))
                named = log.with_suffix('.gbw' if log.suffix.casefold() == '.out' else '.chk')
                wave = named if named.is_file() else checkpoints[0] if len(checkpoints) == 1 else None
                self.wavefunction.set(str(wave) if wave else '')
            self.vmd_recent.configure(values=logs)
            self.vmd_source.set(logs[0] if logs else '')
            if self.vmd_auto.get() and total == 1 and logs:
                self.root.after_idle(self._open_vmd)
                return
        elif label == "HOMO / LUMO":
            self._show_frontier(result)
            return
        elif label == "ESP 表面图":
            self.last_map_folder = result['folder']
            self.map_open.configure(state="normal")
            warning = result.get('viewer_warning')
            self.mascot.set_state("warning" if warning else "done", "ESP 表面图已保存", self.last_map_folder)
            info = f"PNG 和 TIFF 已保存至\n{self.last_map_folder}"
            if warning:
                messagebox.showwarning(label, info + "\nVMD 交互窗口未打开：" + warning, parent=self.root)
            else:
                messagebox.showinfo(label, info, parent=self.root)
            return
        elif label == '电子密度差':
            self.vmd_source.set(result['cube'])
            self.mascot.set_state('done', '差分密度已保存', '可在 VMD 可视化页直接打开。')
        elif label == "日志分析" and not result["normal_termination"]:
            self.mascot.set_state("warning", "日志报告已生成", "计算日志未记录正常结束，请查看报告。")
        else:
            self.mascot.set_state("done", f"{label}已完成", "结果已就绪，可在详情窗口查看。")
        detail_value = result
        detail = json.dumps(detail_value, ensure_ascii=False, indent=2) if isinstance(detail_value, (dict, list)) else str(detail_value)
        messagebox.showinfo(label, detail, parent=self.root)

    def _start_run(self):
        if self.busy:
            return
        inputs = [p for p in self.run_inputs.get().split(";") if p]
        output, executable = self.run_output.get(), self.gaussian.get()
        orca = self.orca.get()
        try:
            cores, memory, slots = self.cores.get(), self.memory.get(), self.slots.get()
        except tk.TclError:
            self.mascot.set_state("error", "运行参数有误", "每任务核心数、内存和并行槽位须填写整数。")
            return
        self._submit("量化计算", lambda report: run_jobs(inputs, output, gaussian=executable, orca=orca,
                                                                 cores=cores, memory=memory, slots=slots,
                                                                 on_progress=report), batch=True)

    def _analyze(self):
        path = self.log_path.get().strip()
        if not path:
            messagebox.showwarning("需要日志", "请选择 计算日志文件。")
            return
        out = self.report_dir.get().strip() or str(Path(path).resolve().parent / "analysis")
        def action(stage):
            stage("解析 计算日志")
            result = parse_log(path)
            stage("生成分析报告及图表")
            report = write_report(result, out)
            return {"report": str(report), "normal_termination": result.normal,
                    "energy_hartree": result.energy_hartree, "warnings": result.warnings}
        self._submit("日志分析", action)

    def _binding(self):
        paths = [p for p in self.binding_logs.get().split(";") if p]
        if len(paths) < 3:
            messagebox.showwarning("需要日志", "结合能至少选择 AB 和两个片段日志。")
            return
        free_energy = self.energy_free.get()
        def action(stage):
            stage("读取复合物及片段日志")
            complex_log, *fragments = [parse_log(path) for path in paths]
            stage("计算结合能")
            return energy_difference(complex_log, fragments, free_energy=free_energy)
        self._submit("结合能", action)

    def _solvation(self):
        gas, solution = self.gas_log.get().strip(), self.solution_log.get().strip()
        if not gas or not solution:
            messagebox.showwarning("需要日志", "请分别选择气相和溶液 计算日志。")
            return
        free_energy = self.energy_free.get()
        def action(stage):
            stage("读取气相与溶液日志")
            gas_log, solution_log = parse_log(gas), parse_log(solution)
            stage("计算溶剂化差")
            return solvation_difference(gas_log, solution_log, free_energy=free_energy)
        self._submit("溶剂化差", action)

    def _density(self):
        values = [self.cube_ab.get().strip(), self.cube_a.get().strip(), self.cube_b.get().strip(), self.cube_out.get().strip()]
        if not all(values):
            messagebox.showwarning("需要文件", "请指定三份 cube 输入和差分 cube 输出路径。")
            return
        def action(stage):
            stage("校验 cube 网格并计算密度差")
            return {**difference_density(*values), "cube": str(Path(values[-1]).resolve())}
        self._submit("电子密度差", action)

    def _multiwfn(self, analysis):
        path = self.wavefunction.get().strip()
        if not path:
            messagebox.showwarning("需要波函数", "请选择 Gaussian 或 ORCA 波函数。")
            return
        output = self.esp_output.get().strip() or None
        self._submit("Multiwfn ESP 分析", lambda stage: run_multiwfn_analysis(path, analysis, output, on_stage=stage))


def main():
    configure_tk_paths()
    root = tk.Tk()
    WorkflowApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
