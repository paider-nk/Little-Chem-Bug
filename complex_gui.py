"""Graphical single-centre coordination-complex editor."""
import json
import tkinter as tk
from tkinter import messagebox, ttk

from complex_core import GEOMETRIES, PRESETS, load_spec, validate_spec


def open_complex_editor(parent, variables):
    window = tk.Toplevel(parent)
    window.title("配合物构建 · molSimplify")
    window.geometry("960x640")
    window.minsize(900, 610)
    frame = ttk.Frame(window, padding=18)
    frame.pack(fill="both", expand=True)
    frame.columnconfigure(0, weight=1)
    ttk.Label(frame, text="金属中心 + 配体 → 配合物三维结构", style="Title.TLabel").grid(row=0, sticky="w")
    preset = tk.StringVar(value=next(iter(PRESETS)))
    bar = ttk.Frame(frame)
    bar.grid(row=1, sticky="ew", pady=12)
    ttk.Label(bar, text="示例配方").pack(side="left")
    ttk.Combobox(bar, textvariable=preset, values=list(PRESETS), state="readonly", width=25).pack(side="left", padx=8)
    fields = {k: tk.StringVar() for k in ("metal", "oxidation_state", "geometry", "charge", "multiplicity")}
    geometry_labels = {f"{label} ({key}, {sites} 位点)": key for key, (label, sites) in GEOMETRIES.items()}
    details = ttk.Frame(frame)
    details.grid(row=2, sticky="ew")
    for index, (key, label, width) in enumerate((("metal", "金属", 6), ("oxidation_state", "氧化态", 5),
                                               ("geometry", "构型", 24), ("charge", "总电荷", 5),
                                               ("multiplicity", "多重度", 5))):
        ttk.Label(details, text=label).grid(row=0, column=index * 2, padx=(8 if index else 0, 4))
        if key == "geometry":
            widget = ttk.Combobox(details, textvariable=fields[key], values=list(geometry_labels), state="readonly", width=width)
        else:
            widget = ttk.Entry(details, textvariable=fields[key], width=width)
        widget.grid(row=0, column=index * 2 + 1)
    table = ttk.Frame(frame)
    table.grid(row=3, sticky="nsew", pady=15)
    table.columnconfigure(1, weight=1)
    for col, label in enumerate(("配体来源", "配体库名称 / SMILES", "数量", "配位原子（1起始）")):
        ttk.Label(table, text=label).grid(row=0, column=col, sticky="w", padx=4, pady=5)
    rows = []
    for row in range(8):
        values = [tk.StringVar(value=v) for v in ("配体库", "", "1", "")]
        rows.append(values)
        ttk.Combobox(table, textvariable=values[0], values=("配体库", "SMILES"), state="readonly", width=9).grid(row=row + 1, column=0, padx=4, pady=3)
        ttk.Entry(table, textvariable=values[1], width=45).grid(row=row + 1, column=1, sticky="ew", padx=4, pady=3)
        ttk.Spinbox(table, textvariable=values[2], from_=1, to=8, width=5).grid(row=row + 1, column=2, padx=4, pady=3)
        ttk.Entry(table, textvariable=values[3], width=21).grid(row=row + 1, column=3, padx=4, pady=3)
    ttk.Label(frame, wraplength=880, justify="left", text=(
        "配体库：ammonia（氨）、water（水）、bipy（联吡啶）、en（乙二胺）、chloride（氯）、cyclopentadienyl（Cp）。\n"
        "自定义示例：SMILES = NCCN，配位原子 = 1,4；数量 3，八面体构型，表示 3 个双齿乙二胺。\n"
        "编号按 SMILES 原子读取顺序，从 1 开始。普通配体按供体数占位；库内 π 配体按一个方向位点占位。\n"
        "空白行忽略。配方用于单金属初始结构；总电荷、多重度及目标构型需自行指定。"
    )).grid(row=4, sticky="w")

    def populate(spec, charge, multiplicity):
        fields["metal"].set(spec["metal"])
        fields["oxidation_state"].set(str(spec["oxidation_state"]))
        fields["geometry"].set(next(label for label, key in geometry_labels.items() if key == spec["geometry"]))
        fields["charge"].set(str(charge))
        fields["multiplicity"].set(str(multiplicity))
        for row in rows:
            for variable, value in zip(row, ("配体库", "", "1", "")):
                variable.set(value)
        for row, ligand in zip(rows, spec["ligands"]):
            for variable, value in zip(row, ("SMILES" if ligand["kind"] == "smiles" else "配体库", ligand["value"],
                                            str(ligand["count"]), ",".join(map(str, ligand.get("donors", []))))):
                variable.set(value)

    def apply():
        try:
            ligands = []
            for kind, value, count, donors in rows:
                if not value.get().strip():
                    continue
                ligand = {"kind": "smiles" if kind.get() == "SMILES" else "library", "value": value.get().strip(), "count": int(count.get())}
                if ligand["kind"] == "smiles":
                    ligand["donors"] = [int(i.strip()) for i in donors.get().replace("，", ",").split(",") if i.strip()]
                ligands.append(ligand)
            spec = validate_spec({"metal": fields["metal"].get().strip(), "oxidation_state": int(fields["oxidation_state"].get()),
                                  "geometry": geometry_labels[fields["geometry"].get()], "ligands": ligands})
            charge, mult = int(fields["charge"].get()), int(fields["multiplicity"].get())
            if not -20 <= charge <= 20 or not 1 <= mult <= 20:
                raise ValueError("总电荷范围 -20–20，多重度范围 1–20。")
        except (ValueError, KeyError) as exc:
            messagebox.showerror("配方有误", str(exc), parent=window)
            return
        variables["complex_spec"].set(json.dumps(spec, ensure_ascii=False))
        variables["input"].set(spec["metal"] + " + " + " + ".join(f"{li['count']} × {li['value']}" for li in spec["ligands"]))
        variables["charge"].set(str(charge))
        variables["multiplicity"].set(str(mult))
        variables["basis"].set("Def2SVP")
        variables["mode"].set("配合物 molSimplify")
        window.destroy()

    ttk.Button(bar, text="载入示例", command=lambda: populate(*PRESETS[preset.get()])).pack(side="left")
    ttk.Button(frame, text="应用配方，返回生成界面", command=apply, style="Accent.TButton").grid(row=5, sticky="e", pady=12)
    populate(*PRESETS[preset.get()])
    if variables["complex_spec"].get():
        try:
            populate(load_spec(variables["complex_spec"].get()), variables["charge"].get(), variables["multiplicity"].get())
        except ValueError:
            pass
    window.transient(parent)
    window.grab_set()
