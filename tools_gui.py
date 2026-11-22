"""Tools and installation page of the unified workbench."""
from pathlib import Path
import os
import queue
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import webbrowser

from tool_setup import (TOOLS, PROBE_MOLSIMPLIFY, checked_run, data_dir, detect_tool,
                        find_in_folder, install_molsimplify, install_zip)


class ToolsPage(ttk.Frame):
    def __init__(self, parent, settings, on_change):
        super().__init__(parent, padding=20)
        self.settings = settings
        self.on_change = on_change
        self.events = queue.Queue()
        self.busy = False
        self.started_at = 0.0
        self.stage = ''
        self.table = ttk.Treeview(self, columns=('state', 'role'), show='tree headings', height=7, selectmode='browse')
        self.table.heading('#0', text='工具')
        self.table.heading('state', text='状态')
        self.table.heading('role', text='用途')
        self.table.column('#0', width=165, stretch=False)
        self.table.column('state', width=110, stretch=False)
        self.table.column('role', width=460)
        self.table.pack(fill='x')
        for key, (name, _, _, role) in TOOLS.items():
            self.table.insert('', 'end', iid=key, text=name, values=('', role))

        self.path = tk.StringVar()
        ttk.Label(self, textvariable=self.path, wraplength=820).pack(anchor='w', pady=(12, 8))
        self.actions = ttk.Frame(self)
        self.actions.pack(fill='x')
        entries = [('选择程序', self.pick_file), ('从目录查找', self.pick_dir), ('打开官网', self.official),
                   ('安装本地包', self.local_install), ('自动安装', self.auto_install), ('清除设置', self.clear)]
        self.buttons = {}
        for name, callback in entries:
            self.buttons[name] = ttk.Button(self.actions, text=name, command=callback)
            self.buttons[name].pack(side='left', padx=(0, 5))
        self.status = tk.StringVar(value='选择工具，连接已有程序或安装所需组件。')
        ttk.Label(self, textvariable=self.status, wraplength=820).pack(anchor='w', pady=(12, 6))
        ttk.Label(self, text='xTB 与 molSimplify 支持自动安装。使用其他工具时，可从官网取得 Windows 安装包，或连接已有目录。',
                  wraplength=820).pack(anchor='w', pady=(10, 6))
        ttk.Button(self, text='打开配置与计算目录', command=self.open_data).pack(anchor='w', pady=4)
        self.table.bind('<<TreeviewSelect>>', lambda _: self.refresh())
        self.table.selection_set('XTB_EXE')
        self.refresh()
        self.after(100, self.poll)

    def selected(self):
        selection = self.table.selection()
        return selection[0] if selection else 'XTB_EXE'

    def refresh(self):
        for key, (_, _, _, role) in TOOLS.items():
            path = detect_tool(key)
            status = '已找到' if path else ('路径失效' if self.settings.paths.get(key) else '待配置')
            self.table.item(key, values=(status, role))
        key = self.selected()
        self.path.set(str(detect_tool(key) or self.settings.paths.get(key) or f'选择 {TOOLS[key][1]}'))
        choices = {'打开官网': bool(TOOLS[key][2]), '自动安装': key in ('XTB_EXE', 'MOLSIMPLIFY_PYTHON'),
                   '安装本地包': key != 'FORMCHK_EXE',
                   '清除设置': bool(self.settings.paths.get(key))}
        for name, button in self.buttons.items():
            button.configure(state='normal' if choices.get(name, True) and not self.busy else 'disabled')

    def save(self, key, path):
        self.settings.save(key, path)
        self.on_change(key)
        self.status.set(f'已连接 {TOOLS[key][0]}' if path else '已清除手动设置。')
        self.refresh()

    def use(self, key, path):
        if key == 'MOLSIMPLIFY_PYTHON':
            def probe(_):
                checked_run([path, '-c', PROBE_MOLSIMPLIFY])
                return path
            self.start(key, probe)
        else:
            try:
                self.save(key, path)
            except OSError as exc:
                messagebox.showerror('保存失败', str(exc), parent=self)

    def pick_file(self):
        key = self.selected()
        path = filedialog.askopenfilename(parent=self, title=f'选择 {TOOLS[key][1]}', filetypes=[('Windows 程序', '*.exe')])
        if path:
            if Path(path).name.casefold() != TOOLS[key][1].casefold():
                messagebox.showerror('程序名称不匹配', f'请选择 {TOOLS[key][1]}。', parent=self)
                return
            self.use(key, Path(path).resolve())

    def pick_dir(self):
        key = self.selected()
        folder = filedialog.askdirectory(parent=self, title=f'选择 {TOOLS[key][0]} 的目录')
        if folder:
            self.start(key, lambda _: self.find_and_probe(folder, key))

    @staticmethod
    def find_and_probe(folder, key):
        path = find_in_folder(folder, key)
        if key == 'MOLSIMPLIFY_PYTHON':
            checked_run([path, '-c', PROBE_MOLSIMPLIFY])
        return path

    def official(self):
        webbrowser.open(TOOLS[self.selected()][2])

    def clear(self):
        try:
            self.save(self.selected(), '')
        except OSError as exc:
            messagebox.showerror('保存失败', str(exc), parent=self)

    def open_data(self):
        try:
            data_dir().mkdir(parents=True, exist_ok=True)
            os.startfile(data_dir())
        except OSError as exc:
            messagebox.showerror('无法打开目录', str(exc), parent=self)

    def local_install(self):
        key = self.selected()
        if key == 'MOLSIMPLIFY_PYTHON':
            wheel = filedialog.askopenfilename(parent=self, title='选择 molSimplify 2.0.0 的官方 wheel',
                                                filetypes=[('Python wheel', '*.whl')])
            if wheel:
                self.start(key, lambda report: install_molsimplify(report, local_wheel=wheel))
            return
        picked = filedialog.askopenfilename(parent=self, title='选择从官网取得的 Windows 安装包',
                                             filetypes=[('安装包', '*.zip *.exe *.msi')])
        if not picked:
            return
        package = Path(picked).resolve()
        if package.suffix.lower() == '.zip':
            self.start(key, lambda report: install_zip(key, package, report))
        elif package.suffix.lower() in ('.exe', '.msi'):
            try:
                os.startfile(str(package), cwd=str(package.parent))
                self.status.set('请完成厂商安装向导，再用「从目录查找」选择安装目录。')
            except OSError as exc:
                messagebox.showerror('安装程序无法启动', str(exc), parent=self)

    def auto_install(self):
        key = self.selected()
        action = install_molsimplify if key == 'MOLSIMPLIFY_PYTHON' else lambda report: install_zip(key, report=report)
        self.start(key, action)

    def start(self, key, action):
        self.busy = True
        self.started_at = time.monotonic()
        self.stage = f'正在准备 {TOOLS[key][0]}…'
        self.status.set(self.stage)
        self.refresh()

        def worker():
            try:
                path = action(lambda stage: self.events.put(('progress', stage)))
                self.events.put(('done', (key, path)))
            except Exception as exc:
                self.events.put(('error', str(exc)))
        threading.Thread(target=worker, daemon=True).start()

    def poll(self):
        try:
            while True:
                event, value = self.events.get_nowait()
                if event == 'progress':
                    self.stage = value
                else:
                    self.busy = False
                    try:
                        if event == 'error':
                            raise RuntimeError(value)
                        self.save(*value)
                    except (RuntimeError, OSError) as exc:
                        self.status.set(f'{self.stage.rstrip("…")}。设置未完成。')
                        messagebox.showerror('工具设置失败', str(exc), parent=self)
                    self.refresh()
        except queue.Empty:
            pass
        if self.busy:
            elapsed = int(time.monotonic() - self.started_at)
            self.status.set(f'{self.stage} 已用 {elapsed} 秒')
        self.after(100, self.poll)
