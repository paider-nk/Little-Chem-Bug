"""Single-window application for generating, running and analysing jobs."""
from __future__ import annotations

import tkinter as tk
from tkinter import messagebox, ttk

from gaussian_input import run_gui
from tool_setup import ToolSettings
from tools_gui import ToolsPage
from workflow_gui import WorkflowApp, configure_tk_paths


class Workbench:
    def __init__(self, root):
        self.root = root
        root.title('ChemWorkbench · 计算化学工作台')
        screen_w, screen_h = root.winfo_screenwidth(), root.winfo_screenheight()
        root.geometry(f'{min(screen_w - 30, 1100)}x{min(screen_h - 85, 920)}')
        root.minsize(850, 570)
        ttk.Style(root).theme_use('clam')
        self.tabs = ttk.Notebook(root, padding=(6, 5))
        self.tabs.pack(fill='both', expand=True)
        self.input_page, self.workflow_page = ttk.Frame(self.tabs), ttk.Frame(self.tabs)
        self.settings = ToolSettings()
        self.workflow = WorkflowApp(root, self.workflow_page)
        self.tools = ToolsPage(self.tabs, self.settings, self.workflow.update_tool_path)
        self.tabs.add(self.input_page, text='生成输入')
        self.tabs.add(self.workflow_page, text='计算与后处理')
        self.tabs.add(self.tools, text='工具与安装')
        canvas = tk.Canvas(self.input_page, background='#f3f6fb', highlightthickness=0)
        scrollbar = ttk.Scrollbar(self.input_page, orient='vertical', command=canvas.yview)
        scrollbar.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        canvas.configure(yscrollcommand=scrollbar.set)
        editor = ttk.Frame(canvas)
        panel = canvas.create_window(0, 0, anchor='nw', window=editor)
        run_gui(editor, self.send_to_workflow)

        def resize(_=None):
            canvas.itemconfigure(panel, width=canvas.winfo_width(),
                                 height=max(canvas.winfo_height(), editor.winfo_reqheight()))
            canvas.configure(scrollregion=canvas.bbox(panel))
        canvas.bind('<Configure>', resize)
        editor.bind('<Configure>', resize)
        root.after_idle(resize)
        root.protocol('WM_DELETE_WINDOW', self.close)

    def send_to_workflow(self, path):
        self.workflow.use_input(path)
        self.tabs.select(self.workflow_page)

    def close(self):
        if self.tools.busy or self.workflow.busy:
            if not messagebox.askyesno('任务仍在运行', '关掉窗口会中断应用中的任务进度。仍要退出吗？', parent=self.root):
                return
        self.root.destroy()


def main():
    configure_tk_paths()
    root = tk.Tk()
    Workbench(root)
    root.mainloop()


if __name__ == '__main__':
    main()
