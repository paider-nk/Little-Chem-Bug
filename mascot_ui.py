"""Shared animated chemistry companion for the two Tkinter windows."""
from __future__ import annotations

import math
import tkinter as tk
from pathlib import Path

NAME = "小化区"
CARD = "#fcfaff"
INK = "#2b2349"
MUTED = "#605677"
LINE = "#e0d8f1"
TRACK = "#ebe4f7"
FONT = "Microsoft YaHei UI"
STATES = {
    "idle": ("待命", "#61547d", "#eee8fa"),
    "running": ("进行中", "#6b45b6", "#ede4fb"),
    "done": ("已完成", "#226450", "#ddf4ec"),
    "warning": ("有异常", "#85521a", "#fff0d5"),
    "error": ("失败", "#a43757", "#fce5ee"),
}
SCENE_CENTER = (161, 78)
ORBIT_RADIUS = (126, 55)
SCENE_INTERVAL = 60
ORBIT_TICKS = 250
BLINK_TICKS = 75
MOTIFS = (
    ("benzene", "#f8e9f5", "#754889"),
    ("cell", "#e3f3f0", "#2d736c"),
    ("chromatogram", "#fcece5", "#975435"),
    ("orbital", "#eae7fc", "#624ca1"),
    ("coordination", "#e7eefc", "#426596"),
    ("energy", "#f9edf9", "#8a548c"),
)


class MascotDock(tk.Frame):
    """Show real task stages, and animate the portrait independently of task progress."""

    def __init__(self, parent, idle_hint):
        super().__init__(parent, bg=CARD, highlightbackground=LINE, highlightthickness=1)
        self.columnconfigure(1, weight=1)
        self.art = tk.Canvas(self, width=322, height=156, bg="#f1ecfa", highlightthickness=0)
        self.art.grid(row=0, column=0, sticky="nsw", padx=(5, 0), pady=5)
        self.collapsed = False
        self.motion = True
        self.frames = self._load_frames()
        self.state, self.counts = "idle", None
        self._progress_phase = 0
        self._scene_tick = 0
        self._progress_timer = None
        self._scene_timer = None

        body = tk.Frame(self, bg=CARD)
        body.grid(row=0, column=1, sticky="ew", padx=(16, 14), pady=8)
        body.columnconfigure(0, weight=1)
        header = tk.Frame(body, bg=CARD)
        header.grid(row=0, column=0, sticky="ew")
        tk.Label(header, text=NAME, font=(FONT, 11, "bold"), fg=INK, bg=CARD).pack(side="left")
        self.badge = tk.Label(header, font=(FONT, 9, "bold"), padx=7, pady=0)
        self.badge.pack(side="left", padx=(9, 0))
        self.toggle = self._header_button(header, "收起立绘", self._toggle)
        self.motion_button = self._header_button(header, "暂停动画", self._toggle_motion)
        self.title = tk.Label(body, font=(FONT, 13, "bold"), bg=CARD, fg=INK, anchor="w")
        self.title.grid(row=1, column=0, sticky="ew", pady=(3, 0))
        self.detail = tk.Label(body, font=(FONT, 9), bg=CARD, fg=MUTED, justify="left", anchor="w")
        self.detail.grid(row=2, column=0, sticky="ew", pady=(0, 5))
        self.track = tk.Canvas(body, height=5, width=1, bg=CARD, highlightthickness=0)
        self.track.grid(row=3, column=0, sticky="ew")
        self.track.bind("<Configure>", lambda _: self._draw_track())
        body.bind("<Configure>", lambda event: self.detail.configure(wraplength=max(200, event.width - 1)))
        self._draw_art()
        self._update_motion()
        self.set_state("idle", "等待任务", idle_hint)

    def _header_button(self, parent, label, command):
        button = tk.Button(parent, text=label, command=command, font=(FONT, 9),
                           bg=CARD, fg=MUTED, activebackground=TRACK, activeforeground=INK,
                           relief="flat", bd=0, highlightthickness=1, highlightbackground=LINE,
                           cursor="hand2", padx=5, pady=0)
        button.pack(side="right", padx=(6, 0))
        return button

    def _load_frames(self):
        assets = Path(__file__).with_name("assets")
        try:
            return [tk.PhotoImage(master=self, file=str(assets / f"xiaohuaqu-{name}.png"))
                    for name in ("open", "half", "closed")]
        except tk.TclError:
            try:
                return [tk.PhotoImage(master=self, file=str(assets / "chromia-104.png"))]
            except tk.TclError:
                return []

    def _toggle(self):
        self.collapsed = not self.collapsed
        if self.collapsed:
            self.art.grid_remove()
            self.motion_button.pack_forget()
        else:
            self.art.grid()
            self.motion_button.pack(side="right", padx=(6, 0))
        self.toggle.configure(text="展开立绘" if self.collapsed else "收起立绘")
        self._update_motion()

    def _toggle_motion(self):
        self.motion = not self.motion
        self.motion_button.configure(text="暂停动画" if self.motion else "播放动画")
        self._update_motion()

    def _update_motion(self):
        if self.motion and not self.collapsed:
            if self._scene_timer is None:
                self._scene_timer = self.after(SCENE_INTERVAL, self._animate_scene)
        elif self._scene_timer is not None:
            self.after_cancel(self._scene_timer)
            self._scene_timer = None
            if self.frames:
                self.art.itemconfigure(self.portrait, image=self.frames[0])

    def _animate_scene(self):
        loop_ticks = math.lcm(ORBIT_TICKS, BLINK_TICKS)
        self._scene_tick = (self._scene_tick + 1) % loop_ticks
        self._render_scene(self._scene_tick)
        self._scene_timer = self.after(SCENE_INTERVAL, self._animate_scene)

    def _render_scene(self, tick):
        center_x, center_y = SCENE_CENTER
        radius_x, radius_y = ORBIT_RADIUS
        angle = math.tau * tick / ORBIT_TICKS
        for index, (tag, previous) in enumerate(self._motif_positions):
            # Start with three motifs on each side of the portrait.
            theta = angle + math.tau * (index / len(MOTIFS) + 1 / 24)
            x = center_x + radius_x * math.cos(theta)
            y = center_y + radius_y * math.sin(theta)
            self.art.move(tag, x - previous[0], y - previous[1])
            self._motif_positions[index] = (tag, (x, y))
        if self.frames:
            phase = tick % BLINK_TICKS
            offset = 3 * math.sin(phase * math.tau / (BLINK_TICKS / 2))
            blink = {65: 1, 66: 2, 67: 2, 68: 1}.get(phase, 0) if len(self.frames) == 3 else 0
            self.art.coords(self.portrait, center_x, center_y + offset)
            self.art.itemconfigure(self.portrait, image=self.frames[blink])

    def set_state(self, state, title, detail="", *, counts=None):
        self.state, self.counts = state, counts
        label, ink, fill = STATES[state]
        self.badge.configure(text=label, fg=ink, bg=fill)
        self.title.configure(text=title)
        self.detail.configure(text=detail)
        self._draw_track()
        animate_progress = state == "running" and counts is None
        if animate_progress and self._progress_timer is None:
            self._progress_timer = self.after(70, self._animate_progress)
        elif not animate_progress and self._progress_timer is not None:
            self.after_cancel(self._progress_timer)
            self._progress_timer = None

    def _animate_progress(self):
        self._progress_phase = (self._progress_phase + 1) % 42
        self._draw_track()
        self._progress_timer = self.after(70, self._animate_progress)

    def destroy(self):
        for timer in (self._progress_timer, self._scene_timer):
            if timer is not None:
                self.after_cancel(timer)
        super().destroy()

    def _draw_track(self):
        bar = self.track
        bar.delete("all")
        width = bar.winfo_width()
        bar.create_rectangle(0, 0, width, 5, outline="", fill=TRACK)
        color = STATES[self.state][1]
        if self.counts:
            done, total = self.counts
            if total > 0 and done > 0:
                bar.create_rectangle(0, 0, width * min(done / total, 1), 5, outline="", fill=color)
        elif self.state == "running":
            fraction = self._progress_phase / 21
            x = min(fraction, 2 - fraction) * width * .73
            bar.create_rectangle(x, 0, x + width * .27, 5, outline="", fill=color)
        elif self.state == "done":
            bar.create_rectangle(0, 0, width, 5, outline="", fill=color)

    def _draw_art(self):
        art = self.art
        center_x, center_y = SCENE_CENTER
        rx, ry = ORBIT_RADIUS
        orbit = (center_x - rx, center_y - ry, center_x + rx, center_y + ry)
        art.create_oval(101, 22, 221, 134, outline="#e0d2f2", fill="#eee5f9")
        art.create_oval(116, 34, 206, 122, outline="#e2d5f5", dash=(2, 4))
        art.create_arc(*orbit, start=0, extent=180, style="arc", outline="#d4c4ec", dash=(3, 4))
        art.create_oval(119, 118, 203, 138, outline="#dfd2f2", fill="#e8def7")
        for x, y, r in ((108, 43, 3), (210, 50, 4), (112, 112, 2), (208, 105, 3)):
            art.create_polygon(x - r, y, x, y - r, x + r, y, x, y + r, fill="#a68acc")
        if self.frames:
            self.portrait = art.create_image(center_x, center_y, image=self.frames[0])
        else:
            self.portrait = art.create_text(center_x, center_y, text="化", font=(FONT, 28, "bold"), fill=INK)
        art.create_arc(*orbit, start=180, extent=180, style="arc", outline="#c6b4e5", dash=(3, 4))

        self._motif_positions = []
        for icon, fill, ink in MOTIFS:
            tag = f"motif-{icon}"
            shape = art.create_polygon(-28, -9, -28, -17, -19, -17, 19, -17, 28, -17, 28, -9,
                                       28, 9, 28, 17, 19, 17, -19, 17, -28, 17, -28, 9,
                                       smooth=True, fill=fill, outline=ink, tags=(tag,))
            self._draw_icon(icon, 0, 0, ink)
            icon_items = range(shape + 1, art.find_all()[-1] + 1)
            left, top, right, bottom = art.bbox(*icon_items)
            # Center the drawn bounds before grouping the icon with its orbiting card.
            for item in icon_items:
                art.move(item, -(left + right) / 2, -(top + bottom) / 2)
                art.addtag_withtag(tag, item)
            self._motif_positions.append((tag, (0, 0)))
        self._render_scene(0)

    def _draw_icon(self, icon, x, y, ink):
        art = self.art

        def coords(points):
            return [(x + a, y + b) for a, b in points]

        def line(*points, smooth=False):
            art.create_line(*coords(points), fill=ink, width=1.4, smooth=smooth)

        def dot(a, b, r, fill):
            art.create_oval(x + a - r, y + b - r, x + a + r, y + b + r, fill=fill, outline=ink)

        if icon == "benzene":
            art.create_polygon(*coords(((12, 2), (20, 7), (20, 16), (12, 21), (4, 16), (4, 7))),
                               fill="", outline=ink, width=1.4)
            dot(12, 11.5, 4.7, "")
        elif icon == "cell":
            art.create_polygon(*coords(((3, 12), (21, 12), (19.5, 20), (4.5, 20))), fill="#bde2dc", outline="")
            line((2, 8), (4, 20), (20, 20), (22, 8))
            line((8, 15), (8, 4), (16, 4), (16, 15))
            line((6, 16), (6, 9))
            line((18, 16), (18, 9))
        elif icon == "chromatogram":
            line((3, 3), (3, 20), (23, 20))
            line((4, 17), (7, 17), (9, 5), (11, 17), (15, 17), (18, 9), (20, 17), (23, 17))
        elif icon == "orbital":
            art.create_oval(x + 3, y + 5, x + 12, y + 17, fill="#bba6ec", outline=ink)
            art.create_oval(x + 12, y + 5, x + 21, y + 17, fill="#dfd3fa", outline=ink)
            dot(12, 11, 1.35, CARD)
        elif icon == "coordination":
            nodes = ((12, 2), (12, 21), (3, 11), (21, 11), (6, 17), (18, 5))
            for a, b in nodes:
                line((12, 11), (a, b))
            for a, b in nodes:
                dot(a, b, 1.7, CARD)
            dot(12, 11, 2.5, "#a9c7ef")
        elif icon == "energy":
            line((3, 3), (3, 20), (23, 20))
            line((5, 13), (8, 13), (10, 5), (13, 4), (16, 16), (20, 16), (23, 16), smooth=True)
