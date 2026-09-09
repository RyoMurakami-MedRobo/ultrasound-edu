"""Interactive DAS beamforming explorer (Tkinter + Matplotlib).

Python port of ``main_gui.m``. MATLAB is the source of truth; keep the two in
sync (see ``CLAUDE.md``). The MATLAB tool uses ``uifigure`` / ``uitabgroup`` /
``uitable``; this port uses the structural equivalents ``tk.Tk`` /
``ttk.Notebook`` / ``ttk.Treeview`` so no heavy GUI dependency is added to a
teaching tool.

Scope: MUST is MATLAB-only, so the "MUST das()" algorithm option is absent
here (see the module docstring of :mod:`ultrasound_das.pipeline`).

Run with ``python -m ultrasound_das`` or ``ultrasound-das-gui``.
"""

from __future__ import annotations

import importlib
import time
import tkinter as tk
from tkinter import messagebox, ttk

import numpy as np

try:
    import matplotlib

    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.colors import ListedColormap
    from matplotlib.figure import Figure
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "The GUI needs matplotlib with the TkAgg backend. Install it with "
        "`pip install matplotlib`.\n" + str(exc)
    )

from .das_reference import das_reference
from .metrics import compare_images, to_db
from .pipeline import ALGORITHMS, run_algorithm, run_all_tx
from .sim_engine import SimConfig, sim_engine
from .wave_animator import (
    WavePropagation,
    alignment_bundles,
    blue_white_red,
    delay_curve_image,
)

SWEEP_SECONDS = 20.0  # wall-clock duration of one full animation sweep at 1x

_PRESETS = ("Single point (PSF)", "Wire phantom", "Anechoic cyst", "Custom")


def _preset_data(name: str) -> np.ndarray:
    if name == "Single point (PSF)":
        return np.array([[0.0, 20.0, 1.0]])
    if name == "Wire phantom":
        XW, ZW = np.meshgrid(np.arange(-8, 9, 4), np.arange(8, 39, 6))
        return np.column_stack([XW.ravel(order="F"), ZW.ravel(order="F"), np.ones(XW.size)])
    if name == "Anechoic cyst":
        rng = np.random.default_rng(0)
        n = 500
        xb = (rng.random(n) * 2 - 1) * 12
        zb = 5 + rng.random(n) * 33
        rb = rng.standard_normal(n) * 0.35
        rb[np.hypot(xb - 0, zb - 22) < 5] = 0.0
        D = np.column_stack([xb, zb, rb])
        return np.vstack([D, [-8, 12, 1], [8, 32, 1]])
    return np.array([[0.0, 20.0, 1.0]])


class _Axes:
    """A matplotlib Figure embedded in a Tk frame."""

    def __init__(self, parent, n_axes=1, layout=(1, 1)):
        self.fig = Figure(figsize=(5, 4), constrained_layout=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.ax = [self.fig.add_subplot(layout[0], layout[1], i + 1) for i in range(n_axes)]

    def draw(self):
        self.canvas.draw_idle()


class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("Ultrasound DAS Beamforming Explorer (Python)")
        root.geometry("1500x920")

        self.S = None
        self.imgA = self.imgB = None
        self.nameA = self.nameB = ""
        self.msA = self.msB = float("nan")
        self.gx = self.gz = None
        self.pt = np.array([0.0, 20e-3])
        self.playing = False
        self.anim_t = 0.0
        self.wave: WavePropagation | None = None
        self._bwr = ListedColormap(blue_white_red())

        self.v = {}
        self._build()
        self._on_scheme_changed()
        self._apply_mode()
        self._apply_preset()
        self._draw_phantom()
        self._set_status("Ready. Press [1] Simulate, then [2] Beamform / Compare.  Backend: mock (MUST is MATLAB-only)")

    # ------------------------------------------------------------------
    #  construction
    # ------------------------------------------------------------------
    def _build(self):
        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True)

        self.ctrl = ttk.Frame(outer, width=380)
        self.ctrl.pack(side="left", fill="y")
        self.ctrl.pack_propagate(False)
        cc = tk.Canvas(self.ctrl, borderwidth=0, highlightthickness=0, width=360)
        sb = ttk.Scrollbar(self.ctrl, orient="vertical", command=cc.yview)
        self.cframe = ttk.Frame(cc)
        self.cframe.bind("<Configure>", lambda e: cc.configure(scrollregion=cc.bbox("all")))
        cc.create_window((0, 0), window=self.cframe, anchor="nw")
        cc.configure(yscrollcommand=sb.set)
        cc.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

        self.nb = ttk.Notebook(outer)
        self.nb.pack(side="right", fill="both", expand=True)

        self._build_controls(self.cframe)
        self._build_tab_image()
        self._build_tab_compare()
        self._build_tab_anim()
        self._build_tab_delay()

        self.status = ttk.Label(self.root, text="", relief="sunken", anchor="w")
        self.status.pack(side="bottom", fill="x")

    def _row(self, parent, label, widget):
        fr = ttk.Frame(parent)
        fr.pack(fill="x", pady=1)
        ttk.Label(fr, text=label, width=22, anchor="e").pack(side="left")
        widget.pack(side="left", fill="x", expand=True)
        return widget

    def _num(self, parent, label, key, default, width=10):
        var = tk.StringVar(value=str(default))
        self.v[key] = var
        return self._row(parent, label, ttk.Entry(parent, textvariable=var, width=width))

    def _combo(self, parent, label, key, values, default, cb=None):
        var = tk.StringVar(value=default)
        self.v[key] = var
        w = ttk.Combobox(parent, textvariable=var, values=list(values), state="readonly")
        if cb:
            w.bind("<<ComboboxSelected>>", lambda e: cb())
        return self._row(parent, label, w)

    def _check(self, parent, label, key, default, cb=None):
        var = tk.BooleanVar(value=default)
        self.v[key] = var
        w = ttk.Checkbutton(parent, variable=var, command=cb or (lambda: None))
        return self._row(parent, label, w)

    _SECTION_PACK = dict(fill="x", padx=4, pady=3)

    def _section(self, parent, title):
        lf = ttk.LabelFrame(parent, text=title)
        lf.pack(**self._SECTION_PACK)
        self._sections_in_order.append(lf)
        return lf

    def _show_section(self, sec):
        """(Re)pack ``sec`` in its original position, not appended at the end."""
        anchor = self._section_anchor.get(sec)
        opts = dict(self._SECTION_PACK)
        if anchor is not None and anchor.winfo_manager() == "pack":
            opts["before"] = anchor
        sec.pack(**opts)

    def _build_controls(self, p):
        self._sections_in_order = []
        self._section_anchor = {}
        hdr = self._section(p, "Control panel mode")
        self._combo(hdr, "Mode", "mode", ("Simple", "Detailed"), "Simple", self._apply_mode)

        s = self._section(p, "Transducer")
        self._combo(s, "Elements", "nel", ("16", "32", "64", "128", "192"), "64", self._draw_phantom)
        self._num(s, "Centre frequency [MHz]", "fc", 5)

        self.adv_probe = self._section(p, "Transducer (advanced)")
        self._num(self.adv_probe, "Pitch [mm]", "pitch", 0.30)
        self._num(self.adv_probe, "Bandwidth -6 dB [%]", "bw", 75)
        self._num(self.adv_probe, "Speed of sound [m/s]", "c", 1540)
        self._combo(self.adv_probe, "Sampling fs / fc", "fsfac", ("4", "6", "8", "12"), "4")

        s = self._section(p, "Transmit scheme")
        self._combo(s, "Transmit mode", "scheme",
                    ("Plane wave", "Focused", "Diverging"), "Plane wave", self._on_scheme_changed)
        self.w_angle = self._num(s, "Steering angle [deg]", "angle", 0)
        self.w_focus = self._num(s, "Focal depth F [mm]", "focus", "20")

        self.adv_tx = self._section(p, "Transmit (advanced)")
        self.w_srcdepth = self._num(self.adv_tx, "Diverging source [mm]", "srcdepth", 10)
        self._check(self.adv_tx, "Single-element transmit", "singleel", False, self._on_scheme_changed)
        self.w_elidx = self._num(self.adv_tx, "Element index", "elidx", 32)

        s = self._section(p, "Targets (phantom)")
        self._combo(s, "Preset", "preset", _PRESETS, "Single point (PSF)")
        brow = ttk.Frame(s)
        brow.pack(fill="x", pady=1)
        ttk.Button(brow, text="Apply preset", command=self._apply_preset_refresh).pack(side="left")
        self.v["clickadd"] = tk.BooleanVar(value=True)
        ttk.Checkbutton(brow, text="Left-click to add", variable=self.v["clickadd"]).pack(side="left")
        cols = ("x", "z", "amp")
        self.tbl = ttk.Treeview(s, columns=cols, show="headings", height=7, selectmode="extended")
        for c, t in zip(cols, ("X [mm]", "Z [mm]", "Amplitude")):
            self.tbl.heading(c, text=t)
            self.tbl.column(c, width=90, anchor="center")
        self.tbl.pack(fill="x", pady=2)
        self.tbl.bind("<Double-1>", self._edit_cell)
        brow2 = ttk.Frame(s)
        brow2.pack(fill="x")
        ttk.Button(brow2, text="Add row", command=self._add_row).pack(side="left")
        ttk.Button(brow2, text="Delete selected", command=self._del_rows).pack(side="left")
        ttk.Button(brow2, text="Clear all", command=self._clear_rows).pack(side="left")

        s = self._section(p, "Reconstruction")
        self._num(s, "Receive f-number (0=full)", "fnum", 1.5)
        self._num(s, "Dynamic range [dB]", "dr", 50)

        self.adv_recon = self._section(p, "Reconstruction grid (advanced)")
        self._num(self.adv_recon, "Pixels X", "nx", 161)
        self._num(self.adv_recon, "Pixels Z", "nz", 221)
        self._num(self.adv_recon, "Lateral half-width [mm]", "xhalf", 12)
        self._num(self.adv_recon, "Depth min [mm]", "zmin", 3)
        self._num(self.adv_recon, "Depth max [mm]", "zmax", 40)
        self._combo(self.adv_recon, "Receive apodisation", "rxapod", ("rect", "hann"), "rect")

        s = self._section(p, "DAS algorithm comparison")
        algs = tuple(ALGORITHMS)
        self._combo(s, "Algorithm A", "algA", algs, "Reference DAS")
        self._combo(s, "Algorithm B", "algB", ("None",) + algs, "Custom DAS")

        self.adv_algo = self._section(p, "Algorithm (advanced)")
        self._num(self.adv_algo, "Custom module:function", "customfn",
                  "ultrasound_das.das_custom_template:das_custom_template", width=32)

        s = self._section(p, "Run")
        ttk.Button(s, text="[1] Simulate (generate RF)", command=self._on_simulate).pack(fill="x", pady=2)
        ttk.Button(s, text="[2] Beamform / Compare", command=self._on_beamform).pack(fill="x", pady=2)

        self._adv_sections = [self.adv_probe, self.adv_tx, self.adv_recon, self.adv_algo]
        # Anchor = the next section in build order, so a hidden advanced panel
        # re-appears in its original slot rather than at the bottom.
        order = self._sections_in_order
        for sec in self._adv_sections:
            i = order.index(sec)
            self._section_anchor[sec] = order[i + 1] if i + 1 < len(order) else None

    def _build_tab_image(self):
        t = ttk.Frame(self.nb)
        self.nb.add(t, text="Phantom & B-mode")
        self.ax_img = _Axes(t, 3, (1, 3))
        self.ax_img.fig.canvas.mpl_connect("button_press_event", self._on_canvas_click)

    def _build_tab_compare(self):
        t = ttk.Frame(self.nb)
        self.nb.add(t, text="Compare & metrics")
        top = ttk.Frame(t)
        top.pack(fill="both", expand=True)
        self.ax_cmp = _Axes(top, 2, (1, 2))
        self.met = ttk.Treeview(t, columns=("m", "a", "b"), show="headings", height=8)
        for c, w, txt in (("m", 320, "Metric"), ("a", 180, "A"), ("b", 180, "B")):
            self.met.heading(c, text=txt)
            self.met.column(c, width=w, anchor="w")
        self.met.pack(fill="x")

    def _build_tab_anim(self):
        t = ttk.Frame(self.nb)
        self.nb.add(t, text="Wave animation")
        self.ax_anim = _Axes(t, 1)
        bar = ttk.Frame(t)
        bar.pack(fill="x")
        self.btn_play = ttk.Button(bar, text="Play", command=self._play_toggle)
        self.btn_play.pack(side="left")
        self.tslider = ttk.Scale(bar, from_=0, to=1, orient="horizontal", command=self._scrub)
        self.tslider.pack(side="left", fill="x", expand=True, padx=6)
        self.lbl_t = ttk.Label(bar, text="t = 0.00 us", width=14)
        self.lbl_t.pack(side="left")
        self.v["speed"] = tk.StringVar(value="1x")
        ttk.Combobox(bar, textvariable=self.v["speed"], width=6, state="readonly",
                     values=("0.25x", "0.5x", "1x", "2x", "4x")).pack(side="left")
        self.v["txsel"] = tk.StringVar(value="Transmit 1")
        self.txsel = ttk.Combobox(bar, textvariable=self.v["txsel"], width=12, state="readonly",
                                  values=("Transmit 1",))
        self.txsel.bind("<<ComboboxSelected>>", lambda e: (self._setup_animation(), self._update_delay_views()))
        self.txsel.pack(side="left")

    def _build_tab_delay(self):
        t = ttk.Frame(self.nb)
        self.nb.add(t, text="Delay curve & alignment")
        self.ax_dl = _Axes(t, 3, (1, 3))
        bar = ttk.Frame(t)
        bar.pack(fill="x")
        ttk.Label(bar, text="Reconstruction point (or click B-mode A):  X [mm]").pack(side="left")
        self.v["ptx"] = tk.StringVar(value="0")
        e1 = ttk.Entry(bar, textvariable=self.v["ptx"], width=7)
        e1.pack(side="left")
        e1.bind("<Return>", lambda e: self._on_point_edited())
        ttk.Label(bar, text="Z [mm]").pack(side="left")
        self.v["ptz"] = tk.StringVar(value="20")
        e2 = ttk.Entry(bar, textvariable=self.v["ptz"], width=7)
        e2.pack(side="left")
        e2.bind("<Return>", lambda e: self._on_point_edited())
        self.lbl_delay = ttk.Label(bar, text="")
        self.lbl_delay.pack(side="left", padx=8)

    # ------------------------------------------------------------------
    #  events / state
    # ------------------------------------------------------------------
    def _apply_mode(self):
        simple = self.v["mode"].get() == "Simple"
        for sec in self._adv_sections:
            if simple:
                sec.pack_forget()
            else:
                self._show_section(sec)

    def _on_scheme_changed(self):
        """Enable/disable transmit controls by scheme + single-element, like
        ``onSchemeChanged`` in main_gui.m."""
        k = self._scheme_key()
        single = bool(self.v["singleel"].get())

        def _en(w, on):
            try:
                w.configure(state="normal" if on else "disabled")
            except tk.TclError:
                pass

        _en(self.w_focus, k == "focused" and not single)
        _en(self.w_srcdepth, k == "diverging" and not single)
        _en(self.w_elidx, single)
        _en(self.w_angle, not single)
        self._draw_phantom()

    def _scheme_key(self):
        return {"Focused": "focused", "Diverging": "diverging"}.get(self.v["scheme"].get(), "plane")

    def _f(self, key, default=0.0):
        try:
            return float(self.v[key].get())
        except (TypeError, ValueError):
            return default

    def _parse_list(self, s):
        out = []
        for tok in s.replace(";", ",").replace(" ", ",").split(","):
            if tok:
                try:
                    x = float(tok)
                    if x != 0:
                        out.append(x)
                except ValueError:
                    pass
        return out

    def _table_data(self):
        rows = []
        for iid in self.tbl.get_children():
            rows.append([float(x) for x in self.tbl.item(iid, "values")])
        return np.array(rows).reshape(-1, 3) if rows else np.zeros((0, 3))

    def _set_table(self, D):
        self.tbl.delete(*self.tbl.get_children())
        for r in np.atleast_2d(D):
            self.tbl.insert("", "end", values=(f"{r[0]:.6g}", f"{r[1]:.6g}", f"{r[2]:.6g}"))

    def _edit_cell(self, event):
        iid = self.tbl.identify_row(event.y)
        col = self.tbl.identify_column(event.x)
        if not iid or not col:
            return
        ci = int(col[1:]) - 1
        x0, y0, w, h = self.tbl.bbox(iid, col)
        cur = self.tbl.item(iid, "values")[ci]
        ed = ttk.Entry(self.tbl)
        ed.insert(0, cur)
        ed.select_range(0, "end")
        ed.focus_set()
        ed.place(x=x0, y=y0, width=w, height=h)

        def commit(_=None):
            vals = list(self.tbl.item(iid, "values"))
            try:
                vals[ci] = f"{float(ed.get()):.6g}"
            except ValueError:
                pass
            self.tbl.item(iid, values=vals)
            ed.destroy()
            self.v["preset"].set("Custom")
            self._draw_phantom()

        ed.bind("<Return>", commit)
        ed.bind("<FocusOut>", commit)

    def _add_row(self):
        self.tbl.insert("", "end", values=("0", "20", "1"))
        self.v["preset"].set("Custom")
        self._draw_phantom()

    def _del_rows(self):
        for iid in self.tbl.selection():
            self.tbl.delete(iid)
        self.v["preset"].set("Custom")
        self._draw_phantom()

    def _clear_rows(self):
        self.tbl.delete(*self.tbl.get_children())
        self.v["preset"].set("Custom")
        self._draw_phantom()

    def _apply_preset(self):
        name = self.v["preset"].get()
        if name != "Custom":
            self._set_table(_preset_data(name))

    def _apply_preset_refresh(self):
        self._apply_preset()
        self._draw_phantom()
        self._set_status("Preset applied. Press [1] Simulate.")

    def _delete_nearest_scatterer(self, x, z):
        D = self._table_data()
        if D.shape[0] == 0:
            return
        k = int(np.argmin(np.hypot(D[:, 0] - x, D[:, 1] - z)))
        if np.hypot(D[k, 0] - x, D[k, 1] - z) > 3.0:      # 3 mm radius, like main_gui.m
            self._set_status("No scatterer within 3 mm of the cursor.")
            return
        self.tbl.delete(self.tbl.get_children()[k])
        self.v["preset"].set("Custom")
        self._draw_phantom()
        self._set_status("Nearest scatterer deleted.")

    def _on_canvas_click(self, event):
        if event.inaxes is None:
            return
        ax_list = self.ax_img.ax
        if event.inaxes is ax_list[0]:          # phantom axes
            x, z = event.xdata, event.ydata
            if x is None or z is None:
                return
            if event.button == 3:              # right-click: delete nearest
                self._delete_nearest_scatterer(x, z)
                return
            if not self.v["clickadd"].get() or z <= 0:
                return
            self.tbl.insert("", "end", values=(f"{x:.6g}", f"{z:.6g}", "1"))
            self.v["preset"].set("Custom")
            self._draw_phantom()
            self._set_status(f"Scatterer added at ({x:.2f}, {z:.2f}) mm")
        elif event.inaxes is ax_list[1] and self.imgA is not None:   # B-mode A: move point
            self.pt = np.array([event.xdata * 1e-3, event.ydata * 1e-3])
            self.v["ptx"].set(f"{event.xdata:.2f}")
            self.v["ptz"].set(f"{event.ydata:.2f}")
            self._refresh_images()
            self._update_delay_views()

    def _on_point_edited(self):
        self.pt = np.array([self._f("ptx") * 1e-3, self._f("ptz") * 1e-3])
        self._refresh_images()
        self._update_delay_views()

    # ------------------------------------------------------------------
    #  config
    # ------------------------------------------------------------------
    def _read_cfg(self) -> SimConfig:
        try:
            nel = int(round(float(self.v["nel"].get())))
        except ValueError:
            nel = 64
            self.v["nel"].set("64")
        focus = self._parse_list(self.v["focus"].get()) or [20.0]
        D = self._table_data()
        return SimConfig(
            n_elements=nel,
            pitch=self._f("pitch", 0.30) * 1e-3,
            fc=self._f("fc", 5) * 1e6,
            bandwidth=self._f("bw", 75),
            c=self._f("c", 1540),
            fs_factor=float(self.v["fsfac"].get()),
            scheme=self._scheme_key(),
            angle_deg=self._f("angle", 0),
            focus_mm=focus,
            src_mm=self._f("srcdepth", 10),
            single_element=bool(self.v["singleel"].get()),
            element_index=int(round(self._f("elidx", 32))),
            scat_x=D[:, 0] * 1e-3,
            scat_z=D[:, 1] * 1e-3,
            scat_rc=D[:, 2],
            zmax=self._f("zmax", 40) * 1e-3,
            fnumber=self._f("fnum", 1.5),
            rx_apod=self.v["rxapod"].get(),
        )

    def _recon_grid(self):
        xh = self._f("xhalf", 12) * 1e-3
        gx = np.linspace(-xh, xh, int(round(self._f("nx", 161))))
        z0 = min(self._f("zmin", 3), self._f("zmax", 40) - 1) * 1e-3
        gz = np.linspace(max(z0, 1e-4), self._f("zmax", 40) * 1e-3, int(round(self._f("nz", 221))))
        return gx, gz

    # ------------------------------------------------------------------
    #  run
    # ------------------------------------------------------------------
    def _on_simulate(self):
        cfg = self._read_cfg()
        self._set_status("Simulating...")
        self.root.update_idletasks()
        try:
            self.S = sim_engine(cfg)
        except Exception as exc:
            messagebox.showerror("Simulation failed", str(exc))
            self._set_status(f"Simulation failed: {exc}")
            return
        self.imgA = self.imgB = None
        n = len(self.S.tx)
        self.txsel.configure(values=[f"Transmit {i+1}" for i in range(n)])
        self.v["txsel"].set("Transmit 1")
        self._draw_phantom()
        self._setup_animation()
        self._update_delay_views()
        self._set_status(
            f"RF generated | backend = {self.S.backend} | {n} transmit event(s) | "
            f"RF {self.S.RF.shape[0]} x {self.S.RF.shape[1]} | fs = {self.S.fs/1e6:.1f} MHz | "
            f"{self.S.elapsed:.2f} s"
        )

    def _current_tx(self):
        try:
            return max(0, int(self.v["txsel"].get().split()[-1]) - 1)
        except (ValueError, IndexError):
            return 0

    def _algo_callable(self, key):
        if key in ALGORITHMS:
            return ALGORITHMS[key], key
        spec = self.v["customfn"].get().strip()
        mod, _, fn = spec.partition(":")
        if not fn:
            raise ValueError("Custom function must be 'module:function'.")
        m = importlib.import_module(mod)
        return getattr(m, fn), fn

    def _on_beamform(self):
        if self.S is None:
            self._on_simulate()
            if self.S is None:
                return
        gx, gz = self._recon_grid()
        self.gx, self.gz = gx, gz
        for tx in self.S.tx:
            tx.fnumber = self._f("fnum", 1.5)
            tx.rx_apod = self.v["rxapod"].get()

        self._set_status("Beamforming...")
        self.root.update_idletasks()
        try:
            fhA, self.nameA = self._algo_callable(self.v["algA"].get())
            self.imgA, self.msA, self.nameA = run_algorithm(fhA, self.S, gx, gz)
        except Exception as exc:
            messagebox.showerror("Algorithm A failed", str(exc))
            self._set_status(f"A failed: {exc}")
            return

        self.imgB, self.msB, self.nameB = None, float("nan"), ""
        if self.v["algB"].get() != "None":
            try:
                fhB, self.nameB = self._algo_callable(self.v["algB"].get())
                self.imgB, self.msB, self.nameB = run_algorithm(fhB, self.S, gx, gz)
            except Exception as exc:
                messagebox.showerror("Algorithm B failed", str(exc))
                self.imgB = None

        self._refresh_images()
        self._update_compare()
        self._update_delay_views()
        msg = f"A: {self.nameA} = {self.msA:.1f} ms"
        if self.imgB is not None:
            msg += f"   |   B: {self.nameB} = {self.msB:.1f} ms"
        self._set_status("Beamforming done.   " + msg)

    # ------------------------------------------------------------------
    #  drawing
    # ------------------------------------------------------------------
    def _draw_phantom(self):
        ax = self.ax_img.ax[0]
        ax.clear()
        D = self._table_data()
        try:
            nel = int(round(float(self.v["nel"].get())))
        except ValueError:
            nel = 64
        ex = (np.arange(nel) - (nel - 1) / 2) * self._f("pitch", 0.30)
        ax.plot(ex, np.zeros(nel), "s", ms=3, color="0.7", mec="0.3")
        if D.size:
            pos = D[:, 2] >= 0
            ax.scatter(D[pos, 0], D[pos, 1], s=30, marker="o", edgecolors="g", facecolors="#88e088")
            ax.scatter(D[~pos, 0], D[~pos, 1], s=30, marker="o", edgecolors="r", facecolors="#ffb0b0")
        xh = max(self._f("xhalf", 12), np.max(np.abs(ex)) * 1.1 if nel else 12)
        ax.set_xlim(-xh, xh)
        ax.set_ylim(self._f("zmax", 40), -2)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("z [mm]")
        ax.set_title(f"Phantom layout ({D.shape[0]} scatterers)")
        ax.grid(True)
        self.ax_img.draw()

    def _show_bmode(self, ax, img, name):
        ax.clear()
        if img is None:
            ax.set_title("(not computed)")
            ax.axis("off")
            return
        ax.axis("on")
        dr = self._f("dr", 50)
        L = 20 * np.log10(img / np.max(img) + 1e-12)
        ax.imshow(L, extent=[self.gx[0] * 1e3, self.gx[-1] * 1e3, self.gz[-1] * 1e3, self.gz[0] * 1e3],
                  cmap="gray", vmin=-dr, vmax=0, aspect="equal")
        D = self._table_data()
        if D.size and D.shape[0] <= 60:
            ax.plot(D[:, 0], D[:, 1], "o", ms=6, mec="#33ff33", mfc="none")
        ax.plot(self.pt[0] * 1e3, self.pt[1] * 1e3, "+", color="#ff6600", ms=12, mew=1.5)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("z [mm]")
        ax.set_title(f"{name}  (DR {round(dr)} dB)")

    def _refresh_images(self):
        self._show_bmode(self.ax_img.ax[1], self.imgA, self.nameA or "Algorithm A")
        self._show_bmode(self.ax_img.ax[2], self.imgB, self.nameB or "Algorithm B")
        self.ax_img.draw()

    def _update_compare(self):
        axd, axp = self.ax_cmp.ax
        axd.clear()
        axp.clear()
        self.met.delete(*self.met.get_children())
        if self.imgA is None:
            self.ax_cmp.draw()
            return
        cmp = compare_images(self.imgA, self.imgB, self.gx, self.gz)

        if cmp["diff_db"] is not None:
            axd.imshow(cmp["diff_db"],
                       extent=[self.gx[0] * 1e3, self.gx[-1] * 1e3, self.gz[-1] * 1e3, self.gz[0] * 1e3],
                       cmap="viridis", vmin=-80, vmax=0, aspect="equal")
            axd.set_title(f"|A - B| [dB]   MSE = {cmp['mse']:.3e}")
            axd.set_xlabel("x [mm]")
            axd.set_ylabel("z [mm]")
        else:
            axd.text(0.5, 0.5, "Select algorithm B to see the difference", ha="center")
            axd.axis("off")

        axp.plot(self.gx * 1e3, to_db(cmp["prof_a"]), "-", lw=1.6, color="#0059cc", label=f"A: {self.nameA}")
        if cmp["prof_b"] is not None:
            axp.plot(self.gx * 1e3, to_db(cmp["prof_b"]), "--", lw=1.6, color="#d9331a", label=f"B: {self.nameB}")
        axp.axhline(-6, ls=":", color="0.4", label="-6 dB (FWHM)")
        axp.set_ylim(-max(self._f("dr", 50), 40), 2)
        axp.set_xlabel("x [mm]")
        axp.set_ylabel("Normalised amplitude [dB]")
        axp.set_title(f"Lateral profile at z = {cmp['peak_a_mm'][1]:.2f} mm")
        axp.legend(loc="lower left", fontsize=8)
        axp.grid(True)
        self.ax_cmp.draw()

        def dash(v, fmt="{:.1f}"):
            return "-" if v is None or (isinstance(v, float) and np.isnan(v)) else fmt.format(v)

        rows = [
            ("Run time [ms]", dash(self.msA), dash(self.msB)),
            ("Peak position x [mm]", f"{cmp['peak_a_mm'][0]:.3f}", dash(cmp["peak_b_mm"][0], "{:.3f}")),
            ("Peak position z [mm]", f"{cmp['peak_a_mm'][1]:.3f}", dash(cmp["peak_b_mm"][1], "{:.3f}")),
            ("Lateral FWHM (-6 dB) [mm]", dash(cmp["fwhm_a_mm"], "{:.4f}"), dash(cmp["fwhm_b_mm"], "{:.4f}")),
            ("MSE (normalised by max(A))", "-", dash(cmp["mse"], "{:.4e}")),
            ("Max absolute error (vs A)", "-", dash(cmp["max_abs_err"], "{:.4e}")),
            ("Algorithm", self.nameA, self.nameB or "-"),
            ("Backend", self.S.backend, ""),
        ]
        for r in rows:
            self.met.insert("", "end", values=r)

    # ------------------------------------------------------------------
    #  wave animation
    # ------------------------------------------------------------------
    def _setup_animation(self):
        if self.S is None:
            return
        xh = max(self._f("xhalf", 12) * 1e-3, np.max(np.abs(self.S.rx_pos)) * 1.05)
        self.wave = WavePropagation.setup(
            self.S, self._current_tx(), (-xh, xh), (0.0, self._f("zmax", 40) * 1e-3)
        )
        self.tslider.configure(to=self.wave.tmax)
        self.tslider.set(0)
        self.anim_t = 0.0
        self._draw_wave(0.0)

    def _draw_wave(self, t):
        if self.wave is None:
            return
        ax = self.ax_anim.ax[0]
        ax.clear()
        w = self.wave
        field, hit_x = w.frame(t)
        ax.imshow(field, extent=[w.xv[0] * 1e3, w.xv[-1] * 1e3, w.zv[-1] * 1e3, w.zv[0] * 1e3],
                  cmap=self._bwr, vmin=-1, vmax=1, aspect="equal")
        ax.plot(w.elem_x * 1e3, w.elem_z * 1e3, "s", ms=3, color="0.75", mec="0.25")
        if w.act_elem_x.size:
            ax.plot(w.act_elem_x * 1e3, np.zeros(w.act_elem_x.size), "s", ms=3, color="#ff6666", mec="#990000")
        if hit_x.size:
            ax.plot(hit_x * 1e3, np.zeros(hit_x.size), "s", ms=6, color="#4d99ff", mec="#000099")
        if w.xs.size:
            ax.plot(w.xs * 1e3, w.zs * 1e3, "o", ms=6, mec="g", mfc="none")
        if np.all(np.isfinite(w.src_xz)) and w.scheme != "plane":
            ax.plot(w.src_xz[0] * 1e3, w.src_xz[1] * 1e3, "*", ms=12, mec="#804d00", mfc="#ffcc33")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("z [mm]")
        ax.set_title(f"t = {t*1e6:.2f} us   (red: transmit wavefront / blue: scattered echo)")
        self.ax_anim.draw()
        self.lbl_t.configure(text=f"t = {t*1e6:.2f} us")

    def _scrub(self, val):
        if self.S is None or self.playing:
            return
        self.anim_t = float(val)
        self._draw_wave(self.anim_t)

    def _play_toggle(self):
        if self.playing:
            self.playing = False
            return
        if self.wave is None:
            self._set_status("Press [1] Simulate first.")
            return
        self.playing = True
        self.btn_play.configure(text="Stop")
        spd = float(self.v["speed"].get().replace("x", ""))
        tmax = self.wave.tmax
        dur = SWEEP_SECONDS / spd
        t_start = self.anim_t
        clk = time.perf_counter()

        def step():
            if not self.playing:
                self.btn_play.configure(text="Play")
                return
            # wall-clock pacing (mirrors main_gui.m:988): simulated time is
            # derived from elapsed real time, never from a frame count.
            t = (t_start + tmax * (time.perf_counter() - clk) / dur) % tmax
            self.anim_t = t
            self.tslider.set(t)
            self._draw_wave(t)
            self.root.after(20, step)

        step()

    # ------------------------------------------------------------------
    #  delay curve / alignment
    # ------------------------------------------------------------------
    def _update_delay_views(self):
        if self.S is None:
            return
        k = self._current_tx()
        name = "das_reference"
        try:
            fh, name = self._algo_callable(self.v["algA"].get())
        except Exception:
            fh = das_reference
        try:
            _, tau = fh(self.S.RF[:, :, k], self.S.tx[k], self.S.rx_pos,
                        np.array([self.pt[0]]), np.array([self.pt[1]]), self.S.c, self.S.fs,
                        want_delays=True)
        except Exception as exc:
            _, tau = das_reference(self.S.RF[:, :, k], self.S.tx[k], self.S.rx_pos,
                                   np.array([self.pt[0]]), np.array([self.pt[1]]),
                                   self.S.c, self.S.fs, want_delays=True)
            name = f"das_reference (fallback: {exc})"
        tau = tau[0, :]

        axrf, axpre, axpost = self.ax_dl.ax
        for a in self.ax_dl.ax:
            a.clear()

        disp_rf, t_us, tau_us, ok = delay_curve_image(self.S, k, tau)
        nel = disp_rf.shape[1]
        axrf.imshow(disp_rf, extent=[0.5, nel + 0.5, t_us[-1], t_us[0]], cmap="gray",
                    aspect="auto", vmin=-np.max(np.abs(disp_rf)) * 0.6, vmax=np.max(np.abs(disp_rf)) * 0.6)
        idx = np.where(ok)[0]
        axrf.plot(idx + 1, tau_us[ok], "-", color="#ff4019", lw=2)
        axrf.plot(idx + 1, tau_us[ok], ".", color="#ffcc00", ms=6)
        axrf.invert_yaxis()
        axrf.set_xlabel("Receive element index")
        axrf.set_ylabel("Time [us]")
        axrf.set_title(f"Delay curve  ({self.pt[0]*1e3:.2f}, {self.pt[1]*1e3:.2f}) mm  {ok.sum()}/{nel} el.")

        bundles = alignment_bundles(self.S, k, tau)
        if bundles is not None:
            self._draw_bundle(axpre, bundles, bundles["pre"], show_tau=True)
            axpre.set_title(f"Before alignment (window at t = {bundles['tc']*1e6:.2f} us)")
            self._draw_bundle(axpost, bundles, bundles["post"], show_tau=False)
            st = bundles["sum_trace"]
            ss = np.max(np.abs(st)) + np.finfo(float).eps
            nel_b = bundles["pre"].shape[1]
            ybase = -0.10 * nel_b - 2
            axpost.plot(bundles["trel_us"], ybase + 0.9 * (0.06 * nel_b + 2) * st / ss, "-",
                        color="#d91a1a", lw=1.8)
            axpost.set_title(f"After alignment + sum  ({self.pt[0]*1e3:.2f}, {self.pt[1]*1e3:.2f}) mm")
        self.ax_dl.draw()
        self.lbl_delay.configure(
            text=f"Delays from: {name}   {np.isfinite(tau).sum()} / {tau.size} elements used"
        )

    def _draw_bundle(self, ax, b, W, show_tau):
        nel = W.shape[1]
        g, sc = b["gain"], b["scale"]
        for e in range(nel):
            col = "#0059bf" if b["ok"][e] else "0.75"
            ax.plot(b["trel_us"], e + 1 + g * W[:, e] / sc, "-", color=col, lw=0.6)
        if show_tau:
            d = (np.asarray(b["tau"]) - b["tc"]) * 1e6
            idx = np.where(b["ok"])[0]
            ax.plot(d[b["ok"]], idx + 1, "-", color="#ff4d00", lw=1.2)
            ax.plot(d[b["ok"]], idx + 1, ".", color="#ff4d00", ms=8)
        else:
            ax.axvline(0, color="#ff4d00", lw=1.2)
        ax.set_xlabel("Relative time [us]")
        ax.set_ylabel("Receive element index")
        ax.grid(True)

    # ------------------------------------------------------------------
    def _set_status(self, msg):
        self.status.configure(text="  " + msg)


def main(argv=None):
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
