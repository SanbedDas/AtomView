"""The AtomView main window."""
from __future__ import annotations

import os
import threading
import time
import traceback

os.environ.setdefault("QT_API", "pyside6")

import numpy as np
from PySide6.QtCore import QObject, QSettings, Qt, QTimer, Signal
from PySide6.QtGui import QAction, QActionGroup, QKeySequence
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDockWidget,
                               QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
                               QListWidget, QListWidgetItem, QMainWindow, QMenu, QMessageBox, QProgressBar, QProgressDialog,
                               QPushButton, QScrollArea, QSlider, QSpinBox, QSplitter, QStyle, QToolBar,
                               QToolButton, QVBoxLayout, QWidget)
from pyvistaqt import QtInteractor

from .. import modifiers as _mods  # noqa: F401  (registers modifiers)
from ..data import scalar_columns
from ..io import load as load_trajectory
from ..io import write as write_frame
from ..pipeline import REGISTRY, Cancelled, Modifier, Pipeline, is_session
from ..render import (RenderSettings, apply_quality, build_scene, default_camera, fit_camera, render_animation, render_image)
from .inspector import Inspector
from .param_editor import ColorButton, ModifierEditor

DATA_FILTER = ("Particle files (*.dump *.lammpstrj *.lammps *.data *.lmp *.xyz *.extxyz *.gz *.mat POSCAR* CONTCAR* "
               "XDATCAR* *.vasp *.cif *.pdb *.gro *.xtc *.trr *.dcd *.nc *.tng *.tpr *.mol2 *.psf *.mmtf *.xsf *.cfg "
               "*.traj OUTCAR*);;All files (*)")
CATEGORY_ORDER = ["Coloring", "Selection", "Modification", "Analysis", "Visual", "Scripting"]


class _Signals(QObject):
    computed = Signal(object, int, int, float)  # frame, frame_index, revision, seconds
    failed = Signal(str)
    progress = Signal(str)


class MainWindow(QMainWindow):
    def __init__(self, path=None):
        super().__init__()
        self.setWindowTitle("AtomView")
        self.resize(1500, 950)
        self.settings_store = QSettings("AtomView", "AtomView")
        self.pipeline: Pipeline | None = None
        self.frame_index = 0
        self.current = None
        self.render_settings = RenderSettings()
        self._busy = False
        self._pending = False
        self._first_view = True
        self._signals = _Signals()
        self._signals.computed.connect(self._on_computed)
        self._signals.failed.connect(lambda msg: self.statusBar().showMessage(msg, 10000))
        self._signals.progress.connect(lambda msg: self.busy_label.setText(msg))

        # 3D view
        self.view = QtInteractor(self)
        self.view.set_background(self.render_settings.background)
        self.view.add_axes()
        apply_quality(self.view, self.render_settings)
        self.view.track_click_position(self._on_pick, side="right")
        center = QWidget()
        cl = QVBoxLayout(center)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.setSpacing(0)
        cl.addWidget(self.view.interactor, 1)
        cl.addWidget(self._make_timeline())
        self.setCentralWidget(center)

        self._make_pipeline_dock()
        self._make_render_dock()
        self.inspector = Inspector()
        dock = QDockWidget("Data inspector", self)
        dock.setObjectName("inspector")
        dock.setWidget(self.inspector)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)
        self.inspector_dock = dock
        self._make_menus()
        self._make_toolbar()
        self.resizeDocks([self.inspector_dock], [230], Qt.Vertical)
        self.resizeDocks([self.pipeline_dock, self.render_dock], [380, 300], Qt.Horizontal)
        state = self.settings_store.value("window_state")
        if state is not None:
            self.restoreState(state)

        self.recompute_timer = QTimer(self, singleShot=True, interval=120)
        self.recompute_timer.timeout.connect(self._start_compute)
        self.play_timer = QTimer(self)
        self.play_timer.timeout.connect(self._advance)
        self.busy_label = QLabel()
        self.busy_bar = QProgressBar()
        self.busy_bar.setRange(0, 0)
        self.busy_bar.setMaximumWidth(120)
        self.busy_cancel = QPushButton("Cancel")
        self.busy_cancel.clicked.connect(self.cancel_compute)
        for w in (self.busy_label, self.busy_bar, self.busy_cancel):
            self.statusBar().addPermanentWidget(w)
            w.hide()
        self.busy_timer = QTimer(self, singleShot=True, interval=500)
        self.busy_timer.timeout.connect(lambda: [w.show() for w in (self.busy_label, self.busy_bar, self.busy_cancel)] if self._busy else None)
        self.setAcceptDrops(True)
        self.statusBar().showMessage("Open a file (⌘O) or drop one onto the window. Right-click an atom to inspect it.")
        self.session_path = None
        self._saved_state = None
        self.title_timer = QTimer(self, interval=1000)
        self.title_timer.timeout.connect(self._update_title)
        self.title_timer.start()
        if path:
            QTimer.singleShot(0, lambda: self.open_any(path))

    # ------------------------------------------------------------------ layout
    def _make_timeline(self):
        bar = QWidget()
        h = QHBoxLayout(bar)
        h.setContentsMargins(8, 4, 8, 4)
        st = self.style()
        self.btn_first = QToolButton(icon=st.standardIcon(QStyle.SP_MediaSkipBackward))
        self.btn_prev = QToolButton(icon=st.standardIcon(QStyle.SP_MediaSeekBackward))
        self.btn_play = QToolButton(icon=st.standardIcon(QStyle.SP_MediaPlay))
        self.btn_next = QToolButton(icon=st.standardIcon(QStyle.SP_MediaSeekForward))
        self.btn_last = QToolButton(icon=st.standardIcon(QStyle.SP_MediaSkipForward))
        self.btn_first.clicked.connect(lambda: self.set_frame(0))
        self.btn_prev.clicked.connect(lambda: self.set_frame(self.frame_index - 1))
        self.btn_next.clicked.connect(lambda: self.set_frame(self.frame_index + 1))
        self.btn_last.clicked.connect(lambda: self.set_frame(self.n_frames - 1))
        self.btn_play.clicked.connect(self.toggle_play)
        for b in (self.btn_first, self.btn_prev, self.btn_play, self.btn_next, self.btn_last):
            h.addWidget(b)
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0, 0)
        self.slider.valueChanged.connect(self.set_frame)
        h.addWidget(self.slider, 1)
        self.frame_spin = QSpinBox()
        self.frame_spin.setRange(0, 0)
        self.frame_spin.valueChanged.connect(self.set_frame)
        h.addWidget(self.frame_spin)
        self.frame_label = QLabel("/ 0")
        h.addWidget(self.frame_label)
        h.addWidget(QLabel("  fps"))
        self.fps = QSpinBox()
        self.fps.setRange(1, 60)
        self.fps.setValue(10)
        h.addWidget(self.fps)
        self.loop = QCheckBox("loop")
        self.loop.setChecked(True)
        h.addWidget(self.loop)
        return bar

    def _make_pipeline_dock(self):
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(4, 4, 4, 4)
        row = QHBoxLayout()
        self.add_btn = QPushButton("＋ Add modifier")
        menu = QMenu(self)
        by_cat = {}
        for name, cls in REGISTRY.items():
            if name == "Function":
                continue
            by_cat.setdefault(cls.category, []).append((cls.title, name))
        for cat in CATEGORY_ORDER + sorted(set(by_cat) - set(CATEGORY_ORDER)):
            if cat not in by_cat:
                continue
            sub = menu.addMenu(cat)
            for title, name in sorted(by_cat[cat]):
                act = sub.addAction(title)
                act.setToolTip(REGISTRY[name].description)
                act.triggered.connect(lambda _=False, n=name: self.add_modifier(n))
        menu.setToolTipsVisible(True)
        self.add_btn.setMenu(menu)
        row.addWidget(self.add_btn, 1)
        for text, tip, fn in (("↑", "Move up", lambda: self.move_modifier(-1)),
                              ("↓", "Move down", lambda: self.move_modifier(1)),
                              ("⧉", "Duplicate", self.duplicate_modifier),
                              ("✕", "Delete", self.remove_modifier)):
            b = QToolButton(text=text)
            b.setToolTip(tip)
            b.clicked.connect(fn)
            row.addWidget(b)
        v.addLayout(row)
        self.list = QListWidget()
        self.list.currentRowChanged.connect(self._show_editor)
        self.list.itemChanged.connect(self._item_toggled)
        split = QSplitter(Qt.Vertical)
        split.addWidget(self.list)
        self.editor_area = QScrollArea()
        self.editor_area.setWidgetResizable(True)
        split.addWidget(self.editor_area)
        split.setSizes([220, 520])
        v.addWidget(split, 1)
        dock = QDockWidget("Pipeline", self)
        dock.setObjectName("pipeline")
        dock.setWidget(w)
        dock.setMinimumWidth(360)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.pipeline_dock = dock

    def _make_render_dock(self):
        s = self.render_settings
        w = QWidget()
        f = QFormLayout(w)
        self.r_bg = ColorButton(s.background)
        self.r_bg.changed.connect(lambda c: self._rset("background", c))
        f.addRow("Background", self.r_bg)
        presets = QHBoxLayout()
        for label, c in (("White", (1, 1, 1)), ("Black", (0, 0, 0)), ("Gray", (0.32, 0.34, 0.38))):
            b = QPushButton(label)
            b.clicked.connect(lambda _=False, c=c: (self.r_bg.set_color(c), self._rset("background", c)))
            presets.addWidget(b)
        f.addRow("", presets)
        self.r_radius = QDoubleSpinBox(minimum=0.05, maximum=10, singleStep=0.05, value=s.radius_scale)
        self.r_radius.valueChanged.connect(lambda x: self._rset("radius_scale", x))
        f.addRow("Radius scale", self.r_radius)
        self.r_res = QSpinBox(minimum=4, maximum=64, value=s.sphere_resolution)
        self.r_res.valueChanged.connect(lambda x: self._rset("sphere_resolution", x))
        f.addRow("Sphere quality", self.r_res)
        self.r_shading = QComboBox()
        self.r_shading.addItems(["smooth", "flat", "pbr"])
        self.r_shading.currentTextChanged.connect(lambda x: self._rset("shading", x))
        f.addRow("Shading", self.r_shading)
        self.r_mode = QComboBox()
        self.r_mode.addItems(["auto", "mesh", "imposter"])
        self.r_mode.setToolTip("mesh: geometric spheres (best for small systems); imposter: GPU spheres, fast for millions; "
                               "auto: imposters above 30k particles")
        self.r_mode.currentTextChanged.connect(lambda x: self._rset("render_mode", x))
        f.addRow("Spheres", self.r_mode)
        self.r_spec = QDoubleSpinBox(minimum=0, maximum=1, singleStep=0.05, value=s.specular)
        self.r_spec.valueChanged.connect(lambda x: self._rset("specular", x))
        f.addRow("Shininess", self.r_spec)
        self.r_roll = QDoubleSpinBox(minimum=-180, maximum=180, singleStep=90, value=s.view_roll, suffix=" deg")
        self.r_roll.setToolTip("Rotates the view presets in the image plane, e.g. 90 puts y horizontal in the top view")
        self.r_roll.valueChanged.connect(self._set_roll)
        f.addRow("View rotation", self.r_roll)
        for attr, label in (("show_cell", "Show cell"), ("show_bonds", "Show bonds"), ("show_legend", "Show legend"),
                            ("show_axes", "Axes in exported images"), ("ambient_occlusion", "Ambient occlusion"),
                            ("anti_aliasing", "Anti-aliasing"), ("parallel_projection", "Orthographic")):
            cb = QCheckBox()
            cb.setChecked(getattr(s, attr))
            cb.toggled.connect(lambda x, a=attr: self._rset(a, x))
            setattr(self, f"r_{attr}", cb)
            f.addRow(label, cb)
        self.r_cellw = QDoubleSpinBox(minimum=0.5, maximum=20, singleStep=0.5, value=s.cell_width)
        self.r_cellw.valueChanged.connect(lambda x: self._rset("cell_width", x))
        f.addRow("Cell line width", self.r_cellw)
        self.r_cellc = ColorButton(s.cell_color)
        self.r_cellc.changed.connect(lambda c: self._rset("cell_color", c))
        f.addRow("Cell color", self.r_cellc)
        self.r_font = QSpinBox(minimum=6, maximum=72, value=s.legend_font)
        self.r_font.valueChanged.connect(lambda x: self._rset("legend_font", x))
        f.addRow("Legend font", self.r_font)
        self.r_title = QLineEdit(s.title)
        self.r_title.setPlaceholderText("Optional title on the image")
        self.r_title.editingFinished.connect(lambda: self._rset("title", self.r_title.text()))
        f.addRow("Title", self.r_title)
        self.r_overlay = QLineEdit(s.overlay)
        self.r_overlay.setPlaceholderText("e.g. t = {Timestep*0.002:.1f} ps  N = {N}")
        self.r_overlay.setToolTip("Text updated every frame. {expression:format} uses numeric attributes (dots become "
                                  "underscores: CNA_counts_FCC), plus N, Frame and Timestep.")
        self.r_overlay.editingFinished.connect(lambda: self._rset("overlay", self.r_overlay.text()))
        f.addRow("Overlay", self.r_overlay)
        self.r_overlay_pos = QComboBox()
        self.r_overlay_pos.addItems(["lower_right", "lower_left", "upper_right", "upper_left", "upper_edge", "lower_edge"])
        self.r_overlay_pos.currentTextChanged.connect(lambda x: self._rset("overlay_position", x))
        f.addRow("Overlay position", self.r_overlay_pos)
        b1 = QPushButton("Render image…")
        b1.clicked.connect(self.export_image)
        b2 = QPushButton("Render animation…")
        b2.clicked.connect(self.export_animation)
        f.addRow(b1)
        f.addRow(b2)
        dock = QDockWidget("Render", self)
        dock.setObjectName("render")
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(w)
        dock.setWidget(scroll)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        self.render_dock = dock

    def _make_menus(self):
        mb = self.menuBar()
        fm = mb.addMenu("&File")
        self._act(fm, "Open file…", self.open_dialog, QKeySequence.Open)
        self._act(fm, "Open trajectory + topology…", self.open_with_topology)
        self.recent_menu = fm.addMenu("Open recent")
        self._refresh_recent()
        fm.addSeparator()
        self._act(fm, "Open session…", self.load_session, "Ctrl+Shift+O")
        self._act(fm, "Save session", self.save_session, QKeySequence.Save)
        self._act(fm, "Save session as…", lambda: self.save_session(True), QKeySequence.SaveAs)
        fm.addSeparator()
        self._act(fm, "Render image…", self.export_image, "Ctrl+E")
        self._act(fm, "Render animation…", self.export_animation, "Ctrl+Shift+E")
        self._act(fm, "Export data…", self.export_data)
        self._act(fm, "Batch render folder…", self.batch_render, "Ctrl+B")
        fm.addSeparator()
        self._act(fm, "Quit", self.close, QKeySequence.Quit)
        am = mb.addMenu("&Analysis")
        self._act(am, "Trajectory dynamics (MSD, alpha2, Fs, VACF)…", self.trajectory_msd, "Ctrl+M")
        vm = mb.addMenu("&View")
        for label, key, view in (("Isometric", "0", "iso"), ("Front (−Y)", "1", "-y"), ("Back (+Y)", "Ctrl+1", "+y"),
                                 ("Right (+X)", "3", "+x"), ("Left (−X)", "Ctrl+3", "-x"), ("Top (+Z)", "7", "+z"),
                                 ("Bottom (−Z)", "Ctrl+7", "-z")):
            self._act(vm, label, lambda v=view: self.set_view(v), key)
        vm.addSeparator()
        self._act(vm, "Zoom to fit", self.view.reset_camera, "F")
        vm.addSeparator()
        for d in (self.pipeline_dock, self.render_dock, self.inspector_dock):
            vm.addAction(d.toggleViewAction())
        hm = mb.addMenu("&Help")
        self._act(hm, "Scripting guide", self.show_help)
        self._act(hm, "About AtomView", lambda: QMessageBox.about(
            self, "AtomView", "<b>AtomView 0.1</b><br>Free, open-source visualization and analysis of atomistic "
            "simulations.<br>Built on VTK/PyVista, freud, ASE, NumPy and Qt."))

    def _act(self, menu, text, fn, key=None):
        a = QAction(text, self)
        if key:
            a.setShortcut(QKeySequence(key))
        a.triggered.connect(fn)
        menu.addAction(a)
        return a

    def _make_toolbar(self):
        tb = QToolBar("View")
        tb.setObjectName("viewbar")
        self.addToolBar(tb)
        tb.addAction("Open", self.open_dialog)
        tb.addSeparator()
        for label, view in (("Iso", "iso"), ("X", "+x"), ("Y", "-y"), ("Z", "+z")):
            tb.addAction(label, lambda v=view: self.set_view(v))
        tb.addAction("Fit", self.view.reset_camera)
        tb.addSeparator()
        tb.addAction("📷 Image", self.export_image)
        tb.addAction("🎞 Movie", self.export_animation)

    # ------------------------------------------------------------------ files
    def open_dialog(self):
        start = self.settings_store.value("last_dir", os.path.expanduser("~"))
        path, _ = QFileDialog.getOpenFileName(self, "Open particle file", start, DATA_FILTER)
        if path:
            self.open_path(path)

    def open_with_topology(self):
        start = self.settings_store.value("last_dir", os.path.expanduser("~"))
        path, _ = QFileDialog.getOpenFileName(self, "Trajectory", start, "Trajectories (*.xtc *.trr *.dcd *.nc *.tng *.trj);;All files (*)")
        if not path:
            return
        top, _ = QFileDialog.getOpenFileName(self, "Topology", os.path.dirname(path), "Topology (*.gro *.pdb *.psf *.tpr *.mol2);;All files (*)")
        if top:
            self.open_path(path, topology=top)

    def open_path(self, path, keep_modifiers=True, topology=None):
        try:
            QApplication.setOverrideCursor(Qt.WaitCursor)
            traj = load_trajectory(path, topology)
        except Exception as e:
            QMessageBox.critical(self, "Could not open file", f"{path}\n\n{e}")
            return
        finally:
            QApplication.restoreOverrideCursor()
        mods = self.pipeline.modifiers if (self.pipeline and keep_modifiers) else []
        self.pipeline = Pipeline(traj, mods)
        self.settings_store.setValue("last_dir", os.path.dirname(os.path.abspath(path)))
        self._remember(path)
        if not self.session_path:
            self._saved_state = self._state()
        self._update_title()
        self._first_view = True
        self._sync_timeline()
        self._rebuild_list()
        self.set_frame(0, force=True)

    def _refresh_recent(self):
        self.recent_menu.clear()
        for p in self.settings_store.value("recent", []) or []:
            self.recent_menu.addAction(p, lambda p=p: self.open_any(p))

    def _remember(self, path):
        path = os.path.abspath(path)
        recent = [p for p in (self.settings_store.value("recent", []) or []) if p != path]
        self.settings_store.setValue("recent", [path] + recent[:9])
        self._refresh_recent()

    def open_any(self, path):
        """Session files open as sessions, everything else as data (double-click, drop, recent, command line)."""
        if is_session(path):
            self.load_session(path)
        else:
            self.open_path(path)

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls():
            e.acceptProposedAction()

    def dropEvent(self, e):
        urls = e.mimeData().urls()
        if urls:
            self.open_any(urls[0].toLocalFile())

    def _session_extra(self):
        return {"render": self.render_settings.to_dict(), "camera": [list(c) for c in self.view.camera_position],
                "frame": self.frame_index, "fps": self.fps.value(), "selected": self.list.currentRow()}

    def _state(self):
        # Everything a session stores except the camera, to tell whether there are unsaved changes.
        import json
        if not self.pipeline:
            return None
        return json.dumps([self.pipeline.to_dict(), self.render_settings.to_dict(), self.frame_index, self.fps.value()],
                          sort_keys=True, default=str)

    def dirty(self):
        return self.pipeline is not None and self._state() != self._saved_state

    def _update_title(self):
        if not self.pipeline:
            return self.setWindowTitle("AtomView")
        name = os.path.basename(self.session_path) if self.session_path else "Untitled session"
        self.setWindowTitle(f"{name}{' *' if self.dirty() else ''} — {os.path.basename(self.pipeline.source.path)} — AtomView")
        self.setWindowModified(self.dirty())

    def save_session(self, as_new=False):
        if not self.pipeline:
            return False
        path = self.session_path
        if as_new or not path:
            start = os.path.join(os.path.dirname(os.path.abspath(self.pipeline.source.path)),
                                 os.path.splitext(os.path.basename(self.pipeline.source.path))[0] + ".atomview")
            path, _ = QFileDialog.getSaveFileName(self, "Save session", path or start, "AtomView session (*.atomview)")
            if not path:
                return False
            if not path.lower().endswith((".atomview", ".json")):
                path += ".atomview"
        try:
            self.pipeline.save(path, extra=self._session_extra())
        except OSError as e:
            QMessageBox.critical(self, "Could not save session", f"{path}\n\n{e}")
            return False
        self.session_path, self._saved_state = path, self._state()
        self._remember(path)
        self._update_title()
        self.statusBar().showMessage(f"Saved {path}", 5000)
        return True

    def load_session(self, path=None):
        if not path:
            start = self.settings_store.value("last_dir", os.path.expanduser("~"))
            path, _ = QFileDialog.getOpenFileName(self, "Open session", start, "AtomView session (*.atomview *.json)")
        if not path:
            return
        try:
            p = Pipeline.load(path)
        except FileNotFoundError as e:
            QMessageBox.information(self, "Data file moved", f"{e}\n\nPlease locate the data file.")
            data, _ = QFileDialog.getOpenFileName(self, "Locate the data file", os.path.dirname(path), DATA_FILTER)
            if not data:
                return
            try:
                p = Pipeline.load(path, source=load_trajectory(data))
            except Exception as e:
                return QMessageBox.critical(self, "Could not open session", str(e))
        except Exception as e:
            return QMessageBox.critical(self, "Could not open session", str(e))
        self.pipeline = p
        extra = getattr(p, "session_extra", {})
        self.render_settings = RenderSettings.from_dict(extra.get("render"))
        self.fps.setValue(int(extra.get("fps", self.fps.value())))
        self._sync_render_controls()
        self._sync_timeline()
        sel = int(extra.get("selected", -1))
        self._rebuild_list(select=sel if 0 <= sel <= len(p.modifiers) else None)
        self._first_view = "camera" not in extra
        self._session_camera = extra.get("camera")
        self.set_frame(int(extra.get("frame", 0)), force=True)
        self.session_path = path
        self.settings_store.setValue("last_dir", os.path.dirname(os.path.abspath(path)))
        self._remember(path)
        self._saved_state = self._state()
        self._update_title()

    def export_data(self):
        if not self.current:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export data", "output.xyz",
                                              "Extended XYZ (*.xyz);;LAMMPS dump (*.dump);;Other via ASE (*.*)")
        if not path:
            return
        frames = range(self.n_frames)
        if self.n_frames > 1 and QMessageBox.question(self, "Export", "Export all frames?\n(No = only the current frame)") \
                != QMessageBox.Yes:
            frames = [self.frame_index]
        for k, f in enumerate(frames):
            write_frame(self.pipeline.compute(f), path, append=k > 0)
        self.statusBar().showMessage(f"Wrote {len(frames)} frame(s) to {path}", 6000)

    # ------------------------------------------------------------------ pipeline UI
    def _rebuild_list(self, select=None):
        self.list.blockSignals(True)
        self.list.clear()
        if self.pipeline:
            src = QListWidgetItem(f"📂 {os.path.basename(self.pipeline.source.path)}")
            src.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled)
            self.list.addItem(src)
            for m in self.pipeline.modifiers:
                it = QListWidgetItem(self._item_text(m))
                it.setFlags(Qt.ItemIsSelectable | Qt.ItemIsEnabled | Qt.ItemIsUserCheckable)
                it.setCheckState(Qt.Checked if m.enabled else Qt.Unchecked)
                self.list.addItem(it)
        self.list.blockSignals(False)
        if select is not None:
            self.list.setCurrentRow(select)
        elif self.list.count():
            self.list.setCurrentRow(self.list.count() - 1)

    def _item_text(self, m: Modifier):
        return ("⚠ " if m.error else "") + m.label()

    def _refresh_items(self):
        if not self.pipeline:
            return
        self.list.blockSignals(True)
        for k, m in enumerate(self.pipeline.modifiers):
            it = self.list.item(k + 1)
            if it:
                it.setText(self._item_text(m))
                it.setForeground(Qt.red if m.error else self.list.palette().text().color())
        self.list.blockSignals(False)
        w = self.editor_area.widget()
        if isinstance(w, ModifierEditor):
            w.refresh_status()

    def _show_editor(self, row):
        if not self.pipeline or row < 0:
            self.editor_area.setWidget(QWidget())
            return
        if row == 0:
            self.editor_area.setWidget(self._source_panel())
            return
        m = self.pipeline.modifiers[row - 1]
        try:
            upstream = self.pipeline.compute(self.frame_index, stop=row - 1)
            names = list(scalar_columns(upstream))
        except Exception:
            names = []
        ed = ModifierEditor(m, names)
        ed.changed.connect(self._modifier_changed)
        self.editor_area.setWidget(ed)

    def _source_panel(self):
        w = QWidget()
        v = QVBoxLayout(w)
        src = self.pipeline.source
        f0 = src.get(self.frame_index)
        text = (f"<b>Data source</b><br>{src.path}<br><br>Format: {src.format}<br>Frames: {src.n_frames}<br>"
                f"{f0.summary()}<br>Periodic: {' '.join('yes' if p else 'no' for p in f0.pbc)}")
        if f0.type_names:
            text += "<br>Types: " + ", ".join(f"{k}={v}" for k, v in sorted(f0.type_names.items()))
        lab = QLabel(text)
        lab.setWordWrap(True)
        lab.setTextInteractionFlags(Qt.TextSelectableByMouse)
        v.addWidget(lab)
        b = QPushButton("Replace input file…")
        b.clicked.connect(self.open_dialog)
        v.addWidget(b)
        v.addStretch(1)
        return w

    def add_modifier(self, name):
        if not self.pipeline:
            QMessageBox.information(self, "AtomView", "Open a particle file first.")
            return
        m = REGISTRY[name]()
        row = self.list.currentRow()
        pos = len(self.pipeline.modifiers) if row < 0 else row  # insert after the selected item
        self.pipeline.modifiers.insert(pos, m)
        self.pipeline.invalidate()
        self._rebuild_list(select=pos + 1)
        self._schedule()

    def remove_modifier(self):
        row = self.list.currentRow()
        if self.pipeline and row > 0:
            del self.pipeline.modifiers[row - 1]
            self.pipeline.invalidate()
            self._rebuild_list(select=min(row, len(self.pipeline.modifiers)))
            self._schedule()

    def duplicate_modifier(self):
        row = self.list.currentRow()
        if self.pipeline and row > 0:
            m = Modifier.from_dict(self.pipeline.modifiers[row - 1].to_dict())
            self.pipeline.modifiers.insert(row, m)
            self.pipeline.invalidate()
            self._rebuild_list(select=row + 1)
            self._schedule()

    def move_modifier(self, d):
        row = self.list.currentRow()
        if not self.pipeline or row <= 0:
            return
        i, j = row - 1, row - 1 + d
        mods = self.pipeline.modifiers
        if 0 <= j < len(mods):
            mods[i], mods[j] = mods[j], mods[i]
            self.pipeline.invalidate()
            self._rebuild_list(select=j + 1)
            self._schedule()

    def _item_toggled(self, item):
        row = self.list.row(item)
        if row > 0:
            self.pipeline.modifiers[row - 1].enabled = item.checkState() == Qt.Checked
            self.pipeline.invalidate()
            self._schedule()

    def _modifier_changed(self):
        self.pipeline.invalidate()
        self._schedule()

    # ------------------------------------------------------------------ frames & compute
    @property
    def n_frames(self):
        return self.pipeline.n_frames if self.pipeline else 0

    def _sync_timeline(self):
        n = max(self.n_frames - 1, 0)
        for w in (self.slider, self.frame_spin):
            w.blockSignals(True)
            w.setRange(0, n)
            w.blockSignals(False)
        self.frame_label.setText(f"/ {n}")

    def set_frame(self, i, force=False):
        if not self.pipeline:
            return
        i = int(np.clip(i, 0, self.n_frames - 1))
        if i == self.frame_index and not force:
            return
        self.frame_index = i
        for w in (self.slider, self.frame_spin):
            w.blockSignals(True)
            w.setValue(i)
            w.blockSignals(False)
        self._schedule(immediate=True)

    def _schedule(self, immediate=False):
        if immediate:
            self._start_compute()
        else:
            self.recompute_timer.start()

    def _start_compute(self):
        if not self.pipeline:
            return
        if self._busy:
            self._pending = True
            return
        self._busy = True
        self._pending = False
        pipe, idx, rev = self.pipeline, self.frame_index, self.pipeline.revision
        pipe.cancel.clear()
        pipe.progress = lambda msg: self._signals.progress.emit(f"Frame {idx}: {msg}")
        self.busy_timer.start()

        def work():
            t0 = time.perf_counter()
            try:
                frame = pipe.compute(idx)
                self._signals.computed.emit(frame, idx, rev, time.perf_counter() - t0)
            except Cancelled:
                self._signals.failed.emit("Cancelled")
                self._signals.computed.emit(None, idx, -1, 0.0)
            except Exception as e:
                traceback.print_exc()
                self._signals.failed.emit(f"Error: {e}")
                self._signals.computed.emit(None, idx, rev, 0.0)

        threading.Thread(target=work, daemon=True).start()

    def cancel_compute(self):
        if self.pipeline:
            self.pipeline.cancel.set()
            self.busy_label.setText("Cancelling… waits for the running step to finish")
        if self.play_timer.isActive():
            self.toggle_play()

    def _on_computed(self, frame, idx, rev, seconds):
        self._busy = False
        self.busy_timer.stop()
        for w in (self.busy_label, self.busy_bar, self.busy_cancel):
            w.hide()
        if rev == -1:  # cancelled: do not recompute until something changes
            self._pending = False
            return
        if frame is not None and self.pipeline is not None:
            self.current = frame
            self._draw()
            self.inspector.set_frame(frame, idx)
            self._refresh_items()
            ts = frame.attributes.get("Timestep")
            msg = f"Frame {idx}" + (f" · timestep {ts}" if ts is not None else "") + f" · {frame.n} particles"
            msg += f" · {seconds * 1000:.0f} ms"
            errs = self.pipeline.errors()
            if errs:
                msg += f" · ⚠ {errs[0][0]}: {errs[0][1]}"
            self.statusBar().showMessage(msg)
        if self._pending or (self.pipeline and (rev != self.pipeline.revision or idx != self.frame_index)):
            self._start_compute()

    def _draw(self):
        cam = self.view.camera_position if not self._first_view else None
        build_scene(self.view, self.current, self.render_settings, interactive=True)
        if self._first_view:
            cam_saved = getattr(self, "_session_camera", None)
            if cam_saved:
                self.view.camera_position = [tuple(c) for c in cam_saved]
                self._session_camera = None
            else:
                default_camera(self.view, self.current, roll=self.render_settings.view_roll, rscale=self.render_settings.radius_scale)
            self._first_view = False
        else:
            self.view.camera_position = cam
        self.view.render()

    def toggle_play(self):
        if self.play_timer.isActive():
            self.play_timer.stop()
            self.btn_play.setIcon(self.style().standardIcon(QStyle.SP_MediaPlay))
        else:
            self.play_timer.start(int(1000 / self.fps.value()))
            self.btn_play.setIcon(self.style().standardIcon(QStyle.SP_MediaPause))

    def _advance(self):
        if self._busy:
            return
        self.play_timer.setInterval(int(1000 / self.fps.value()))
        nxt = self.frame_index + 1
        if nxt >= self.n_frames:
            if not self.loop.isChecked():
                self.toggle_play()
                return
            nxt = 0
        self.set_frame(nxt)

    # ------------------------------------------------------------------ view
    def set_view(self, view):
        if self.current is not None:
            default_camera(self.view, self.current, view, self.render_settings.view_roll, self.render_settings.radius_scale)
            self.view.render()

    def _set_roll(self, x):
        self.view.camera.Roll(x - self.render_settings.view_roll)
        self.render_settings.view_roll = x
        if self.current is not None: fit_camera(self.view, self.current, self.render_settings.radius_scale)
        self.view.render()

    def _rset(self, key, value):
        setattr(self.render_settings, key, value)
        if key in ("ambient_occlusion", "anti_aliasing"):
            apply_quality(self.view, self.render_settings)
        if self.current is not None:
            self._draw()

    def _sync_render_controls(self):
        s = self.render_settings
        self.r_bg.set_color(s.background)
        self.r_cellc.set_color(s.cell_color)
        for w, v in ((self.r_radius, s.radius_scale), (self.r_res, s.sphere_resolution), (self.r_spec, s.specular),
                     (self.r_cellw, s.cell_width), (self.r_font, s.legend_font), (self.r_roll, s.view_roll)):
            w.blockSignals(True)
            w.setValue(v)
            w.blockSignals(False)
        self.r_shading.setCurrentText(s.shading)
        self.r_mode.setCurrentText(s.render_mode)
        self.r_title.setText(s.title)
        self.r_overlay.setText(s.overlay)
        self.r_overlay_pos.setCurrentText(s.overlay_position)
        for attr in ("show_cell", "show_bonds", "show_legend", "show_axes", "ambient_occlusion", "anti_aliasing",
                     "parallel_projection"):
            cb = getattr(self, f"r_{attr}")
            cb.blockSignals(True)
            cb.setChecked(getattr(s, attr))
            cb.blockSignals(False)
        apply_quality(self.view, s)

    def _on_pick(self, point):
        if self.current is None or self.current.n == 0 or point is None:
            return
        from ..styles import particle_radii
        d = np.linalg.norm(self.current.positions - np.asarray(point), axis=1)
        k = int(np.argmin(d))
        if d[k] > 2.5 * particle_radii(self.current)[k] * self.render_settings.radius_scale:
            return
        cols = scalar_columns(self.current)
        parts = []
        for name in ("id", "type", "x", "y", "z"):
            if name in cols:
                parts.append(f"{name}={cols[name][k]:.4g}" if name in "xyz" else f"{name}={cols[name][k]}")
        for name, v in cols.items():
            if name not in ("id", "type", "x", "y", "z") and len(parts) < 12:
                val = v[k]
                parts.append(f"{name}={val:.4g}" if isinstance(val, (float, np.floating)) else f"{name}={val}")
        tname = self.current.type_names.get(int(cols.get("type", [0] * (k + 1))[k]), "")
        self.statusBar().showMessage(f"Particle {k}{' (' + tname + ')' if tname else ''}: " + ", ".join(parts))
        self.inspector.show_particle(k)

    # ------------------------------------------------------------------ export
    def _size_dialog(self, title, animation=False):
        dlg = QDialog(self)
        dlg.setWindowTitle(title)
        f = QFormLayout(dlg)
        w = QSpinBox(minimum=64, maximum=16000, value=int(self.settings_store.value("exp_w", 2400)))
        h = QSpinBox(minimum=64, maximum=16000, value=int(self.settings_store.value("exp_h", 1800)))
        f.addRow("Width (px)", w)
        f.addRow("Height (px)", h)
        match = QPushButton("Match viewport aspect")
        vp = self.view.interactor.size()
        match.clicked.connect(lambda: h.setValue(int(w.value() * vp.height() / max(vp.width(), 1))))
        f.addRow("", match)
        extras = {}
        if animation:
            extras["start"] = QSpinBox(minimum=0, maximum=max(self.n_frames - 1, 0), value=0)
            extras["end"] = QSpinBox(minimum=0, maximum=max(self.n_frames - 1, 0), value=max(self.n_frames - 1, 0))
            extras["every"] = QSpinBox(minimum=1, maximum=10000, value=1)
            extras["fps"] = QSpinBox(minimum=1, maximum=120, value=self.fps.value())
            for k, label in (("start", "First frame"), ("end", "Last frame"), ("every", "Every n-th"), ("fps", "fps")):
                f.addRow(label, extras[k])
        else:
            extras["transparent"] = QCheckBox("Transparent background (PNG)")
            f.addRow("", extras["transparent"])
            note = QLabel("Tip: save as .svg/.pdf/.eps for vector output (best without ambient occlusion).")
            note.setWordWrap(True)
            note.setStyleSheet("color: gray;")
            f.addRow(note)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        f.addRow(bb)
        if dlg.exec() != QDialog.Accepted:
            return None
        self.settings_store.setValue("exp_w", w.value())
        self.settings_store.setValue("exp_h", h.value())
        return (w.value(), h.value()), extras

    def export_image(self):
        if self.current is None:
            return
        res = self._size_dialog("Render image")
        if not res:
            return
        size, extras = res
        path, _ = QFileDialog.getSaveFileName(self, "Save image", "figure.png",
                                              "PNG (*.png);;JPEG (*.jpg);;TIFF (*.tif);;SVG (*.svg);;PDF (*.pdf);;EPS (*.eps)")
        if not path:
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            render_image(self.current, path, size=size, settings=self.render_settings,
                         camera=self.view.camera_position, transparent=extras["transparent"].isChecked())
            self.statusBar().showMessage(f"Saved {path}", 8000)
        except Exception as e:
            QMessageBox.critical(self, "Render failed", str(e))
        finally:
            QApplication.restoreOverrideCursor()

    def export_animation(self):
        if self.current is None:
            return
        res = self._size_dialog("Render animation", animation=True)
        if not res:
            return
        size, ex = res
        path, _ = QFileDialog.getSaveFileName(self, "Save animation", "movie.mp4", "MP4 (*.mp4);;GIF (*.gif);;AVI (*.avi)")
        if not path:
            return
        frames = list(range(ex["start"].value(), ex["end"].value() + 1, ex["every"].value()))
        prog = QProgressDialog("Rendering animation…", "Cancel", 0, len(frames), self)
        prog.setWindowModality(Qt.WindowModal)
        prog.show()

        def progress(i, n):
            prog.setValue(i)
            QApplication.processEvents()
            return not prog.wasCanceled()

        try:
            render_animation(self.pipeline, path, frames=frames, size=size, fps=ex["fps"].value(),
                             settings=self.render_settings, camera=self.view.camera_position, progress=progress)
            self.statusBar().showMessage(f"Saved {path}", 8000)
        except Exception as e:
            QMessageBox.critical(self, "Render failed", str(e))
        finally:
            prog.close()

    def trajectory_msd(self):
        from ..trajectory import dynamics
        if not self.pipeline or self.n_frames < 2:
            QMessageBox.information(self, "MSD", "Open a trajectory with at least two frames first.")
            return
        dlg = QDialog(self)
        dlg.setWindowTitle("Trajectory dynamics")
        f = QFormLayout(dlg)
        n = self.n_frames
        start, end = QSpinBox(minimum=0, maximum=n - 2, value=0), QSpinBox(minimum=1, maximum=n - 1, value=n - 1)
        every = QSpinBox(minimum=1, maximum=max(n - 1, 1), value=1)
        dt = QDoubleSpinBox(minimum=0, maximum=1e12, decimals=6, value=0)
        dt.setSpecialValueText("from Timestep")
        sel = QCheckBox("Only particles selected in the first frame")
        kval = QDoubleSpinBox(minimum=0, maximum=1e6, decimals=4, value=0)
        kval.setSpecialValueText("off")
        kval.setToolTip("Wave number for Fs(k,t), typically the first peak of S(k) (Analysis: Structure factor)")
        for label, w in (("First frame", start), ("Last frame", end), ("Every n-th", every), ("Time per frame", dt),
                         ("k for Fs(k,t)", kval), ("", sel)):
            f.addRow(label, w)
        note = QLabel("Uses the pipeline output of every frame. Periodic coordinates are unwrapped between frames, "
                      "which requires particles to move less than half the cell per frame.")
        note.setWordWrap(True)
        note.setStyleSheet("color: gray;")
        f.addRow(note)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        f.addRow(bb)
        if dlg.exec() != QDialog.Accepted:
            return
        frames = range(start.value(), end.value() + 1, every.value())
        prog = QProgressDialog("Computing dynamics…", "Cancel", 0, len(frames), self)
        prog.setWindowModality(Qt.WindowModal)
        prog.setMinimumDuration(300)

        def progress(k, n):
            prog.setValue(k)
            prog.setLabelText(f"Frame {k + 1} of {n}")
            QApplication.processEvents()
            return not prog.wasCanceled()

        try:
            res = dynamics(self.pipeline, frames, dt.value() or None, kval.value() or None, sel.isChecked(), progress=progress)
        except Cancelled:
            self.statusBar().showMessage("Cancelled", 5000)
            return
        except Exception as e:
            QMessageBox.critical(self, "Dynamics failed", str(e))
            return
        finally:
            prog.close()
        for name, tab in res.items():
            self.inspector.add_table(name, tab)
        self.inspector.table_choice.setCurrentText("MSD(t)")
        r = res["MSD(t)"]
        self.statusBar().showMessage(f"MSD over {r['frames']} frames, {r['particles']} particles: "
                                     f"D = {r['D']:.6g} (fit over the second half, {r['dim']}D)"
                                     + (f"; {r['dropped']} particles not present in every frame were left out" if r['dropped'] else ""), 15000)

    def batch_render(self):
        from .. import batch
        dlg = QDialog(self)
        dlg.setWindowTitle("Batch render folder")
        f = QFormLayout(dlg)
        info = QLabel("Goes through the folder and all sub-folders and renders images of the snapshots in each one "
                      "with the current pipeline and render settings.")
        info.setWordWrap(True)
        f.addRow(info)

        def folder_row(placeholder):
            line = QLineEdit()
            line.setPlaceholderText(placeholder)
            btn = QPushButton("…")
            btn.setFixedWidth(28)
            btn.clicked.connect(lambda: line.setText(QFileDialog.getExistingDirectory(self, "Choose folder",
                                                                                     line.text() or os.path.expanduser("~")) or line.text()))
            box = QWidget()
            h = QHBoxLayout(box)
            h.setContentsMargins(0, 0, 0, 0)
            h.addWidget(line, 1)
            h.addWidget(btn)
            return box, line

        root_box, root = folder_row("Folder to scan")
        root.setText(self.settings_store.value("batch_root", ""))
        f.addRow("Folder", root_box)
        out_box, out = folder_row("(empty = save each image next to its file)")
        out.setText(self.settings_store.value("batch_out", ""))
        f.addRow("Output folder", out_box)
        patterns = QLineEdit(self.settings_store.value("batch_patterns", "*.xyz, *.mat"))
        f.addRow("File patterns", patterns)
        pick = QComboBox()
        pick.addItems(["last by file name", "newest by modification time"])
        f.addRow("Which file", pick)
        files = QComboBox()
        files.addItems(["last matching file", "all matching files, in order"])
        f.addRow("Snapshots from", files)
        frames = QLineEdit(self.settings_store.value("batch_frames", "last"))
        frames.setToolTip("last, first, all, every5, ::5 (every 5th), -10: (last 10), 0,50,100")
        f.addRow("Snapshots", frames)
        use_pipe = QCheckBox("Apply the current pipeline's modifiers")
        use_pipe.setChecked(True)
        f.addRow("", use_pipe)
        recursive = QCheckBox("Include sub-folders")
        recursive.setChecked(True)
        f.addRow("", recursive)
        w = QSpinBox(minimum=64, maximum=16000, value=int(self.settings_store.value("exp_w", 2400)))
        h = QSpinBox(minimum=64, maximum=16000, value=int(self.settings_store.value("exp_h", 1800)))
        f.addRow("Width (px)", w)
        f.addRow("Height (px)", h)
        view = QComboBox()
        view.addItems(["iso", "+x", "-x", "+y", "-y", "+z", "-z"])
        f.addRow("View", view)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        f.addRow(bb)
        if dlg.exec() != QDialog.Accepted or not root.text().strip():
            return
        for k, v in (("batch_root", root.text()), ("batch_out", out.text()), ("batch_patterns", patterns.text()),
                     ("batch_frames", frames.text())):
            self.settings_store.setValue(k, v)
        pats = [p.strip() for p in patterns.text().split(",") if p.strip()]
        mods = [m.to_dict() for m in self.pipeline.modifiers] if (self.pipeline and use_pipe.isChecked()) else []
        n = sum(1 for d in batch.folders(root.text(), recursive.isChecked()) if batch.latest_file(d, pats))
        one = batch.single(frames.text()) and not files.currentIndex()
        prog = QProgressDialog("Rendering…", "Cancel", 0, max(n, 1) if one else 0, self)
        prog.setWindowModality(Qt.WindowModal)
        prog.setMinimumWidth(520)
        prog.show()
        lines = []

        def progress(msg):
            lines.append(msg)
            prog.setLabelText(msg[-120:])
            if one: prog.setValue(min(len(lines), n))
            QApplication.processEvents()

        results = batch.run(root.text(), patterns=pats, out_dir=out.text().strip() or None, modifiers=mods,
                            settings=self.render_settings, pick="mtime" if pick.currentIndex() else "name",
                            size=(w.value(), h.value()), view=view.currentText(), recursive=recursive.isChecked(),
                            progress=progress, cancel=prog.wasCanceled, frames=frames.text().strip() or "last",
                            files="all" if files.currentIndex() else "last")
        prog.close()
        failed = [r for r in results if r["status"].startswith("error")]
        msg = f"{len(results) - len(failed)} image(s) written, {len(failed)} failed.\n\n"
        msg += "\n".join(f"{os.path.basename(r['file'])}: {r['status']}" for r in failed[:15])
        msg += f"\n\nSummary: {os.path.join(out.text().strip() or root.text(), 'batch_summary.csv')}"
        QMessageBox.information(self, "Batch render", msg)

    def show_help(self):
        from ..modifiers.scripting import MATLAB_TEMPLATE, PYTHON_TEMPLATE
        text = ("<h3>Scripting in AtomView</h3>"
                "<p><b>Python script modifier</b> — define <code>modify(frame, data)</code>. It runs on every "
                "frame; change <code>data</code> in place or return a new frame.</p><pre>" + PYTHON_TEMPLATE +
                "</pre><p><b>MATLAB / Octave script modifier</b> — runs a .m script with MATLAB if installed, "
                "otherwise the free GNU Octave (<code>brew install octave</code>).</p><pre>" + MATLAB_TEMPLATE +
                "</pre><p><b>Expressions</b> (selection, compute property, inspector filter) accept "
                "<code>type == 2 &amp;&amp; z &gt; 10</code>, <code>abs(c_pe) &gt; 3</code>, "
                "<code>sqrt(vx^2+vy^2)</code>.</p><p><b>Batch mode</b>: <code>atomview render file.dump -s "
                "session.atomview -o fig.png</code>, or <code>import atomview</code> in your own Python.</p>")
        box = QMessageBox(self)
        box.setWindowTitle("Scripting guide")
        box.setTextFormat(Qt.RichText)
        box.setText(text)
        box.exec()

    def closeEvent(self, e):
        if self.dirty():
            b = QMessageBox.question(self, "Unsaved changes", "Save the session before closing?",
                                     QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel, QMessageBox.Save)
            if b == QMessageBox.Cancel or (b == QMessageBox.Save and not self.save_session()):
                return e.ignore()
        self.play_timer.stop()
        self.settings_store.setValue("window_state", self.saveState())
        self.view.close()
        super().closeEvent(e)


class App(QApplication):
    """Receives the files macOS hands over on double-click or drop onto the Dock icon."""
    def __init__(self, argv):
        super().__init__(argv)
        self.window, self.queued = None, []

    def event(self, e):
        from PySide6.QtCore import QEvent
        if e.type() == QEvent.FileOpen:
            f = e.file()
            self.window.open_any(f) if self.window else self.queued.append(f)
            return True
        return super().event(e)


def icon_path():
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "resources", "icon.png")


def run(path=None):
    import sys
    from PySide6.QtGui import QIcon
    app = QApplication.instance() or App(sys.argv)
    app.setApplicationName("AtomView")
    if os.path.exists(icon_path()):
        app.setWindowIcon(QIcon(icon_path()))
    from ..plugins import ERRORS, load_plugins
    load_plugins()
    win = MainWindow(path)
    win.show()
    if isinstance(app, App):
        app.window = win
        for f in app.queued:
            QTimer.singleShot(0, lambda f=f: win.open_any(f))
    if ERRORS:
        QMessageBox.warning(win, "Plugins", "Some plugins failed to load:\n\n" + "\n\n".join(f"{n}\n{e}" for n, e in ERRORS)[:3000])
    return app.exec()
