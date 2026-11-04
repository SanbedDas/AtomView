"""Auto-generated parameter forms for modifiers."""
from __future__ import annotations

import os

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QCheckBox, QColorDialog, QComboBox, QCompleter, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit, QPushButton,
                               QSpinBox, QVBoxLayout, QWidget)

from ..pipeline import Modifier
from .highlighter import CodeEditor


def _color_style(c):
    r, g, b = [int(255 * v) for v in c[:3]]
    return f"background-color: rgb({r},{g},{b}); border: 1px solid #888; border-radius: 3px;"


class ColorButton(QPushButton):
    changed = Signal(tuple)

    def __init__(self, color):
        super().__init__()
        self.setFixedSize(48, 22)
        self.set_color(color)
        self.clicked.connect(self._pick)

    def set_color(self, c):
        self.color = tuple(float(v) for v in c[:3])
        self.setStyleSheet(_color_style(self.color))

    def _pick(self):
        q = QColorDialog.getColor(QColor.fromRgbF(*self.color), self, "Choose color")
        if q.isValid():
            self.set_color((q.redF(), q.greenF(), q.blueF()))
            self.changed.emit(self.color)


class ModifierEditor(QWidget):
    """Edits one modifier's values. Emits `changed` after any edit is committed."""

    changed = Signal()

    def __init__(self, modifier: Modifier, property_names=(), parent=None):
        super().__init__(parent)
        self.m = modifier
        outer = QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        title = QLabel(f"<b>{modifier.title}</b>")
        outer.addWidget(title)
        if modifier.description:
            d = QLabel(modifier.description)
            d.setWordWrap(True)
            d.setStyleSheet("color: gray;")
            outer.addWidget(d)
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        outer.addLayout(form)
        self.code_editor = None
        self.log_view = None
        for p in modifier.params:
            w = self._make(p, property_names)
            if w is None:
                continue
            if p.kind in ("code", "matlab"):
                outer.addWidget(QLabel(p.label + ":"))
                outer.addWidget(w, 1)
            else:
                form.addRow(p.label, w)
            if p.help:
                w.setToolTip(p.help)
        if self.code_editor is not None:
            row = QHBoxLayout()
            apply = QPushButton("Apply script  (⌘↵)")
            apply.clicked.connect(self._commit_code)
            row.addWidget(apply)
            if any(p.kind == "matlab" for p in modifier.params):
                load = QPushButton("Load .m…")
                load.clicked.connect(lambda: self._load_into_editor("MATLAB (*.m)"))
                row.addWidget(load)
            else:
                load = QPushButton("Load .py…")
                load.clicked.connect(lambda: self._load_into_editor("Python (*.py)"))
                row.addWidget(load)
            save = QPushButton("Save as…")
            save.clicked.connect(self._save_code)
            row.addWidget(save)
            outer.addLayout(row)
            self.log_view = QPlainTextEdit()
            self.log_view.setReadOnly(True)
            self.log_view.setMaximumHeight(110)
            self.log_view.setPlaceholderText("Script output and errors appear here")
            outer.addWidget(self.log_view)
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.status.setTextInteractionFlags(Qt.TextSelectableByMouse)
        outer.addWidget(self.status)
        if self.code_editor is None:
            outer.addStretch(1)
        self.refresh_status()

    # ---- widgets ------------------------------------------------------------
    def _set(self, name, value):
        if self.m.values.get(name) != value:
            self.m.values[name] = value
            self.changed.emit()

    def _make(self, p, property_names):
        v = self.m.values[p.name]
        if p.kind == "float":
            w = QDoubleSpinBox()
            w.setDecimals(4)
            w.setRange(p.min if p.min is not None else -1e9, p.max if p.max is not None else 1e9)
            w.setSingleStep(p.step or 0.1)
            w.setValue(float(v))
            w.setKeyboardTracking(False)
            w.valueChanged.connect(lambda x, n=p.name: self._set(n, float(x)))
            if p.name in ("start", "end"):
                self._range_widgets = getattr(self, "_range_widgets", {})
                self._range_widgets[p.name] = w
            return w
        if p.kind == "int":
            w = QSpinBox()
            w.setRange(int(p.min if p.min is not None else -10**9), int(p.max if p.max is not None else 10**9))
            w.setSingleStep(int(p.step or 1))
            w.setValue(int(v))
            w.setKeyboardTracking(False)
            w.valueChanged.connect(lambda x, n=p.name: self._set(n, int(x)))
            return w
        if p.kind == "bool":
            w = QCheckBox()
            w.setChecked(bool(v))
            w.toggled.connect(lambda x, n=p.name: self._set(n, bool(x)))
            return w
        if p.kind == "choice":
            w = QComboBox()
            w.addItems([str(c) for c in p.choices])
            w.setCurrentText(str(v))
            w.currentTextChanged.connect(lambda x, n=p.name: self._set(n, x))
            return w
        if p.kind == "property":
            w = QComboBox()
            w.setEditable(True)
            w.addItems(list(property_names))
            w.setCurrentText(str(v))
            w.currentTextChanged.connect(lambda x, n=p.name: self._set(n, x) if x else None)
            return w
        if p.kind in ("str", "expr"):
            w = QLineEdit(str(v))
            if p.kind == "expr":
                w.setPlaceholderText("e.g. type == 2 && z > 10")
                comp = QCompleter(list(property_names) + ["N", "Frame", "CellSize_X", "CellSize_Y", "CellSize_Z",
                                                          "abs", "sqrt", "min", "max"])
                comp.setCaseSensitivity(Qt.CaseInsensitive)
                w.setCompleter(comp)
            w.editingFinished.connect(lambda n=p.name, w=w: self._set(n, w.text()))
            return w
        if p.kind == "color":
            w = ColorButton(v)
            w.changed.connect(lambda c, n=p.name: self._set(n, c))
            return w
        if p.kind == "vector":
            box = QWidget()
            h = QHBoxLayout(box)
            h.setContentsMargins(0, 0, 0, 0)
            spins = []
            for i in range(3):
                s = QDoubleSpinBox()
                s.setRange(-1e6, 1e6)
                s.setDecimals(3)
                s.setValue(float(v[i]))
                s.setKeyboardTracking(False)
                spins.append(s)
                h.addWidget(s)
            for s in spins:
                s.valueChanged.connect(lambda _x, n=p.name: self._set(n, tuple(sp.value() for sp in spins)))
            return box
        if p.kind == "file":
            box = QWidget()
            h = QHBoxLayout(box)
            h.setContentsMargins(0, 0, 0, 0)
            line = QLineEdit(str(v))
            line.setPlaceholderText("(none)")
            btn = QPushButton("…")
            btn.setFixedWidth(28)
            clear = QPushButton("×")
            clear.setFixedWidth(28)
            h.addWidget(line, 1)
            h.addWidget(btn)
            h.addWidget(clear)
            filt = "MATLAB (*.m)" if any(q.kind == "matlab" for q in self.m.params) else "Python (*.py)"

            def browse(n=p.name):
                path, _ = QFileDialog.getOpenFileName(self, "Choose script", os.path.expanduser("~"), filt)
                if path:
                    line.setText(path)
                    self._set(n, path)

            btn.clicked.connect(browse)
            clear.clicked.connect(lambda n=p.name: (line.setText(""), self._set(n, "")))
            line.editingFinished.connect(lambda n=p.name: self._set(n, line.text()))
            return box
        if p.kind in ("code", "matlab"):
            w = CodeEditor("matlab" if p.kind == "matlab" else "python")
            w.setPlainText(str(v))
            w.setMinimumHeight(220)
            w.apply_requested.connect(self._commit_code)
            self.code_editor = (p.name, w)
            return w
        return None

    def _commit_code(self):
        name, w = self.code_editor
        self.m.values[name] = w.toPlainText()
        if hasattr(self.m, "_compiled_for"):
            self.m._compiled_for = None
        self.changed.emit()

    def _load_into_editor(self, filt):
        path, _ = QFileDialog.getOpenFileName(self, "Load script", os.path.expanduser("~"), filt)
        if path:
            with open(path) as fh:
                self.code_editor[1].setPlainText(fh.read())
            self._commit_code()

    def _save_code(self):
        ext = ".m" if self.code_editor[1].highlighter and any(q.kind == "matlab" for q in self.m.params) else ".py"
        path, _ = QFileDialog.getSaveFileName(self, "Save script", os.path.expanduser(f"~/script{ext}"))
        if path:
            with open(path, "w") as fh:
                fh.write(self.code_editor[1].toPlainText())

    def refresh_status(self):
        m = self.m
        if m.error:
            self.status.setText(f"<span style='color:#d33'>⚠ {m.error}</span>")
        else:
            self.status.setText(f"<span style='color:gray'>{m.status}</span>" if m.status else "")
        if self.log_view is not None:
            self.log_view.setPlainText(m.log or "")
        for name, w in getattr(self, "_range_widgets", {}).items():
            if m.values.get("auto_range"):
                w.blockSignals(True)
                w.setValue(float(m.values[name]))
                w.blockSignals(False)
