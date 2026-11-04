"""Data inspector: particle table, global attributes and plots of data tables."""
from __future__ import annotations

import numpy as np
from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PySide6.QtWidgets import (QComboBox, QFileDialog, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QPushButton,
                               QTableView, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from ..data import Frame, scalar_columns
from ..expressions import evaluate


class ParticleModel(QAbstractTableModel):
    def __init__(self):
        super().__init__()
        self.cols, self.names, self.rows = {}, [], np.zeros(0, dtype=int)

    def set_frame(self, frame: Frame | None, rows=None):
        self.beginResetModel()
        if frame is None:
            self.cols, self.names, self.rows = {}, [], np.zeros(0, dtype=int)
        else:
            self.cols = scalar_columns(frame)
            front = [k for k in ("id", "type", "x", "y", "z") if k in self.cols]
            self.names = front + [k for k in self.cols if k not in front]
            self.rows = np.arange(frame.n) if rows is None else rows
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return len(self.names)

    def data(self, index, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        v = self.cols[self.names[index.column()]][self.rows[index.row()]]
        if isinstance(v, (float, np.floating)):
            return f"{v:.6g}"
        return str(v)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        if orientation == Qt.Horizontal:
            return self.names[section]
        return str(self.rows[section])


class Inspector(QTabWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.frame = None
        self.extra_tables = {}  # trajectory-level results (e.g. MSD) shown alongside per-frame tables
        self.frame_index = 0

        # Particles
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(4, 4, 4, 4)
        row = QHBoxLayout()
        self.filter = QLineEdit()
        self.filter.setPlaceholderText("Filter rows, e.g. selection || c_pe > -3")
        self.filter.editingFinished.connect(self._apply_filter)
        self.count = QLabel()
        row.addWidget(self.filter, 1)
        row.addWidget(self.count)
        csv = QPushButton("CSV…")
        csv.clicked.connect(self._export_particles)
        row.addWidget(csv)
        v.addLayout(row)
        self.model = ParticleModel()
        self.table = QTableView()
        self.table.setModel(self.model)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.table.verticalHeader().setDefaultSectionSize(20)
        v.addWidget(self.table)
        self.addTab(page, "Particles")

        # Attributes
        self.attrs = QTableWidget(0, 2)
        self.attrs.setHorizontalHeaderLabels(["Attribute", "Value"])
        self.attrs.horizontalHeader().setStretchLastSection(True)
        self.attrs.verticalHeader().hide()
        self.addTab(self.attrs, "Attributes")

        # Plots
        page = QWidget()
        v = QVBoxLayout(page)
        v.setContentsMargins(4, 4, 4, 4)
        row = QHBoxLayout()
        self.table_choice = QComboBox()
        self.table_choice.currentTextChanged.connect(self._draw_plot)
        row.addWidget(self.table_choice, 1)
        for label, fn in (("Save plot…", self._save_plot), ("CSV…", self._export_table)):
            b = QPushButton(label)
            b.clicked.connect(fn)
            row.addWidget(b)
        v.addLayout(row)
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure
        self.figure = Figure(figsize=(4, 3), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.figure)
        v.addWidget(self.canvas, 1)
        self.plot_hint = QLabel("Add Coordination analysis (RDF), Histogram or Cluster analysis to get plots.")
        self.plot_hint.setStyleSheet("color: gray;")
        self.plot_hint.setWordWrap(True)
        v.addWidget(self.plot_hint)
        self.addTab(page, "Plots")

    def tables(self):
        return {**(self.frame.tables if self.frame is not None else {}), **self.extra_tables}

    def add_table(self, name, table):
        self.extra_tables[name] = table
        if self.frame is not None:
            self.set_frame(self.frame, self.frame_index)
        self.table_choice.setCurrentText(name)
        self.setCurrentIndex(2)

    def set_frame(self, frame: Frame, frame_index=0):
        self.frame, self.frame_index = frame, frame_index
        self._apply_filter()
        attrs = {k: v for k, v in frame.attributes.items() if not k.startswith("_")}
        attrs = {"Particles": frame.n, "Cell volume": round(frame.volume(), 6), **attrs}
        self.attrs.setRowCount(len(attrs))
        for r, (k, val) in enumerate(attrs.items()):
            self.attrs.setItem(r, 0, QTableWidgetItem(str(k)))
            self.attrs.setItem(r, 1, QTableWidgetItem(f"{val:.8g}" if isinstance(val, float) else str(val)))
        current = self.table_choice.currentText()
        self.table_choice.blockSignals(True)
        self.table_choice.clear()
        self.table_choice.addItems(list(self.tables()))
        if current in self.tables():
            self.table_choice.setCurrentText(current)
        self.table_choice.blockSignals(False)
        self.plot_hint.setVisible(not self.tables())
        self._draw_plot()

    def _apply_filter(self):
        if self.frame is None:
            return
        text = self.filter.text().strip()
        rows = None
        self.filter.setStyleSheet("")
        if text:
            try:
                rows = np.nonzero(evaluate(text, self.frame, self.frame_index).astype(bool))[0]
            except Exception as e:
                self.filter.setStyleSheet("border: 1px solid #d33;")
                self.filter.setToolTip(str(e))
        self.model.set_frame(self.frame, rows)
        self.count.setText(f"{self.model.rowCount()} rows")

    def show_particle(self, index):
        self.setCurrentIndex(0)
        if self.filter.text():
            self.filter.clear()
            self._apply_filter()
        self.table.selectRow(int(index))
        self.table.scrollTo(self.model.index(int(index), 0))

    def _draw_plot(self, *_):
        self.figure.clear()
        name = self.table_choice.currentText()
        if self.frame is None or name not in self.tables():
            self.canvas.draw_idle()
            return
        t = self.tables()[name]
        ax = self.figure.add_subplot(111)
        x, y = np.asarray(t["x"]), np.asarray(t["y"])
        if t.get("kind") == "bar":
            width = (x[1] - x[0]) if len(x) > 1 else 1
            ax.bar(x, y, width=width * 0.9, color="#4c78a8")
        else:
            ax.plot(x, y, color="#4c78a8", lw=1.6)
        ax.set_xlabel(t.get("xlabel", ""))
        ax.set_ylabel(t.get("ylabel", ""))
        ax.set_title(name, fontsize=10)
        ax.grid(alpha=0.25)
        self.canvas.draw_idle()

    def _save_plot(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save plot", "plot.png", "Images (*.png *.pdf *.svg)")
        if path:
            self.figure.savefig(path, dpi=300)

    def _export_table(self):
        name = self.table_choice.currentText()
        if not name or self.frame is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export table", f"{name.split(':')[0]}.csv", "CSV (*.csv)")
        if path:
            t = self.tables()[name]
            np.savetxt(path, np.column_stack([t["x"], t["y"]]), delimiter=",",
                       header=f"{t.get('xlabel', 'x')},{t.get('ylabel', 'y')}", comments="")

    def _export_particles(self):
        if self.frame is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export particles", "particles.csv", "CSV (*.csv)")
        if path:
            names = [n for n in self.model.names if np.asarray(self.model.cols[n]).dtype.kind in "biuf"]
            data = np.column_stack([self.model.cols[n][self.model.rows] for n in names])
            np.savetxt(path, data, delimiter=",", header=",".join(names), comments="", fmt="%.8g")
