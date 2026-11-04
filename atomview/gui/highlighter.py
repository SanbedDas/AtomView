"""Code editor widget with Python / MATLAB syntax highlighting and line numbers."""
from __future__ import annotations

import re

from PySide6.QtCore import QRect, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase, QKeySequence, QPainter, QShortcut, QSyntaxHighlighter, QTextCharFormat
from PySide6.QtWidgets import QPlainTextEdit, QWidget

PY_KEYWORDS = ("and as assert break class continue def del elif else except False finally for from global if import "
               "in is lambda None nonlocal not or pass raise return True try while with yield").split()
M_KEYWORDS = ("break case catch classdef continue else elseif end for function global if otherwise parfor "
              "persistent return spmd switch try while true false").split()


def _fmt(color, bold=False, italic=False):
    f = QTextCharFormat()
    f.setForeground(QColor(color))
    if bold:
        f.setFontWeight(QFont.Bold)
    f.setFontItalic(italic)
    return f


class Highlighter(QSyntaxHighlighter):
    def __init__(self, doc, language="python", dark=False):
        super().__init__(doc)
        kw = PY_KEYWORDS if language == "python" else M_KEYWORDS
        comment = r"#[^\n]*" if language == "python" else r"%[^\n]*"
        pal = {"kw": "#c678dd" if dark else "#7a1fa2", "num": "#d19a66" if dark else "#a05a00",
               "str": "#98c379" if dark else "#2a7d2a", "com": "#7f848e" if dark else "#8a8a8a",
               "fn": "#61afef" if dark else "#1a5fb4"}
        self.rules = [
            (re.compile(r"\b(" + "|".join(kw) + r")\b"), _fmt(pal["kw"], bold=True)),
            (re.compile(r"\b[0-9]+(\.[0-9]*)?([eE][-+]?[0-9]+)?\b"), _fmt(pal["num"])),
            (re.compile(r"\b[A-Za-z_]\w*(?=\()"), _fmt(pal["fn"])),
            (re.compile(r"\"[^\"\n]*\"|'[^'\n]*'"), _fmt(pal["str"])),
            (re.compile(comment), _fmt(pal["com"], italic=True)),
        ]

    def highlightBlock(self, text):
        for rx, fmt in self.rules:
            for m in rx.finditer(text):
                self.setFormat(m.start(), m.end() - m.start(), fmt)


class _LineNumbers(QWidget):
    def __init__(self, editor):
        super().__init__(editor)
        self.editor = editor

    def sizeHint(self):
        return QSize(self.editor.gutter_width(), 0)

    def paintEvent(self, event):
        self.editor.paint_gutter(event)


class CodeEditor(QPlainTextEdit):
    """Plain-text code editor. Emits `apply_requested` on Ctrl/Cmd+Enter."""

    apply_requested = Signal()

    def __init__(self, language="python", parent=None):
        super().__init__(parent)
        font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
        font.setPointSize(12)
        self.setFont(font)
        self.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.setTabStopDistance(4 * self.fontMetrics().horizontalAdvance(" "))
        dark = self.palette().base().color().lightness() < 128
        self.highlighter = Highlighter(self.document(), language, dark)
        self.gutter = _LineNumbers(self)
        self.blockCountChanged.connect(lambda _: self.setViewportMargins(self.gutter_width(), 0, 0, 0))
        self.updateRequest.connect(self._update_gutter)
        self.setViewportMargins(self.gutter_width(), 0, 0, 0)
        for seq in ("Ctrl+Return", "Meta+Return"):
            QShortcut(QKeySequence(seq), self, activated=self.apply_requested.emit)

    def keyPressEvent(self, e):
        if e.key() == Qt.Key_Tab:
            self.insertPlainText("    ")
            return
        if e.key() in (Qt.Key_Return, Qt.Key_Enter) and not e.modifiers():
            line = self.textCursor().block().text()
            indent = re.match(r"\s*", line).group(0)
            if line.rstrip().endswith(":"):
                indent += "    "
            super().keyPressEvent(e)
            self.insertPlainText(indent)
            return
        super().keyPressEvent(e)

    def gutter_width(self):
        digits = len(str(max(1, self.blockCount())))
        return 10 + self.fontMetrics().horizontalAdvance("9") * digits

    def resizeEvent(self, e):
        super().resizeEvent(e)
        cr = self.contentsRect()
        self.gutter.setGeometry(QRect(cr.left(), cr.top(), self.gutter_width(), cr.height()))

    def _update_gutter(self, rect, dy):
        if dy:
            self.gutter.scroll(0, dy)
        else:
            self.gutter.update(0, rect.y(), self.gutter.width(), rect.height())

    def paint_gutter(self, event):
        painter = QPainter(self.gutter)
        base = self.palette().base().color()
        painter.fillRect(event.rect(), base.darker(106) if base.lightness() > 128 else base.lighter(130))
        block = self.firstVisibleBlock()
        top = int(self.blockBoundingGeometry(block).translated(self.contentOffset()).top())
        painter.setPen(QColor("#999999"))
        while block.isValid() and top <= event.rect().bottom():
            if block.isVisible():
                painter.drawText(0, top, self.gutter.width() - 4, self.fontMetrics().height(), Qt.AlignRight,
                                 str(block.blockNumber() + 1))
            top += int(self.blockBoundingRect(block).height())
            block = block.next()
