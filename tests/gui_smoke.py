"""GUI smoke test: opens the window, adds every modifier, takes screenshots, quits.

    python tests/gui_smoke.py OUTDIR
"""
import os
import sys
import time

from PySide6.QtWidgets import QApplication

from atomview.gui.main_window import MainWindow
from atomview.pipeline import REGISTRY

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = sys.argv[1] if len(sys.argv) > 1 else "."
os.makedirs(out, exist_ok=True)
app = QApplication(sys.argv)
win = MainWindow()
win.show()


def settle(timeout=30):
    t0 = time.time()
    while time.time() - t0 < timeout:
        app.processEvents()
        time.sleep(0.02)
        if not win._busy and not win.recompute_timer.isActive():
            for _ in range(10):
                app.processEvents()
                time.sleep(0.02)
            if not win._busy:
                return
    raise TimeoutError("pipeline did not settle")


win.open_path(os.path.join(ROOT, "examples", "data", "cu_void.dump"))
settle()
problems = []
skip = {"Function", "MatlabScript", "Hexatic", "Vectors", "OrientationalOrder", "BondOrderCorrelation"}  # need 2D data or an orientation; covered in test_2d.py and test_softmatter.py
for name in REGISTRY:
    if name in skip:
        continue
    win.add_modifier(name)
    settle(60)
    m = win.pipeline.modifiers[win.list.currentRow() - 1]
    if m.error:
        problems.append((name, m.error))
    win.remove_modifier()
    settle()
print("modifier problems:", problems)

for name in ("CommonNeighborAnalysis", "Slice", "Coordination"):
    win.add_modifier(name)
    settle()
win.pipeline.modifiers[1].values.update(normal=(0, 1, 0), width=8.0)
win.pipeline.invalidate()
win._rebuild_list(select=2)
win._schedule(immediate=True)
settle()
win.set_frame(7)
settle()
from atomview.trajectory import msd
win.inspector.add_table("MSD(t)", msd(win.pipeline, range(0, 20)))
print("plots:", [win.inspector.table_choice.itemText(i) for i in range(win.inspector.table_choice.count())])
win.inspector.setCurrentIndex(2)
app.processEvents()
win.view.screenshot(os.path.join(out, "gui_view.png"))
win.grab().save(os.path.join(out, "gui_window.png"))
import time as _t
from atomview.modifiers import PythonScript
win.pipeline.modifiers.append(PythonScript(code="import time\ndef modify(frame, data):\n    time.sleep(5)\n"))
win.pipeline.invalidate(); win._rebuild_list(); win._schedule(immediate=True)
t0 = _t.time()
while _t.time() - t0 < 1.0:
    app.processEvents(); _t.sleep(0.02)
shown = win.busy_cancel.isVisible() and "Python script" in win.busy_label.text()
win.busy_cancel.click()
for _ in range(20): app.processEvents(); _t.sleep(0.02)
label = win.busy_label.text()
t0 = _t.time()
while win._busy and _t.time() - t0 < 10: app.processEvents(); _t.sleep(0.05)
print("cancel label:", label, "| busy after the script returned:", win._busy, "| cached:", (7, len(win.pipeline.modifiers)) in win.pipeline._cache)
print("busy bar shown:", shown, "| after cancel busy:", win._busy, "| status:", win.statusBar().currentMessage())
win.pipeline.modifiers.pop(); win.pipeline.invalidate(); win._rebuild_list(); win._schedule(immediate=True); settle()
print("frame", win.frame_index, "particles", win.current.n, "errors", win.pipeline.errors())
win.close()
