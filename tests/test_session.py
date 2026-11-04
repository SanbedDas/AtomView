import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

import atomview as av
from atomview import modifiers as m
from atomview.pipeline import Pipeline, is_session

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DUMP = os.path.join(ROOT, "examples", "data", "cu_void.dump")


def session(tmp_path):
    proj = tmp_path / "proj"
    (proj / "data").mkdir(parents=True)
    shutil.copy(DUMP, proj / "data" / "run.dump")
    p = Pipeline(av.load(str(proj / "data" / "run.dump")), [m.ColorCoding(property="z"), m.Slice(width=4.0)])
    p.modifiers[1].enabled = False
    p.save(str(proj / "s.atomview"), extra={"frame": 7, "fps": 25})
    return proj, p


def test_session_round_trip(tmp_path):
    proj, p = session(tmp_path)
    d = json.load(open(proj / "s.atomview"))
    assert d["atomview_session"] == 1 and d["source_relative"] == os.path.join("data", "run.dump") and d["fps"] == 25
    q = Pipeline.load(str(proj / "s.atomview"))
    js = lambda mods: json.dumps([x.to_dict() for x in mods])
    assert js(q.modifiers) == js(p.modifiers) and q.session_extra["frame"] == 7
    a, b = q.compute(7), p.compute(7)
    assert not q.errors() and np.array_equal(a.properties["color"], b.properties["color"])
    assert is_session(str(proj / "s.atomview")) and not is_session(DUMP)


def test_moved_project_folder(tmp_path):
    proj, p = session(tmp_path)
    moved = tmp_path / "elsewhere" / "proj2"
    shutil.move(str(proj), str(moved))
    q = Pipeline.load(str(moved / "s.atomview"))
    assert os.path.samefile(q.source.path, moved / "data" / "run.dump")
    os.remove(moved / "data" / "run.dump")
    with pytest.raises(FileNotFoundError):
        Pipeline.load(str(moved / "s.atomview"))


def test_cli_renders_session_without_data_argument(tmp_path):
    proj, _ = session(tmp_path)
    r = subprocess.run([sys.executable, "-m", "atomview", "render", "-s", str(proj / "s.atomview"), "-o", str(tmp_path / "s.png"),
                        "--size", "200x150"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "s.png").stat().st_size > 1000


def test_gui_session_state(tmp_path):
    import pyvista as pv
    import atomview.render  # noqa: F401  (sets off-screen mode on headless Linux)
    if pv.OFF_SCREEN:
        pytest.skip("GUI test needs an interactive display")
    proj, _ = session(tmp_path)
    from PySide6.QtCore import QEvent, QUrl
    from PySide6.QtGui import QFileOpenEvent
    from PySide6.QtWidgets import QApplication
    from atomview.gui import main_window as mw
    app = QApplication.instance() or mw.App([])
    w = mw.MainWindow()
    app.window = w
    if isinstance(app, mw.App):  # the path a macOS double-click takes
        app.event(QFileOpenEvent(str(proj / "s.atomview")))
    else:
        w.open_any(str(proj / "s.atomview"))
    assert w.session_path == str(proj / "s.atomview") and w.fps.value() == 25 and w.frame_index == 7
    assert not w.dirty() and "s.atomview" in w.windowTitle()
    w.pipeline.add(m.ExpressionSelect(expression="z > 10"))
    assert w.dirty()
    assert w.save_session() and not w.dirty()
    assert [x["type"] for x in json.load(open(proj / "s.atomview"))["modifiers"]] == ["ColorCoding", "Slice", "ExpressionSelect"]
    w._saved_state = w._state()
    w.close()


def test_resave_elsewhere_updates_relative(tmp_path):
    proj, _ = session(tmp_path)
    (tmp_path / "other").mkdir()
    p = Pipeline.load(str(proj / "s.atomview"))
    p.save(str(tmp_path / "other" / "s.atomview"), p.session_extra)
    d = json.load(open(tmp_path / "other" / "s.atomview"))
    assert d["source_relative"] == os.path.join("..", "proj", "data", "run.dump") and d["frame"] == 7
