"""Third-party modifiers: entry points (group 'atomview.modifiers') and plugin folders.

A plugin is any module that defines Modifier subclasses with a `title`; defining the class registers it.
Folders searched: ~/.atomview/plugins and every folder in ATOMVIEW_PLUGIN_PATH (os.pathsep separated).
"""
from __future__ import annotations

import glob
import importlib.util
import os
import traceback
from importlib.metadata import entry_points

LOADED: list[str] = []
ERRORS: list[tuple[str, str]] = []
_done = False


def plugin_dirs():
    extra = [d for d in os.environ.get("ATOMVIEW_PLUGIN_PATH", "").split(os.pathsep) if d]
    return [os.path.join(os.path.expanduser("~"), ".atomview", "plugins"), *extra]


def load_plugins(force=False):
    # Imports every plugin once; a failing plugin is recorded in ERRORS instead of stopping the app.
    global _done
    if _done and not force:
        return LOADED
    _done = True
    from . import modifiers  # noqa: F401  (built-ins first, so plugins can subclass them)
    for ep in entry_points(group="atomview.modifiers"):
        try:
            ep.load()
            LOADED.append(f"{ep.name} ({ep.value})")
        except Exception:
            ERRORS.append((ep.name, traceback.format_exc(limit=3)))
    for d in plugin_dirs():
        for fn in sorted(glob.glob(os.path.join(d, "*.py"))):
            name = "atomview_plugin_" + os.path.splitext(os.path.basename(fn))[0]
            try:
                spec = importlib.util.spec_from_file_location(name, fn)
                spec.loader.exec_module(importlib.util.module_from_spec(spec))
                LOADED.append(fn)
            except Exception:
                ERRORS.append((fn, traceback.format_exc(limit=3)))
    return LOADED
