# PyInstaller build: one folder with AtomView (GUI) and atomview-cli (console) sharing the same libraries.
#   pyinstaller packaging/atomview.spec --noconfirm
# Qt stays as separate shared libraries in the folder (LGPL: users can replace them).
import sys
from PyInstaller.utils.hooks import collect_all, collect_submodules

datas, binaries, hidden = [], [], []
for pkg in ("atomview", "pyvista", "pyvistaqt", "freud", "parsnip", "chemfiles", "imageio_ffmpeg", "h5py", "ase"):
    d, b, h = collect_all(pkg)
    datas += d; binaries += b; hidden += h
hidden += collect_submodules("ase.io") + collect_submodules("vtkmodules") + ["qtpy", "atomview._ptmlib"]
datas += [("../THIRD_PARTY_NOTICES.md", "."), ("../LICENSE", ".")]
excludes = ["tkinter", "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtWebEngineQuick",
            "PySide6.Qt3DCore", "PySide6.Qt3DRender", "PySide6.QtQuick3D", "PySide6.QtMultimedia", "PySide6.QtPdf",
            "PySide6.QtCharts", "PySide6.QtDataVisualization", "PySide6.QtGraphs", "PySide6.QtLocation", "IPython", "pytest"]

a = Analysis(["entry.py"], pathex=[".."], binaries=binaries, datas=datas, hiddenimports=hidden, excludes=excludes)
pyz = PYZ(a.pure)
icon = "icon.ico" if sys.platform == "win32" else "icon.icns"
gui = EXE(pyz, a.scripts, [], exclude_binaries=True, name="AtomView", console=False, icon=icon)
cli = EXE(pyz, a.scripts, [], exclude_binaries=True, name="atomview-cli", console=True, icon=icon)
coll = COLLECT(gui, cli, a.binaries, a.datas, name="AtomView")
if sys.platform == "darwin":
    # Double-clicking a .atomview session opens AtomView; data files list AtomView under Open With
    # (rank Alternate, so it does not take over .mat or .pdb from their usual apps).
    data_ext = ["xyz", "extxyz", "dump", "lammpstrj", "lmp", "mat", "pdb", "gro", "xtc", "trr", "dcd", "cif", "vasp", "data"]
    app = BUNDLE(coll, name="AtomView.app", bundle_identifier="org.atomview.AtomView", icon="icon.icns",
                 info_plist={"NSHighResolutionCapable": True, "CFBundleShortVersionString": "0.2.0",
                             "CFBundleDocumentTypes": [
                                 {"CFBundleTypeName": "AtomView session", "CFBundleTypeRole": "Editor",
                                  "LSHandlerRank": "Owner", "LSItemContentTypes": ["org.atomview.session"],
                                  "CFBundleTypeIconFile": "icon.icns"},
                                 {"CFBundleTypeName": "Particle data", "CFBundleTypeRole": "Viewer",
                                  "LSHandlerRank": "Alternate", "CFBundleTypeExtensions": data_ext}],
                             "UTExportedTypeDeclarations": [
                                 {"UTTypeIdentifier": "org.atomview.session", "UTTypeDescription": "AtomView session",
                                  "UTTypeConformsTo": ["public.json", "public.data"], "UTTypeIconFile": "icon.icns",
                                  "UTTypeTagSpecification": {"public.filename-extension": ["atomview"]}}]})
