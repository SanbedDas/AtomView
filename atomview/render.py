"""Scene building and publication-quality image/animation export (PyVista/VTK)."""
from __future__ import annotations

import os
import sys
from dataclasses import asdict, dataclass

import numpy as np

# Headless Linux (SSH session, cluster node, Docker): no X server or Wayland display.
# Use VTK's software OSMesa window so off-screen rendering works without one
# (needs `libosmesa6` on Ubuntu). Override with VTK_DEFAULT_OPENGL_WINDOW, e.g.
# vtkEGLRenderWindow for GPU rendering through EGL.
HEADLESS = sys.platform.startswith("linux") and not (os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
if HEADLESS:
    os.environ.setdefault("VTK_DEFAULT_OPENGL_WINDOW", "vtkOSOpenGLRenderWindow")

import pyvista as pv  # noqa: E402

if HEADLESS:
    pv.OFF_SCREEN = True
    try:  # PyVista warns about the missing X server even though OSMesa doesn't need one
        from pyvista.plotting import plotter as _pv_plotter
        _pv_plotter._warn_xserver.has_support = True
    except Exception:
        pass

from .data import TRANSPARENCY, Frame
from .styles import particle_colors, particle_radii, type_legend


@dataclass
class RenderSettings:
    background: tuple = (1.0, 1.0, 1.0)
    sphere_resolution: int = 16
    radius_scale: float = 1.0
    show_cell: bool = True
    cell_color: tuple = (0.15, 0.15, 0.15)
    cell_width: float = 2.0
    show_bonds: bool = True
    bond_radius: float = 0.0  # 0 = use the Create bonds setting
    show_legend: bool = True
    show_axes: bool = True
    selection_color: tuple = (1.0, 0.1, 0.1)
    ambient_occlusion: bool = False
    anti_aliasing: bool = True
    shading: str = "smooth"  # smooth | flat | pbr
    specular: float = 0.35
    parallel_projection: bool = False
    legend_font: int = 14
    title: str = ""
    overlay: str = ""  # per-frame text, e.g. "t = {Timestep*0.002:.1f} ps"
    overlay_position: str = "lower_right"
    render_mode: str = "auto"  # auto | mesh | imposter
    imposter_threshold: int = 30_000
    view_roll: float = 0.0  # degrees, rotates the default views in the image plane

    def to_dict(self):
        return asdict(self)

    @staticmethod
    def from_dict(d):
        s = RenderSettings()
        for k, v in (d or {}).items():
            if hasattr(s, k):
                setattr(s, k, tuple(v) if isinstance(getattr(s, k), tuple) else v)
        return s


def _text_color(bg):
    return (0.05, 0.05, 0.05) if np.mean(bg) > 0.5 else (0.95, 0.95, 0.95)


def cell_mesh(frame: Frame, dz=0.0) -> pv.PolyData:
    o, (a, b, c) = frame.origin, frame.cell
    if frame.attributes.get("Dimensions") == 2:  # cell outline in the particle plane, lifted by dz toward +z
        o = np.r_[o[:2], (frame.positions[0, 2] if frame.n else o[2]) + dz]
        return pv.PolyData(np.array([o, o + a, o + a + b, o + b]), lines=np.array([5, 0, 1, 2, 3, 0]))
    corners = np.array([o, o + a, o + b, o + c, o + a + b, o + a + c, o + b + c, o + a + b + c])
    edges = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]
    lines = np.hstack([[2, i, j] for i, j in edges])
    return pv.PolyData(corners, lines=lines)


def use_imposters(frame: Frame, s: RenderSettings) -> bool:
    return s.render_mode == "imposter" or (s.render_mode == "auto" and frame.n > s.imposter_threshold)


def particle_points(frame: Frame, s: RenderSettings) -> pv.PolyData | None:
    # One point per particle with rgba and radius, drawn as GPU sphere imposters.
    if frame.n == 0:
        return None
    colors = particle_colors(frame)
    colors[frame.selection] = s.selection_color
    alpha = 1.0 - np.clip(np.asarray(frame.properties.get(TRANSPARENCY, np.zeros(frame.n)), float), 0, 1)
    vis = alpha > 0.001
    pts = pv.PolyData(frame.positions[vis])
    pts["rgba"] = (np.clip(np.column_stack([colors[vis], alpha[vis]]), 0, 1) * 255).astype(np.uint8)
    pts["radius"] = particle_radii(frame)[vis] * s.radius_scale
    return pts


def particle_mesh(frame: Frame, s: RenderSettings) -> pv.PolyData | None:
    if frame.n == 0:
        return None
    colors = particle_colors(frame)
    sel = frame.selection
    colors[sel] = s.selection_color
    radii = particle_radii(frame) * s.radius_scale
    alpha = 1.0 - np.clip(np.asarray(frame.properties.get(TRANSPARENCY, np.zeros(frame.n)), float), 0, 1)
    visible = alpha > 0.001
    pts = pv.PolyData(frame.positions[visible])
    pts["radius"] = radii[visible]
    rgba = np.column_stack([colors[visible], alpha[visible]])
    pts["rgba"] = (np.clip(rgba, 0, 1) * 255).astype(np.uint8)
    res = max(4, int(s.sphere_resolution))
    sphere = pv.Sphere(radius=1.0, theta_resolution=res, phi_resolution=max(3, res // 2 + 2))
    glyphs = pts.glyph(geom=sphere, scale="radius", orient=False, factor=1.0)
    return glyphs


def bond_mesh(frame: Frame, s: RenderSettings) -> pv.PolyData | None:
    if not s.show_bonds or frame.bonds is None or len(frame.bonds) == 0:
        return None
    i, j = frame.bonds[:, 0], frame.bonds[:, 1]
    vec = frame.minimum_image(frame.positions[j] - frame.positions[i])
    colors = particle_colors(frame)
    # Each bond is drawn as two halves colored like their atoms; this also
    # handles bonds that cross a periodic boundary.
    starts = np.concatenate([frame.positions[i], frame.positions[j]])
    ends = np.concatenate([frame.positions[i] + vec / 2, frame.positions[j] - vec / 2])
    col = np.concatenate([colors[i], colors[j]])
    m = len(starts)
    pts = np.empty((2 * m, 3))
    pts[0::2], pts[1::2] = starts, ends
    lines = np.column_stack([np.full(m, 2), np.arange(0, 2 * m, 2), np.arange(1, 2 * m, 2)]).ravel()
    poly = pv.PolyData(pts, lines=lines)
    poly.cell_data["rgb"] = (col * 255).astype(np.uint8)
    radius = s.bond_radius or frame.attributes.get("_bond_radius", 0.15)
    return poly.tube(radius=radius, n_sides=max(6, s.sphere_resolution // 2))


def build_scene(plotter: pv.Plotter, frame: Frame, s: RenderSettings, interactive=False, font_scale=1.0):
    """Draw a frame into a plotter (replacing what a previous call drew)."""
    for name in ("particles", "bonds", "cell", "title", *[f"vectors{k}" for k in range(8)]):
        plotter.remove_actor(name, render=False)
    for actor in getattr(plotter, "_av_overlay", []):
        plotter.remove_actor(actor, render=False)
    plotter._av_overlay = []
    if getattr(plotter, "_av_scalarbar", None) is not None:
        plotter.remove_actor(plotter._av_scalarbar, render=False)
        plotter._av_scalarbar = None
    plotter.set_background(s.background)
    pbr = s.shading == "pbr"
    kw = dict(smooth_shading=s.shading != "flat", specular=0 if pbr else s.specular, specular_power=30,
              pbr=pbr, metallic=0.1 if pbr else None, roughness=0.45 if pbr else None)
    kw = {k: v for k, v in kw.items() if v is not None}
    imp = use_imposters(frame, s)
    mesh = particle_points(frame, s) if imp else particle_mesh(frame, s)
    if imp and mesh is not None and mesh.n_points:
        a = plotter.add_mesh(mesh, scalars="rgba", rgba=True, name="particles", render=False, show_scalar_bar=False,
                             style="points_gaussian", render_points_as_spheres=True, emissive=False,
                             specular=s.specular, specular_power=30)
        a.mapper.SetScaleArray("radius")
        a.mapper.SetScaleFactor(1.0)
    elif mesh is not None and mesh.n_points:
        transparent = TRANSPARENCY in frame.properties
        plotter.add_mesh(mesh, scalars="rgba", rgba=True, name="particles", render=False,
                         show_scalar_bar=False, **kw)
        if transparent:
            plotter.enable_depth_peeling(number_of_peels=8)
    tr = frame.attributes.get("_trails")
    if tr is not None and tr["path"].shape[0] > 1:
        _add_trails(plotter, frame, tr)
    else:
        plotter.remove_actor("trails", render=False)
    for k, vs in enumerate(frame.attributes.get("_vectors", [])):
        _add_vectors(plotter, frame, vs, f"vectors{k}", kw)
    bonds = bond_mesh(frame, s)
    if bonds is not None:
        plotter.add_mesh(bonds, scalars="rgb", rgb=True, name="bonds", render=False, show_scalar_bar=False, **kw)
    if s.show_cell:
        dz = 1.05 * particle_radii(frame).max() * s.radius_scale if frame.n else 0.0  # above the spheres in 2D
        plotter.add_mesh(cell_mesh(frame, dz), color=s.cell_color, line_width=s.cell_width, name="cell",
                         render=False, render_lines_as_tubes=True, lighting=False)
    text = _text_color(s.background)
    if s.show_legend:
        cb = frame.attributes.get("_colorbar")
        legend = frame.attributes.get("_legend") or (type_legend(frame) if cb is None else None)
        if cb:
            plotter._av_scalarbar = _add_colorbar(plotter, cb, text, int(s.legend_font * font_scale))
        elif legend and len(legend) <= 24:
            _add_legend(plotter, legend, text, int(s.legend_font * font_scale))
    if s.title:
        plotter.add_text(s.title, position="upper_left", font_size=int(s.legend_font * font_scale * 0.8),
                         color=text, name="title")
    if s.overlay:
        from .expressions import format_overlay
        plotter.add_text(format_overlay(s.overlay, frame, int(frame.attributes.get("SourceFrame", 0))),
                         position=s.overlay_position, font_size=int(s.legend_font * font_scale * 0.8), color=text, name="overlay")
    else:
        plotter.remove_actor("overlay", render=False)
    if s.parallel_projection:
        plotter.enable_parallel_projection()
    else:
        plotter.disable_parallel_projection()


def _add_vectors(plotter, frame: Frame, vs, name, kw, max_arrows=200_000):
    # Arrow glyphs (or plain line segments above max_arrows) for one Vectors modifier.
    from .modifiers.basic import vector_field
    if frame.n == 0: return
    v = vector_field(frame, vs["property"])
    mag = np.linalg.norm(v, axis=1)
    if vs.get("unit_length"): v = v / np.where(mag > 0, mag, 1)[:, None]
    v = v * vs["scale"]
    base = frame.positions - (v / 2 if vs.get("centered") else 0)
    ok = np.linalg.norm(v, axis=1) > 0
    base, v = base[ok], v[ok]
    rgb = (particle_colors(frame)[ok] if vs.get("color_by_particle") else np.tile(np.asarray(vs["color"], float)[:3], (len(v), 1)))
    rgb = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)
    if len(v) > max_arrows:
        n = len(v)
        pts = np.empty((2 * n, 3)); pts[0::2], pts[1::2] = base, base + v
        lines = np.column_stack([np.full(n, 2), np.arange(0, 2 * n, 2), np.arange(1, 2 * n, 2)]).ravel()
        poly = pv.PolyData(pts, lines=lines); poly.cell_data["rgb"] = rgb
        plotter.add_mesh(poly, scalars="rgb", rgb=True, name=name, render=False, show_scalar_bar=False, line_width=2)
        return
    pts = pv.PolyData(base); pts["vec"] = v; pts["len"] = np.linalg.norm(v, axis=1); pts["rgb"] = rgb
    w = max(vs.get("width", 0.1), 1e-3)
    arrow = pv.Arrow(start=(0, 0, 0), direction=(1, 0, 0), tip_length=0.3, tip_radius=min(2.2 * w, 0.45), shaft_radius=w,
                     tip_resolution=12, shaft_resolution=12)
    g = pts.glyph(orient="vec", scale="len", factor=1.0, geom=arrow)
    plotter.add_mesh(g, scalars="rgb", rgb=True, name=name, render=False, show_scalar_bar=False, **kw)


def _add_trails(plotter, frame: Frame, tr):
    # One polyline per particle through its previous positions; segments with a missing frame are skipped.
    P = tr["path"]  # (L, M, 3)
    L, M = P.shape[:2]
    ok = np.isfinite(P).all(-1)
    seg = ok[:-1] & ok[1:]  # (L-1, M)
    t, mi = np.nonzero(seg)
    if not len(t):
        return
    pts = np.concatenate([P[t, mi], P[t + 1, mi]])
    n = len(t)
    lines = np.column_stack([np.full(n, 2), np.arange(n), np.arange(n, 2 * n)]).ravel()
    poly = pv.PolyData(pts, lines=lines)
    rows = tr["rows"]
    col = particle_colors(frame)[rows] if rows.max(initial=-1) < frame.n else np.tile([0.4, 0.4, 0.4], (M, 1))
    fade = 0.35 + 0.65 * (t + 1) / (L - 1)  # older segments lighter
    rgb = col[mi] * fade[:, None] + (1 - fade[:, None]) * np.asarray(plotter.background_color.float_rgb)
    poly.cell_data["rgb"] = (np.clip(rgb, 0, 1) * 255).astype(np.uint8)
    plotter.add_mesh(poly, scalars="rgb", rgb=True, name="trails", render=False, show_scalar_bar=False,
                     line_width=tr["width"], lighting=False)


def _add_legend(plotter, entries, text_color, font):
    """A compact legend in the upper right: colored dot (with outline) + label in the text color."""
    import vtk
    W, H = plotter.window_size
    try:  # VTK scales text by DPI on high-resolution (Retina) screens
        dpi_scale = max(plotter.ren_win.GetDPI() / 72.0, 1.0)
    except Exception:
        dpi_scale = 1.0
    r = font * 0.42 * dpi_scale
    for k, (label, color) in enumerate(entries):
        y = 0.94 - k * 1.9 * font * dpi_scale / H
        for col, filled, rad in (((0.4, 0.4, 0.4), True, r + 1.2), (tuple(color), True, r)):
            disc = vtk.vtkRegularPolygonSource()
            disc.SetNumberOfSides(40)
            disc.SetRadius(rad)
            disc.SetCenter(0.80 * W, y * H, 0)
            mapper = vtk.vtkPolyDataMapper2D()
            mapper.SetInputConnection(disc.GetOutputPort())
            actor = vtk.vtkActor2D()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(col)
            plotter.add_actor(actor, reset_camera=False, render=False)
            plotter._av_overlay.append(actor)
        t = vtk.vtkTextActor()
        t.SetInput(str(label))
        tp = t.GetTextProperty()
        tp.SetColor(text_color)
        tp.SetFontSize(font)
        tp.SetVerticalJustificationToCentered()
        t.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        t.SetPosition(0.82, y)
        plotter.add_actor(t, reset_camera=False, render=False)
        plotter._av_overlay.append(t)


def _add_colorbar(plotter, cb, text_color, font):
    import matplotlib
    import vtk
    cmap = matplotlib.colormaps[cb["cmap"]]
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    lut.SetTableRange(cb["start"], cb["end"] if cb["end"] != cb["start"] else cb["start"] + 1)
    for k in range(256):
        r, g, b, _ = cmap(k / 255)
        lut.SetTableValue(k, r, g, b, 1.0)
    lut.Build()
    bar = vtk.vtkScalarBarActor()
    bar.SetLookupTable(lut)
    bar.SetTitle(cb["title"])
    bar.SetNumberOfLabels(5)
    bar.SetLabelFormat("%-#6.3g")
    bar.SetPosition(0.86, 0.2)
    bar.SetWidth(0.1)
    bar.SetHeight(0.6)
    bar.SetBarRatio(0.3)
    bar.SetVerticalTitleSeparation(int(font * 0.8))
    bar.SetUnconstrainedFontSize(True)
    for prop in (bar.GetTitleTextProperty(), bar.GetLabelTextProperty()):
        prop.SetColor(text_color)
        prop.SetFontSize(font)
        prop.SetBold(False)
        prop.SetItalic(False)
        prop.SetShadow(False)
    bar.GetTitleTextProperty().SetFontSize(int(font * 1.15))
    plotter.add_actor(bar, reset_camera=False, render=False)
    return bar


def apply_quality(plotter: pv.Plotter, s: RenderSettings):
    if s.anti_aliasing:
        try:
            plotter.enable_anti_aliasing("ssaa")
        except Exception:
            pass
    else:
        plotter.disable_anti_aliasing()
    if s.ambient_occlusion:
        plotter.enable_ssao(radius=2.0, bias=0.01, kernel_size=128, blur=True)
    else:
        try:
            plotter.disable_ssao()
        except Exception:
            pass


def default_camera(plotter, frame: Frame, view="iso", roll=0.0, rscale=1.0):
    if view == "iso" and frame.attributes.get("Dimensions") == 2:
        view = "+z"
    center = frame.origin + frame.cell.sum(axis=0) / 2
    if frame.attributes.get("Dimensions") == 2 and frame.n:
        center[2] = frame.positions[0, 2]
    diag = np.linalg.norm(frame.cell.sum(axis=0)) or 10.0
    dirs = {"iso": (1, -1.3, 0.9), "+x": (1, 0, 0), "-x": (-1, 0, 0), "+y": (0, 1, 0), "-y": (0, -1, 0),
            "+z": (0, 0, 1), "-z": (0, 0, -1)}
    d = np.asarray(dirs[view], dtype=float)
    d /= np.linalg.norm(d)
    up = (0, 1, 0) if view in ("+z", "-z") else (0, 0, 1)
    plotter.camera_position = [tuple(center + d * diag * 1.6), tuple(center), up]
    plotter.camera.Roll(roll)
    fit_camera(plotter, frame, rscale)


def fit_camera(plotter, frame: Frame, rscale=1.0, margin=1.04):
    # keep the view direction, center and zoom so the cell and particles just fill the window
    cam = plotter.camera
    pos, foc, up = (np.array(v, float) for v in plotter.camera_position)
    b = (pos - foc) / np.linalg.norm(pos - foc); u = up - (up @ b) * b; u /= np.linalg.norm(u); r = np.cross(u, b)
    P = frame.origin + np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float) @ frame.cell
    if frame.n: P = np.vstack([P, frame.positions])
    P = P - foc; x, y, z = P @ r, P @ u, P @ b
    rad = particle_radii(frame).max() * rscale if frame.n else 0.0
    cx, cy = (x.max() + x.min()) / 2, (y.max() + y.min()) / 2
    x, y, foc = np.abs(x - cx) + rad, np.abs(y - cy) + rad, foc + cx * r + cy * u
    W, H = plotter.window_size; ty = np.tan(np.radians(cam.view_angle) / 2); tx = ty * W / H
    cam.parallel_scale = margin * max(y.max(), x.max() * H / W)
    d = margin * max((z + x / tx).max(), (z + y / ty).max())
    plotter.camera_position = [tuple(foc + b * d), tuple(foc), tuple(u)]
    plotter.reset_camera_clipping_range()


def _offscreen(size, s: RenderSettings):
    p = pv.Plotter(off_screen=True, window_size=list(size))
    apply_quality(p, s)
    return p


def render_image(pipeline, path, frame=0, size=(1600, 1200), settings: RenderSettings | None = None,
                 camera=None, transparent=False, view="iso"):
    """Render one frame of a pipeline to a PNG/JPG/TIFF (or SVG/PDF/EPS as vector)."""
    s = settings or RenderSettings()
    data = pipeline.compute(frame) if hasattr(pipeline, "compute") else pipeline
    p = _offscreen(size, s)
    build_scene(p, data, s, font_scale=size[1] / 900)
    if camera is not None:
        p.camera_position = [tuple(c) for c in camera]
        p.reset_camera_clipping_range()
    else:
        default_camera(p, data, view, s.view_roll, s.radius_scale)
    if s.show_axes:
        p.add_axes(color=_text_color(s.background))
    if str(path).lower().endswith((".svg", ".pdf", ".eps", ".ps")):
        p.save_graphic(str(path))
    else:
        p.screenshot(str(path), transparent_background=transparent)
    p.close()
    return path


def render_animation(pipeline, path, frames=None, size=(1280, 960), fps=10, settings: RenderSettings | None = None,
                     camera=None, view="iso", progress=None):
    """Render frames to MP4/GIF/AVI. `progress(i, n)` may return False to cancel."""
    import imageio.v2 as imageio
    s = settings or RenderSettings()
    frames = list(range(pipeline.n_frames)) if frames is None else list(frames)
    p = _offscreen(size, s)
    first = pipeline.compute(frames[0])
    build_scene(p, first, s, font_scale=size[1] / 900)
    if camera is not None:
        p.camera_position = [tuple(c) for c in camera]
    else:
        default_camera(p, first, view, s.view_roll, s.radius_scale)
    if s.show_axes:
        p.add_axes(color=_text_color(s.background))
    gif = str(path).lower().endswith(".gif")
    writer = imageio.get_writer(str(path), **({"duration": 1000 / fps, "loop": 0} if gif else {"fps": fps, "quality": 8, "macro_block_size": 1}))
    try:
        for k, f in enumerate(frames):
            data = pipeline.compute(f)
            cam = p.camera_position
            build_scene(p, data, s, font_scale=size[1] / 900)
            p.camera_position = cam
            p.render()
            writer.append_data(p.screenshot(return_img=True))
            if progress and progress(k + 1, len(frames)) is False:
                break
    finally:
        writer.close()
        p.close()
    return path
