"""AtomView — free, open-source visualization and analysis of atomistic data.

Scripting quick start:

    import atomview as av
    from atomview import modifiers as m

    pipe = av.Pipeline(av.load("dump.lammpstrj"))
    pipe.add(m.CommonNeighborAnalysis())
    pipe.add(lambda frame, data: data.attributes.update(N=data.n))
    data = pipe.compute(0)
    av.render_image(pipe, "figure.png", frame=0, size=(2000, 1500))
"""
from . import modifiers
from .data import Frame
from .io import load, write
from .pipeline import Modifier, Param, Pipeline

__version__ = "0.2.0"


def render_image(*args, **kwargs):
    from .render import render_image as _r
    return _r(*args, **kwargs)


def render_animation(*args, **kwargs):
    from .render import render_animation as _r
    return _r(*args, **kwargs)


__all__ = ["Frame", "Modifier", "Param", "Pipeline", "load", "write", "modifiers", "render_image", "render_animation"]
