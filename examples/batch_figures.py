"""Make publication figures from a script, no GUI needed.

    .venv/bin/python examples/batch_figures.py
"""
import os

import atomview as av
from atomview import modifiers as m
from atomview.render import RenderSettings

HERE = os.path.dirname(os.path.abspath(__file__))
pipe = av.Pipeline(av.load(os.path.join(HERE, "data", "cu_void.dump")))
pipe.add(m.CommonNeighborAnalysis())
pipe.add(m.Slice(normal=(0, 1, 0), width=6))

style = RenderSettings(ambient_occlusion=True, parallel_projection=True, title="Cu with a void")
for frame in (0, 10, 19):
    av.render_image(pipe, f"cna_frame{frame}.png", frame=frame, size=(2400, 1800), settings=style, view="-y")
    data = pipe.compute(frame)
    print(frame, {k: v for k, v in data.attributes.items() if k.startswith("CNA.")})

av.render_animation(pipe, "cna.mp4", size=(1280, 960), fps=10, settings=style, view="-y")
