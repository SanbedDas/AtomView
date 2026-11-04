# App icon: three flat particles on a rounded tile -> icon.png, .icns, .ico
import os
import numpy as np
from PIL import Image, ImageDraw

here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "..", "atomview", "resources")
K = 4; S = 1024 * K; c = S / 2; r = 150 * K; bg = (28, 38, 66)
im = Image.new("RGBA", (S, S), (0, 0, 0, 0)); d = ImageDraw.Draw(im)
d.rounded_rectangle([100 * K, 100 * K, S - 100 * K, S - 100 * K], 185 * K, fill=bg)
for k, col in enumerate(((245, 120, 50), (55, 150, 235), (250, 200, 70))):
    a = np.pi / 2 + k * 2 * np.pi / 3; x, y = c + 1.05 * r * np.cos(a), c + 1.05 * r * np.sin(a) + 20 * K
    g = r + 22 * K; d.ellipse([x - g, y - g, x + g, y + g], fill=bg); d.ellipse([x - r, y - r, x + r, y + r], fill=col)
tile = im.resize((1024, 1024), Image.LANCZOS)
tile.save(os.path.join(out, "icon.png"))
tile.save(os.path.join(here, "icon.icns"))
tile.save(os.path.join(here, "icon.ico"), sizes=[(s, s) for s in (16, 24, 32, 48, 64, 128, 256)])
print("icon written")
