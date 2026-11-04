# Python script modifier example: load this file in a "Python script" modifier.
# Selects the top 10 % of atoms and stores their mean energy as an attribute.
import numpy as np


def modify(frame, data):
    z = data.positions[:, 2]
    top = z > np.percentile(z, 90)
    data['selection'] = top
    if 'c_pe' in data:
        data.attributes['TopLayerPE'] = float(data['c_pe'][top].mean())
    print(f"frame {frame}: {top.sum()} atoms in the top layer")
