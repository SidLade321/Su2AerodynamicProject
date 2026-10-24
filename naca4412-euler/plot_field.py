"""Pressure-coefficient contours and the mesh near the aerofoil.

Reads the cell connectivity from the SU2 mesh file and the conserved variables from
SU2's CSV solution output (restart_flow.csv), so no VTK reader is needed. Cp is
recovered from the conserved variables in the same way as in analyse.py.

Usage:  python plot_field.py runs/medium_a5/restart_flow.csv mesh_naca4412_medium.su2
Writes cp_field.png.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm

flow_csv = sys.argv[1] if len(sys.argv) > 1 else "runs/medium_a5/restart_flow.csv"
mesh_file = sys.argv[2] if len(sys.argv) > 2 else "mesh_naca4412_medium.su2"

with open(mesh_file) as f:
    lines = f.read().splitlines()
start = next(i for i, l in enumerate(lines) if l.startswith("NELEM="))
nelem = int(lines[start].split("=")[1])
quads = np.array([[int(v) for v in l.split()[1:5]] for l in lines[start + 1:start + 1 + nelem]])

d = pd.read_csv(flow_csv)
d.columns = [c.strip().strip('"') for c in d.columns]
d = d.sort_values("PointID")
x, y = d["x"].to_numpy(), d["y"].to_numpy()
GAMMA, MACH = 1.4, 0.15
p = (GAMMA - 1) * (d["Energy"] - 0.5 * (d["Momentum_x"]**2 + d["Momentum_y"]**2) / d["Density"])
cp = ((p - 1.0) / (0.5 * GAMMA * MACH**2)).to_numpy()   # free-stream pressure = 1
tri = mtri.Triangulation(x, y, np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]]))

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1.45, 1]})
ax = axes[0]
levels = np.linspace(-2.0, 1.0, 31)
norm = TwoSlopeNorm(vmin=levels[0], vcenter=0.0, vmax=levels[-1])   # white = free-stream pressure
cs = ax.tricontourf(tri, np.clip(cp, levels[0], levels[-1]), levels=levels, cmap="RdBu", norm=norm)
ax.tricontour(tri, cp, levels=levels[::3], colors="k", linewidths=0.3)
fig.colorbar(cs, ax=ax, label="$C_p$", ticks=np.arange(-2.0, 1.01, 0.5), shrink=0.9)
ax.set_xlim(-0.3, 1.4); ax.set_ylim(-0.45, 0.55); ax.set_aspect("equal")
ax.set_title("Pressure coefficient, M = 0.15, α = 5°")
ax.set_xlabel("x / c"); ax.set_ylabel("y / c")

ax = axes[1]
xy = np.column_stack([x, y])
segs = np.concatenate([xy[quads[:, [i, (i + 1) % 4]]] for i in range(4)])
ax.add_collection(LineCollection(segs, colors="k", linewidths=0.2))
ax.set_xlim(-0.15, 1.15); ax.set_ylim(-0.45, 0.45); ax.set_aspect("equal")
ax.set_title(f"Mesh near the aerofoil ({len(quads):,} cells in total)")
ax.set_xlabel("x / c")
fig.tight_layout()
fig.savefig("cp_field.png", dpi=150)
print("Saved cp_field.png")
