"""Structured O-mesh generator for a NACA 4-digit aerofoil, written in SU2 native format.

The mesh is built by marching outward from the aerofoil surface one layer at a time
(a simple explicit version of hyperbolic grid generation):

  * each new layer moves along the local surface normal, by a step that grows
    geometrically from a small wall spacing out to a far field ~100 chords away;
  * the marching directions are smoothed between neighbours, which fans the grid
    out cleanly around the sharp trailing edge instead of letting lines cross;
  * each layer is lightly smoothed along its length, more strongly far from the wall.

A fine mesh is generated first. The medium and coarse meshes are taken from it by
keeping every 2nd and every 4th grid line in both directions, so the three meshes
form a consistent family with a refinement ratio of exactly 2, which is what the
Richardson extrapolation in mesh_convergence.py assumes.

Usage:  python make_mesh.py            -> mesh_naca4412_{coarse,medium,fine}.su2
"""
import numpy as np

# ---- Parameters (fine level) -------------------------------------------------
NACA = "4412"
N_SIDE = 256          # surface intervals per side (upper and lower) on the fine mesh
N_LAYERS = 192        # layers in the wall-normal direction on the fine mesh
WALL_SPACING = 2.5e-4 # first-layer height, chords (fine mesh)
FARFIELD = 100.0      # approximate far-field distance, chords


def naca4(code, n):
    """Closed-trailing-edge NACA 4-digit section, cosine-spaced (n intervals per side).

    Returned as one closed loop ordered TE -> lower surface -> LE -> upper surface,
    without repeating the trailing-edge point. That ordering runs clockwise, so the
    outward normal is the tangent rotated by +90 degrees.
    """
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    b = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))         # x/c, LE -> TE
    yt = 5 * t * (0.2969 * np.sqrt(b) - 0.1260 * b - 0.3516 * b**2
                  + 0.2843 * b**3 - 0.1036 * b**4)                  # -0.1036 closes the TE
    yc = np.where(b < p, m / p**2 * (2 * p * b - b**2),
                  m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * b - b**2))
    dyc = np.where(b < p, 2 * m / p**2 * (p - b), 2 * m / (1 - p)**2 * (p - b))
    th = np.arctan(dyc)
    xu, yu = b - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = b + yt * np.sin(th), yc - yt * np.cos(th)
    x = np.concatenate([xl[::-1], xu[1:-1]])   # TE(lower) ... LE ... (upper), TE not repeated
    y = np.concatenate([yl[::-1], yu[1:-1]])
    x[0], y[0] = 1.0, 0.0                      # exact closed trailing edge
    return np.column_stack([x, y])


def layer_steps(ds0, total, n):
    """Geometric step sizes: first step ds0, n steps summing to `total`."""
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(200):                       # bisection for the growth ratio
        r = 0.5 * (lo + hi)
        s = ds0 * (r**n - 1) / (r - 1)
        lo, hi = (r, hi) if s < total else (lo, r)
    return ds0 * r ** np.arange(n), r


def normals(P):
    """Outward unit normals at the points of a closed clockwise loop (edge-normal average)."""
    e = np.roll(P, -1, axis=0) - P                         # edge i -> i+1
    en = np.column_stack([-e[:, 1], e[:, 0]])
    en /= np.linalg.norm(en, axis=1)[:, None]
    n = en + np.roll(en, 1, axis=0)                        # edges either side of each point
    return n / np.linalg.norm(n, axis=1)[:, None]


def smooth_loop(A, nu, passes):
    for _ in range(passes):
        A = A + nu * (0.5 * (np.roll(A, 1, axis=0) + np.roll(A, -1, axis=0)) - A)
    return A


def march(surface, steps, blend_dist=2.0, centre=(0.5, 0.0)):
    """Normal marching for the near field, then straight-line blending to a circle.

    Explicit marching becomes unstable once the step is much larger than the point
    spacing along a layer, so it is only used out to `blend_dist` chords. From there
    each grid line runs straight out to a circular far field along the direction from
    the mesh centre, keeping the same geometric spacing.
    """
    grid = [surface]
    P = surface.copy()
    dist = np.cumsum(steps)
    k0 = int(np.searchsorted(dist, blend_dist))
    for k in range(k0 + 1):
        frac = (k + 1) / (k0 + 1)
        n = smooth_loop(normals(P), 0.5, 4 + int(30 * frac))   # direction smoothing
        n /= np.linalg.norm(n, axis=1)[:, None]
        P = P + steps[k] * n
        P = smooth_loop(P, 0.02 + 0.3 * frac, 2)                # along-layer smoothing
        grid.append(P.copy())

    c = np.asarray(centre)
    radial = (P - c) / np.linalg.norm(P - c, axis=1)[:, None]
    R = dist[-1] + 0.5                                          # far-field radius about the centre
    outer = c + R * radial
    rest = steps[k0 + 1:]
    g = np.cumsum(rest) / rest.sum()
    for gk in g:
        grid.append(P + gk * (outer - P))
    return np.array(grid)                                       # shape (K+1, N, 2)


def signed_areas(G):
    """Signed area of every quad, using the same node ordering as write_su2 (must be > 0)."""
    a, a2 = G[:-1], np.roll(G[:-1], -1, axis=1)
    b2, b = np.roll(G[1:], -1, axis=1), G[1:]
    d1, d2 = b2 - a, b - a2
    return 0.5 * (d1[..., 0] * d2[..., 1] - d1[..., 1] * d2[..., 0])


def skew_deg(G):
    """Largest deviation from 90 degrees of the corner angles, per cell."""
    a, a2 = G[:-1], np.roll(G[:-1], -1, axis=1)
    b2, b = np.roll(G[1:], -1, axis=1), G[1:]
    worst = np.zeros(a.shape[:2])
    for p, q, r in ((a, a2, b), (a2, b2, a), (b2, b, a2), (b, a, b2)):
        u, v = q - p, r - p
        cosang = (u * v).sum(-1) / (np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1))
        worst = np.maximum(worst, np.abs(90 - np.degrees(np.arccos(np.clip(cosang, -1, 1)))))
    return worst


def write_su2(G, fname):
    K1, N, _ = G.shape
    pts = G.reshape(-1, 2)
    idx = lambda k, i: k * N + (i % N)
    with open(fname, "w") as f:
        f.write(f"% NACA {NACA} O-mesh, {N} x {K1 - 1} quads, generated by make_mesh.py\n")
        f.write("NDIME= 2\n")
        f.write(f"NELEM= {N * (K1 - 1)}\n")
        e = 0
        for k in range(K1 - 1):
            for i in range(N):
                # VTK quad = 9, counter-clockwise: outward-layer orientation gives positive area
                f.write(f"9 {idx(k, i)} {idx(k, i + 1)} {idx(k + 1, i + 1)} {idx(k + 1, i)} {e}\n")
                e += 1
        f.write(f"NPOIN= {len(pts)}\n")
        for j, (x, y) in enumerate(pts):
            f.write(f"{x:.12e} {y:.12e} {j}\n")
        f.write("NMARK= 2\n")
        # Boundary lines (VTK type 3) are ordered with the fluid on their left, as SU2
        # expects: clockwise round the aerofoil, anticlockwise round the far field.
        f.write(f"MARKER_TAG= airfoil\nMARKER_ELEMS= {N}\n")
        for i in range(N):
            f.write(f"3 {idx(0, i)} {idx(0, i + 1)}\n")
        f.write(f"MARKER_TAG= farfield\nMARKER_ELEMS= {N}\n")
        for i in range(N):
            f.write(f"3 {idx(K1 - 1, i + 1)} {idx(K1 - 1, i)}\n")


if __name__ == "__main__":
    surf = naca4(NACA, N_SIDE)
    steps, ratio = layer_steps(WALL_SPACING, FARFIELD, N_LAYERS)
    fine = march(surf, steps)
    r_ff = np.linalg.norm(fine[-1] - [0.5, 0], axis=1)
    print(f"Fine mesh: {surf.shape[0]} x {N_LAYERS}, growth ratio {ratio:.4f}, "
          f"far-field radius {r_ff.mean():.1f} chords")

    for name, s in (("fine", 1), ("medium", 2), ("coarse", 4)):
        G = fine[::s, ::s]
        sa, sk = signed_areas(G), skew_deg(G)
        if (sa <= 0).any():
            raise SystemExit(f"{name}: {(sa <= 0).sum()} inverted cells")
        write_su2(G, f"mesh_naca4412_{name}.su2")
        print(f"  {name:6s}: {G.shape[1]:4d} x {G.shape[0] - 1:3d} = {sa.size:6d} cells, "
              f"wall spacing {np.linalg.norm(G[1] - G[0], axis=1).min():.2e} c, "
              f"max skew {sk.max():.1f} deg (99th pct {np.percentile(sk, 99):.1f})")
