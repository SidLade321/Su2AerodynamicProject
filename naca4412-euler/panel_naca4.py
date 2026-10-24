"""Linear-strength vortex panel method (Kuethe & Chow) for any NACA 4-digit aerofoil.

Extends the symmetric version in ../naca0012-liftcurve/panel_method.py to cambered
sections, and also returns the surface velocity, so it gives the pressure distribution
and pitching moment as well as the lift. Results are inviscid and incompressible; a
Prandtl-Glauert factor 1/sqrt(1 - M^2) corrects them for low subsonic Mach number.
"""
import numpy as np


def naca4(code, n):
    """Closed-TE NACA 4-digit section, cosine-spaced, ordered TE(lower) -> LE -> TE(upper).

    Returns n panels per side (2n panels in total), with the TE point repeated at both
    ends as the panel method needs.
    """
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    b = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))
    yt = 5 * t * (0.2969 * np.sqrt(b) - 0.1260 * b - 0.3516 * b**2 + 0.2843 * b**3 - 0.1036 * b**4)
    if m == 0:
        yc = dyc = np.zeros_like(b)
    else:
        yc = np.where(b < p, m / p**2 * (2 * p * b - b**2), m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * b - b**2))
        dyc = np.where(b < p, 2 * m / p**2 * (p - b), 2 * m / (1 - p)**2 * (p - b))
    th = np.arctan(dyc)
    xu, yu = b - yt * np.sin(th), yc + yt * np.cos(th)
    xl, yl = b + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xl[::-1], xu[1:]]), np.concatenate([yl[::-1], yu[1:]])


def solve(x, y, alpha_deg):
    """Return a dict with CL (from circulation), Cm about c/4, and panel-midpoint Cp."""
    a = np.radians(alpha_deg)
    N = len(x) - 1
    xc, yc = 0.5 * (x[:-1] + x[1:]), 0.5 * (y[:-1] + y[1:])
    S = np.hypot(np.diff(x), np.diff(y))
    th = np.arctan2(np.diff(y), np.diff(x))

    dx = xc[:, None] - x[None, :-1]
    dy = yc[:, None] - y[None, :-1]
    ti, tj, Sj = th[:, None], th[None, :], S[None, :]
    A = -dx * np.cos(tj) - dy * np.sin(tj)
    B = dx**2 + dy**2
    C, D = np.sin(ti - tj), np.cos(ti - tj)
    E = dx * np.sin(tj) - dy * np.cos(tj)
    with np.errstate(divide="ignore", invalid="ignore"):
        F = np.log(1 + Sj * (Sj + 2 * A) / B)
        G = np.arctan2(E * Sj, B + A * Sj)
        P = dx * np.sin(ti - 2 * tj) + dy * np.cos(ti - 2 * tj)
        Q = dx * np.cos(ti - 2 * tj) - dy * np.sin(ti - 2 * tj)
        CN2 = D + 0.5 * Q * F / Sj - (A * C + D * E) * G / Sj
        CN1 = 0.5 * D * F + C * G - CN2
        CT2 = C + 0.5 * P * F / Sj + (A * D - C * E) * G / Sj
        CT1 = 0.5 * C * F - D * G - CT2
    ii = np.arange(N)
    CN1[ii, ii], CN2[ii, ii] = -1.0, 1.0
    CT1[ii, ii] = CT2[ii, ii] = 0.5 * np.pi

    def assemble(c1, c2):
        M = np.zeros((N, N + 1))
        M[:, 0] = c1[:, 0]
        M[:, N] = c2[:, N - 1]
        M[:, 1:N] = c1[:, 1:] + c2[:, :-1]
        return M

    AN, AT = assemble(CN1, CN2), assemble(CT1, CT2)
    M = np.vstack([AN, np.zeros(N + 1)])
    M[N, 0] = M[N, N] = 1.0                                     # Kutta condition
    rhs = np.concatenate([np.sin(th - a), [0.0]])
    g = np.linalg.solve(M, rhs)

    Vt = np.cos(th - a) + AT @ g
    cp = 1 - Vt**2
    CL = 2 * 2 * np.pi * np.sum(0.5 * (g[:-1] + g[1:]) * S)    # Kutta-Joukowski, c = V = 1

    # Pressure integration: force coefficient on each panel is -Cp * n * S, with n the
    # outward normal (the tangent rotated +90 deg for this clockwise ordering)
    nx, ny = -np.sin(th), np.cos(th)
    fx, fy = -cp * nx * S, -cp * ny * S
    CL_p = np.sum(fy * np.cos(a) - fx * np.sin(a))
    Cm = np.sum(fx * yc - fy * (xc - 0.25))                     # about c/4, nose-up positive
    return {"CL": CL, "CL_pressure": CL_p, "Cm": Cm, "x": xc, "y": yc, "cp": cp}


if __name__ == "__main__":
    for n in (50, 100, 200, 400):
        x, y = naca4("4412", n)
        r = solve(x, y, 5.0)
        print(f"{2 * n:4d} panels: CL(5 deg) = {r['CL']:.4f} (pressure-integrated {r['CL_pressure']:.4f}), "
              f"Cm_c/4 = {r['Cm']:+.4f}")
    # Sanity check against the symmetric-section code used for the NACA 0012 case
    import sys; sys.path.insert(0, "../naca0012-liftcurve")
    from panel_method import naca_symmetric, lift_coefficient
    x0, y0 = naca_symmetric(0.12, 100)
    print(f"NACA 0012 at 4 deg: this code {solve(*naca4('0012', 100), 4)['CL']:.5f}, "
          f"original {lift_coefficient(x0, y0, 4):.5f}")
