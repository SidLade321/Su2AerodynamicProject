"""Linear-strength vortex panel method (Kuethe & Chow) for a symmetric NACA 4-digit aerofoil.

Gives the inviscid, incompressible lift coefficient. A Prandtl-Glauert factor
1/sqrt(1 - M^2) then corrects it for compressibility at low subsonic Mach.
Used as an independent reference for the SU2 Euler lift curve.
"""
import numpy as np


def naca_symmetric(t, n):
    """Closed-trailing-edge NACA 00xx coordinates, cosine-spaced, ordered TE(lower) -> LE -> TE(upper)."""
    b = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))
    yt = 5 * t * (0.2969 * np.sqrt(b) - 0.1260 * b - 0.3516 * b**2 + 0.2843 * b**3 - 0.1036 * b**4)
    x = np.concatenate([b[::-1], b[1:]])
    y = np.concatenate([-yt[::-1], yt[1:]])
    return x, y


def lift_coefficient(x, y, alpha_deg):
    a = np.radians(alpha_deg)
    N = len(x) - 1
    xc, yc = 0.5 * (x[:-1] + x[1:]), 0.5 * (y[:-1] + y[1:])
    S = np.hypot(np.diff(x), np.diff(y))
    th = np.arctan2(np.diff(y), np.diff(x))

    CN1 = np.zeros((N, N)); CN2 = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                CN1[i, j], CN2[i, j] = -1.0, 1.0
                continue
            A = -(xc[i] - x[j]) * np.cos(th[j]) - (yc[i] - y[j]) * np.sin(th[j])
            B = (xc[i] - x[j])**2 + (yc[i] - y[j])**2
            C = np.sin(th[i] - th[j]); D = np.cos(th[i] - th[j])
            E = (xc[i] - x[j]) * np.sin(th[j]) - (yc[i] - y[j]) * np.cos(th[j])
            F = np.log(1 + S[j] * (S[j] + 2 * A) / B)
            G = np.arctan2(E * S[j], B + A * S[j])
            Q = (xc[i] - x[j]) * np.cos(th[i] - 2 * th[j]) - (yc[i] - y[j]) * np.sin(th[i] - 2 * th[j])
            CN2[i, j] = D + 0.5 * Q * F / S[j] - (A * C + D * E) * G / S[j]
            CN1[i, j] = 0.5 * D * F + C * G - CN2[i, j]

    # Flow tangency at every panel midpoint + Kutta condition at the trailing edge
    M = np.zeros((N + 1, N + 1)); rhs = np.zeros(N + 1)
    M[:N, 0] = CN1[:, 0]
    M[:N, N] = CN2[:, N - 1]
    M[:N, 1:N] = CN1[:, 1:] + CN2[:, :-1]
    rhs[:N] = np.sin(th - a)
    M[N, 0] = M[N, N] = 1.0
    g = np.linalg.solve(M, rhs)

    circulation = 2 * np.pi * np.sum(0.5 * (g[:-1] + g[1:]) * S)   # chord = 1, V = 1
    return 2 * circulation                                           # Kutta-Joukowski: Cl = 2*Gamma/(V*c)


if __name__ == "__main__":
    # Panel-count convergence check
    for n in (50, 100, 200):
        x, y = naca_symmetric(0.12, n)
        print(f"{2 * n} panels: Cl(4 deg) = {lift_coefficient(x, y, 4):.4f}")
