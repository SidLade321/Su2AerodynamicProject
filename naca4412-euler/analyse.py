"""NACA 4412 Euler results: mesh convergence, pressure distribution and lift curve.

Reads the SU2 histories and surface files that run_cases.sh copies into results/,
compares them with the vortex panel method (panel_naca4.py, Prandtl-Glauert corrected)
and writes cp_validation.png, liftcurve_naca4412.png and mesh_convergence.png.
"""
import glob
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from panel_naca4 import naca4, solve

MACH, GAMMA, AOA = 0.15, 1.4, 5.0
BETA = np.sqrt(1 - MACH**2)                     # Prandtl-Glauert factor
SU2_BLUE = "#1f77b4"
LEVELS = {"coarse": 128 * 48, "medium": 256 * 96, "fine": 512 * 192}


def history(tag):
    h = pd.read_csv(f"results/history_{tag}.csv")
    h.columns = [c.strip().strip('"') for c in h.columns]
    return {"CL": h["CL"].iloc[-1], "CD": h["CD"].iloc[-1],
            "Cm": -h["CMz"].iloc[-1],           # SU2's CMz is positive nose-down
            "iters": int(h["Inner_Iter"].iloc[-1]) + 1,
            "res_drop": h["rms[Rho]"].iloc[0] - h["rms[Rho]"].iloc[-1],
            "res_final": h["rms[Rho]"].iloc[-1]}


def surface_cp(tag):
    """Cp round the aerofoil, in mesh order (TE -> lower -> LE -> upper -> TE).

    SU2's surface CSV holds the conserved variables, so pressure is recovered from the
    ideal-gas relation. Free-stream pressure is 1 (FREESTREAM_PRESS_EQ_ONE), and the
    dynamic pressure is 0.5*gamma*p*M^2 whatever the non-dimensionalisation.
    """
    s = pd.read_csv(f"results/surface_{tag}.csv")
    s.columns = [c.strip().strip('"') for c in s.columns]
    s = s.sort_values("PointID").reset_index(drop=True)
    rho = s["Density"]
    p = (GAMMA - 1) * (s["Energy"] - 0.5 * (s["Momentum_x"]**2 + s["Momentum_y"]**2) / rho)
    s["Cp"] = (p - 1.0) / (0.5 * GAMMA * MACH**2)
    return s


def split_surfaces(x, y, cp):
    """Walk the closed loop from the leading edge: lower runs TE->LE, upper LE->TE.

    Splitting by the sign of y fails on a cambered section (the lower surface is above
    y = 0 near the leading edge), so the split is made at the minimum-x point instead.
    """
    le = int(np.argmin(x))
    lower = slice(0, le + 1)
    upper = np.r_[le:len(x), 0]                 # close back onto the trailing edge
    return (x[upper], y[upper], cp[upper]), (x[lower], y[lower], cp[lower])


def integrate_cl(s, alpha_deg):
    """Lift from the extracted surface Cp (check on the extraction, cf. the flat-plate case)."""
    x, y, cp = s["x"].to_numpy(), s["y"].to_numpy(), s["Cp"].to_numpy()
    x2, y2, cp2 = np.roll(x, -1), np.roll(y, -1), np.roll(cp, -1)
    dx, dy, cpm = x2 - x, y2 - y, 0.5 * (cp + cp2)
    fx, fy = cpm * dy, -cpm * dx                # -Cp * outward normal * length (clockwise loop)
    a = np.radians(alpha_deg)
    return np.sum(fy * np.cos(a) - fx * np.sin(a))


def richardson(f1, f2, f3, r=2.0):
    """Observed order, extrapolated value and fine-grid GCI (Roache, safety factor 1.25)."""
    p = np.log(abs((f3 - f2) / (f2 - f1))) / np.log(r)
    f_ext = f1 + (f1 - f2) / (r**p - 1)
    gci = 1.25 * abs((f1 - f2) / f1) / (r**p - 1)
    return p, f_ext, gci


def thin_aerofoil_zero_lift(code="4412"):
    m, p = int(code[0]) / 100, int(code[1]) / 10
    t0 = np.linspace(0, np.pi, 20001)
    x = 0.5 * (1 - np.cos(t0))
    dyc = np.where(x < p, 2 * m / p**2 * (p - x), 2 * m / (1 - p)**2 * (p - x))
    return -np.trapezoid(dyc * (np.cos(t0) - 1), t0) / np.pi          # radians


def style(ax):
    ax.grid(alpha=0.3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


if __name__ == "__main__":
    pd.set_option("display.float_format", "{:.5f}".format)

    # ---------- 1. Mesh convergence at AoA = 5 deg ----------
    xr, yr = naca4("4412", 400)
    ref = solve(xr, yr, AOA)
    CL_ref, Cm_ref = ref["CL"] / BETA, ref["Cm"] / BETA

    rows = []
    for lvl, cells in LEVELS.items():
        h = history(f"{lvl}_a5")
        h["CL_from_Cp"] = integrate_cl(surface_cp(f"{lvl}_a5"), AOA)
        rows.append({"mesh": lvl, "cells": cells, **h})
    mc = pd.DataFrame(rows)
    mc["CL_vs_panel_%"] = 100 * (mc["CL"] - CL_ref) / CL_ref
    print("Mesh convergence, AoA = 5 deg")
    print(mc.to_string(index=False))

    f = {k: mc.set_index("mesh")[k] for k in ("CL", "CD", "Cm")}
    print(f"\nPanel method + PG (800 panels): CL = {CL_ref:.4f}, Cm_c/4 = {Cm_ref:+.4f}")
    p, fe, gci = richardson(f["CL"]["fine"], f["CL"]["medium"], f["CL"]["coarse"])
    print(f"CL: observed order {p:.2f}, extrapolated {fe:.4f} "
          f"({100 * (fe - CL_ref) / CL_ref:+.3f}% vs panel), fine-mesh GCI {100 * gci:.3f}%")
    # Cm changes by under 0.1% and not monotonically, so Richardson extrapolation does not
    # apply to it (the observed order would be meaningless); report the spread instead.
    cm = f["Cm"]
    print(f"Cm: {cm['coarse']:.5f}, {cm['medium']:.5f}, {cm['fine']:.5f} "
          f"(spread {100 * (cm.max() - cm.min()) / abs(cm.mean()):.2f}%), "
          f"fine vs panel {100 * (cm['fine'] - Cm_ref) / abs(Cm_ref):+.2f}%")
    cd = f["CD"]
    print(f"CD (exactly 0 in inviscid, shock-free flow), in drag counts: "
          f"{1e4 * cd['coarse']:.1f} -> {1e4 * cd['medium']:.1f} -> {1e4 * cd['fine']:.1f}")

    # ---------- 2. Pressure distribution ----------
    s = surface_cp("fine_a5")
    (xu, yu, cpu), (xl, yl, cpl) = split_surfaces(s["x"].to_numpy(), s["y"].to_numpy(), s["Cp"].to_numpy())
    xp, cpp = ref["x"], ref["cp"] / BETA
    le = int(np.argmin(xp))

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(xu, cpu, "-", color=SU2_BLUE, lw=3.5, alpha=0.6, label="SU2 Euler, fine mesh")
    ax.plot(xl, cpl, "-", color=SU2_BLUE, lw=3.5, alpha=0.6)
    ax.plot(xp[le:], cpp[le:], "k--", lw=1.2, label="Vortex panel method + Prandtl–Glauert")
    ax.plot(xp[:le + 1], cpp[:le + 1], "k--", lw=1.2)
    ax.annotate("upper surface", xy=(0.6, np.interp(0.6, xu, cpu)), xytext=(0.62, -1.25),
                arrowprops=dict(arrowstyle="-", lw=0.6), fontsize=9)
    ax.annotate("lower surface", xy=(0.6, np.interp(0.6, xl[::-1], cpl[::-1])), xytext=(0.62, 0.75),
                arrowprops=dict(arrowstyle="-", lw=0.6), fontsize=9)
    ax.annotate("kink at x/c = 0.4: the NACA 4-digit\ncamber line changes formula here",
                xy=(0.4, np.interp(0.4, xu, cpu)), xytext=(0.12, -0.55),
                arrowprops=dict(arrowstyle="-", lw=0.6), fontsize=8, color="0.3")
    ax.set_ylim(1.25, -2.0)          # (inverted) clips the two cells beside the sharp TE, Cp ~ 1.8
    ax.set_xlabel("x / c")
    ax.set_ylabel("Pressure coefficient $C_p$")
    ax.set_title(f"NACA 4412 pressure distribution, M = {MACH}, α = {AOA:.0f}°")
    ax.legend(loc="lower center", bbox_to_anchor=(0.38, 0.02))
    style(ax)
    fig.tight_layout()
    fig.savefig("cp_validation.png", dpi=150)

    # ---------- 3. Lift curve (medium mesh) ----------
    tags = sorted(glob.glob("results/history_medium_a*.csv"),
                  key=lambda t: float(re.search(r"_a(-?\d+)", t).group(1)))
    lc = []
    for t in tags:
        a = float(re.search(r"_a(-?\d+)", t).group(1))
        h = history(f"medium_a{a:g}")
        lc.append({"AoA": a, "CL_SU2": h["CL"], "Cm_SU2": h["Cm"], "res_drop": h["res_drop"],
                   "CL_panel_PG": solve(xr, yr, a)["CL"] / BETA})
    lc = pd.DataFrame(lc)
    lc["Diff"] = lc["CL_SU2"] - lc["CL_panel_PG"]          # absolute: % is misleading near CL = 0
    lc["Diff_%"] = 100 * lc["Diff"] / lc["CL_panel_PG"]
    print("\nLift curve, medium mesh")
    print(lc.to_string(index=False))

    a_rad = np.radians(lc["AoA"])
    (k_su2, c_su2), (k_ref, c_ref) = np.polyfit(a_rad, lc["CL_SU2"], 1), np.polyfit(a_rad, lc["CL_panel_PG"], 1)
    a0_thin = thin_aerofoil_zero_lift()
    print(f"Lift-curve slope: SU2 {k_su2:.3f} /rad, panel+PG {k_ref:.3f} /rad "
          f"({100 * (k_su2 - k_ref) / k_ref:+.1f}%), thin aerofoil+PG {2 * np.pi / BETA:.3f} /rad")
    print(f"Zero-lift angle: SU2 {np.degrees(-c_su2 / k_su2):.2f} deg, panel {np.degrees(-c_ref / k_ref):.2f} deg, "
          f"thin aerofoil {np.degrees(a0_thin):.2f} deg")

    aa = np.linspace(lc["AoA"].min(), lc["AoA"].max(), 50)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(aa, 2 * np.pi * (np.radians(aa) - a0_thin) / BETA, "k:", label="Thin aerofoil theory + Prandtl–Glauert")
    ax.plot(lc["AoA"], lc["CL_panel_PG"], "k--", label="Vortex panel method + Prandtl–Glauert")
    ax.plot(lc["AoA"], lc["CL_SU2"], "o-", color=SU2_BLUE, markersize=8, label="SU2 Euler, medium mesh")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("Angle of attack (deg)")
    ax.set_ylabel("Lift coefficient $C_L$")
    ax.set_title(f"NACA 4412 inviscid lift curve, M = {MACH}")
    ax.legend()
    style(ax)
    fig.tight_layout()
    fig.savefig("liftcurve_naca4412.png", dpi=150)

    # ---------- 4. Mesh-convergence figure ----------
    h_rel = np.sqrt(LEVELS["fine"] / mc["cells"])               # grid spacing relative to fine
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    ax.plot(h_rel, mc["CL"], "o-", color=SU2_BLUE, markersize=8, label="SU2 Euler")
    ax.plot(0, fe, "s", mfc="white", mec=SU2_BLUE, markersize=9, label=f"Richardson extrapolation (p = {p:.1f})")
    ax.axhline(CL_ref, color="k", ls="--", lw=1.2, label="Panel method + PG")
    ax.set_xlabel("Relative grid spacing h / h$_{fine}$")
    ax.set_ylabel("$C_L$ at α = 5°")
    ax.set_title("Lift")
    ax.legend(fontsize=9)
    style(ax)
    ax = axes[1]
    ax.plot(h_rel, 1e4 * mc["CD"], "o-", color=SU2_BLUE, markersize=8, label="SU2 Euler")
    ax.axhline(0, color="k", ls="--", lw=1.2, label="Exact inviscid value (0)")
    ax.set_xlim(left=-0.2)
    ax.set_xlabel("Relative grid spacing h / h$_{fine}$")
    ax.set_ylabel("Spurious drag (counts, 1 count = 0.0001)")
    ax.set_title("Numerical drag")
    ax.legend(fontsize=9)
    style(ax)
    fig.suptitle(f"NACA 4412 mesh convergence (6k, 25k, 98k cells), M = {MACH}")
    fig.tight_layout()
    fig.savefig("mesh_convergence.png", dpi=150)
    print("\nSaved cp_validation.png, liftcurve_naca4412.png, mesh_convergence.png")
