import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from panel_method import naca_symmetric, lift_coefficient

MACH = 0.3
ANGLES = [0, 2, 4, 6, 8]
beta = np.sqrt(1 - MACH**2)                     # Prandtl-Glauert compressibility factor

rows = []
for a in ANGLES:
    h = pd.read_csv(f"history_aoa{a}.csv")
    h.columns = [c.strip().strip('"') for c in h.columns]
    rho_drop = h["rms[Rho]"].iloc[0] - h["rms[Rho]"].iloc[-1]
    rows.append({"AoA": a, "CL_SU2": h["CL"].iloc[-1], "CD_SU2": h["CD"].iloc[-1],
                 "Iterations": len(h), "Residual_drop": rho_drop})
res = pd.DataFrame(rows)

x, y = naca_symmetric(0.12, 100)
res["CL_panel_PG"] = [lift_coefficient(x, y, a) / beta for a in ANGLES]
res["CL_thin_PG"] = [2 * np.pi * np.radians(a) / beta for a in ANGLES]
ref = res["CL_panel_PG"].where(res["CL_panel_PG"].abs() > 1e-6)      # no % difference at 0 deg
res["Diff_%"] = 100 * (res["CL_SU2"] - ref) / ref

pd.set_option("display.float_format", "{:.4f}".format)
print(res.to_string(index=False))

slope_su2 = np.polyfit(np.radians(res["AoA"]), res["CL_SU2"], 1)[0]
slope_ref = np.polyfit(np.radians(res["AoA"]), res["CL_panel_PG"], 1)[0]
print(f"\nLift-curve slope: SU2 {slope_su2:.3f} /rad, panel+PG {slope_ref:.3f} /rad "
      f"({100 * (slope_su2 - slope_ref) / slope_ref:+.1f}%), thin aerofoil+PG {2 * np.pi / beta:.3f} /rad")

plt.figure(figsize=(8, 5.5))
plt.plot(res["AoA"], res["CL_thin_PG"], "k:", label="Thin aerofoil theory + Prandtl-Glauert")
plt.plot(res["AoA"], res["CL_panel_PG"], "k--", label="Vortex panel method + Prandtl-Glauert")
plt.plot(res["AoA"], res["CL_SU2"], "o-", color="#1f77b4", markersize=8, label="SU2 Euler")
plt.xlabel("Angle of attack (deg)")
plt.ylabel("Lift coefficient $C_L$")
plt.title(f"NACA 0012 inviscid lift curve, M = {MACH}")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("liftcurve_validation.png", dpi=150)
print("Saved liftcurve_validation.png")
