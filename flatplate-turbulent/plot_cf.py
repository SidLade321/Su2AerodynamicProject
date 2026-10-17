import meshio
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RE_PER_M = 5e6          # REYNOLDS_NUMBER in the config, per metre (REYNOLDS_LENGTH = 1)
CHECK_X = 0.97          # standard station used by NASA's turbulence modelling resource

# --- Load SU2's wall data ---
m = meshio.read("surface_flow.vtk")
x = m.points[:, 0]
cf = m.point_data["Skin_Friction_Coefficient"]
cf = cf[:, 0] if cf.ndim > 1 else cf          # x-component = streamwise friction
yplus = m.point_data["Y_Plus"]

order = np.argsort(x)
x, cf, yplus = x[order], cf[order], yplus[order]

# --- Textbook turbulent flat-plate correlations (local Cf) ---
def prandtl(rex):  return 0.0592 * rex ** -0.2
def schultz(rex):  return 0.370 * np.log10(rex) ** -2.584
def white(rex):    return 0.455 / np.log(0.06 * rex) ** 2

xs = np.linspace(0.05, x.max(), 300)
rex = RE_PER_M * xs

# --- Numbers for the README ---
cf_check = np.interp(CHECK_X, x, cf)
print(f"SU2 Cf at x = {CHECK_X} m: {cf_check:.5f}")
for name, f in [("Prandtl 1/5 power", prandtl), ("Schultz-Grunow", schultz), ("White", white)]:
    ref = f(RE_PER_M * CHECK_X)
    print(f"  {name:18s}: {ref:.5f}  (SU2 differs by {100 * (cf_check - ref) / ref:+.1f}%)")

cf_avg = np.trapezoid(cf, x) / (x.max() - x.min())
print(f"\nAverage Cf from wall data: {cf_avg:.6f} (solver's integrated CD was 0.002859)")
print(f"Wall y+: max {yplus.max():.2f}, mean {yplus.mean():.2f}")

# --- Plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                               gridspec_kw={"height_ratios": [3, 1]})
ax1.plot(x, cf, color="#1f77b4", linewidth=2.5, label="SU2 (Spalart-Allmaras RANS)")
ax1.plot(xs, prandtl(rex), "k--", label="Prandtl 1/5 power law")
ax1.plot(xs, schultz(rex), "k:", label="Schultz-Grunow")
ax1.plot(xs, white(rex), "k-.", label="White")
ax1.axvline(CHECK_X, color="grey", linewidth=0.8)
ax1.set_ylim(0, 0.006)
ax1.set_ylabel("Skin friction coefficient $C_f$")
ax1.set_title("Turbulent flat plate, Re = 5×10⁶ per m, M = 0.2: SU2 vs correlations")
ax1.legend()
ax1.grid(alpha=0.3)

ax2.plot(x, yplus, color="#2ca02c")
ax2.set_xlabel("Distance along plate x (m)")
ax2.set_ylabel("Wall $y^+$")
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("cf_validation.png", dpi=150)
print("Saved cf_validation.png")
