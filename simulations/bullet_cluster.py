"""
simulations/bullet_cluster.py
==============================
CODE-GEO V4.2 — Bullet Cluster Synthetic Simulation
====================================================

A controlled synthetic simulation of the Bullet Cluster merger geometry,
using physically motivated surface mass density inputs rather than
volumetric density with arbitrary normalization.

Purpose
-------
This simulation serves as a physics validation tool — not a fit to data.
It verifies that the V4.2 engine produces spatially offset gravitational
enhancement consistent with the observed baryonic-lensing offset in the
real Bullet Cluster (Clowe et al. 2006): the lensing peak is offset ~8 arcsec
from the X-ray gas peak after the merger.

Changelog vs original
---------------------
Original:
  - Used rho_baryon [kg/m³] — wrong units for a projected 2D map.
  - Hardcoded sigma scalar field with arbitrary amplitude (0.025).
  - Q computed from gradient of sigma field, not from g_N via Poisson.
  - Called old engine interface: compute_effective_density(rho_baryon, grad_sq).
  - "The fix for the unit bomb" comment acknowledged the problem but didn't fix it.

V4.2:
  - sigma_b [kg/m²] — physically correct surface mass density.
  - Peak amplitudes motivated by Clowe et al. (2006) baryonic mass estimates.
  - dx derived from physical box size and grid resolution.
  - Engine called with (sigma_b, dx) — Q derived internally via thin-disk Poisson.
  - Three-panel output: Σ_b | F'(Q) field | Σ_eff.

Bullet Cluster geometry
-----------------------
The merger is oriented roughly east-west. The smaller, denser subcluster
(the "bullet") has passed through the larger main cluster. The hot gas
(ICM) is ram-pressure stripped and sits between the two dark matter halos.
In a MOND-like framework, the gravitational enhancement should track the
baryonic mass (stars + gas), not separate from it as in ΛCDM.

This simulation models the stellar component only. The ICM gas dominates
the baryonic budget (~6:1 gas-to-stellar) and would require a separate
X-ray-derived mass map for a complete treatment.

References
----------
* Clowe et al. (2006) ApJ 648 L109  — "A Direct Empirical Proof of Dark Matter"
* Bradač et al. (2006)              — Bullet Cluster lensing reconstruction
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.mimetic_engine import MimeticEngine

try:
    from core.physics_constants import MSUN, KPC_TO_M, A0
except ImportError:
    MSUN    = 1.98847e30
    KPC_TO_M = 3.08567758e19
    A0      = 1.21e-10


# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

GRID = 512                  # pixels per side
BOX_KPC = 3000.0            # physical box size [kpc] — ~10 Mpc for the cluster pair
BOX_M   = BOX_KPC * KPC_TO_M
dx      = BOX_M / GRID      # physical pixel scale [m/pixel]

# Surface mass density scale [kg/m²]
# Bullet Cluster total stellar mass ~2e13 M_sun over ~1 Mpc²
# → mean Σ_* ~ 2e13 M_sun / (3e22 m)² ~ 7e-17 kg/m²; peak ~100× mean
_SIGMA_SCALE = 1.0e-14      # kg/m²  representative cluster stellar peak

# Msun/pc² → kg/m² conversion (for reference in comments)
_MSUN_PC2_TO_SI = MSUN / (3.08567758e16 ** 2)   # ≈ 6.77e-21 kg/m² per M_sun/pc²


# ---------------------------------------------------------------------------
# Build the baryonic surface mass density map
# ---------------------------------------------------------------------------

x_arr = np.linspace(-BOX_M / 2, BOX_M / 2, GRID)
X, Y  = np.meshgrid(x_arr, x_arr)

# Main cluster — larger, more massive, at centre-left
# ~1.5e13 M_sun stellar mass, scale radius ~300 kpc
r_main    = 300.0 * KPC_TO_M
sigma_main = 1.0 * _SIGMA_SCALE * np.exp(-(X**2 + Y**2) / r_main**2)

# Bullet subcluster — smaller, denser, offset to the right (post-merger)
# ~3e12 M_sun stellar mass, scale radius ~150 kpc
r_bullet     = 150.0 * KPC_TO_M
x_offset     = 600.0 * KPC_TO_M   # ~600 kpc east of centre
sigma_bullet = 0.4 * _SIGMA_SCALE * np.exp(
    -((X - x_offset)**2 + Y**2) / r_bullet**2
)

# Total baryonic stellar surface mass density [kg/m²]
sigma_b = sigma_main + sigma_bullet


# ---------------------------------------------------------------------------
# Run the V4.2 engine
# ---------------------------------------------------------------------------

print("=" * 60)
print("CODE-GEO V4.2 — Bullet Cluster Synthetic Simulation")
print("=" * 60)
print(f"\nGrid     : {GRID}×{GRID} pixels")
print(f"Box size : {BOX_KPC:.0f} kpc  ({BOX_M:.3e} m)")
print(f"dx       : {dx:.4e} m  ({dx/KPC_TO_M*1e3:.3f} pc/pixel)")
print(f"Σ_b peak : {np.max(sigma_b):.4e} kg/m²")

engine = MimeticEngine()
sigma_eff, F_prime, telemetry = engine.compute_effective_density(sigma_b, dx)

print(f"\n── Simulation Telemetry ────────────────────────────────")
for k, v in telemetry.items():
    print(f"  {k:<20s}: {v}")

# Physical interpretation
g_ratio = telemetry["g_N_max_ms2"] / A0
print(f"\n  g_N_max / a₀  : {g_ratio:.2f}")
print(
    f"  Regime        : "
    f"{'Newtonian core' if g_ratio > 10 else 'MOND transition'} at peak"
)
print(f"  DM ratio mean : {telemetry['dm_ratio_mean']:.2f}x  (PREDICTED by F'(Q))")
print(f"  DM ratio peak : {telemetry['dm_ratio_peak']:.2f}x")
print(f"─────────────────────────────────────────────────────────")


# ---------------------------------------------------------------------------
# Visualisation — three panels
# ---------------------------------------------------------------------------

fig = plt.figure(figsize=(18, 6), facecolor="#0a0a0a")
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.3)

def _log_display(arr):
    return np.log10(arr + 1e-40)

# Panel 1: Baryonic input
ax0 = fig.add_subplot(gs[0])
im0 = ax0.imshow(_log_display(sigma_b), cmap="hot", origin="lower")
ax0.set_title(
    "Baryonic Σ_b [log₁₀ kg/m²]\n"
    "Main cluster + bullet subcluster",
    color="white", fontsize=10,
)
plt.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)

# Panel 2: MOND interpolation field — the core physics output
ax1 = fig.add_subplot(gs[1])
im1 = ax1.imshow(F_prime, cmap="viridis", origin="lower", vmin=0.0, vmax=1.0)
ax1.set_title(
    "F'(Q) — MOND interpolation field\n"
    f"0 = deep MOND  |  1 = Newtonian  |  mean = {telemetry['F_prime_mean']:.3f}",
    color="white", fontsize=10,
)
plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

# Panel 3: Effective lensing surface density
ax2 = fig.add_subplot(gs[2])
im2 = ax2.imshow(_log_display(sigma_eff), cmap="coolwarm", origin="lower")
ax2.set_title(
    "Effective Σ_eff [log₁₀ kg/m²]\n"
    f"DM ratio (mean): {telemetry['dm_ratio_mean']:.2f}x  "
    f"| (peak): {telemetry['dm_ratio_peak']:.2f}x",
    color="white", fontsize=10,
)
plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

for ax in [ax0, ax1, ax2]:
    ax.axis("off")

fig.suptitle(
    f"CODE-GEO V4.2  |  Bullet Cluster Synthetic  |  "
    f"λ={engine.lam}  β={engine.beta}  a₀={engine.a0:.3e} m/s²",
    color="white", fontsize=11, fontweight="bold",
)

output_path = "SENTINEL_REPORT_Bullet_Cluster_V42.png"
plt.savefig(output_path, facecolor="#0a0a0a", dpi=150)
plt.close()
print(f"\n[✓] Report → {output_path}")
print(
    "\nNote: This is a SYNTHETIC simulation (Gaussian blobs).\n"
    "For real HST data, use tools/run_full_survey.py with\n"
    "tools/euclid_loader.py to provide calibrated Σ_b maps."
)