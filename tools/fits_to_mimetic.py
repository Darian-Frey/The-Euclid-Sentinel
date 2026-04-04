"""
tools/fits_to_mimetic.py
=========================
CODE-GEO V4.2 — Single-Target FITS Inspection Tool
====================================================

A focused diagnostic tool for running the V4.2 pipeline on a single
FITS file and producing a detailed inspection report. Useful for:
  - Verifying calibration of a newly downloaded FITS file
  - Debugging the pixel scale / WCS extraction for a new target
  - Quick single-target runs before committing to a full survey
  - Sanity-checking that Σ_b values are physically plausible

Unlike run_full_survey.py (which loops all targets), this script
takes a single FITS path and target key as arguments and produces
a detailed five-panel diagnostic output plus a console report.

Changelog vs original
---------------------
Original:
  - Hardcoded FITS_PATH to a single HUDF file.
  - Used scipy zoom to resize to 512×512 with no physical scale tracking.
  - Applied gaussian_filter(sigma=4.0) for "cluster-scale gas distribution"
    — an unjustified smoothing that altered the physics input.
  - Q_field peak forced to 2.0 with no physical derivation.
  - Reported DM ratio as if it were a measurement.

V4.2:
  - Accepts FITS path and target key as CLI arguments.
  - No arbitrary resizing or smoothing — data used at native resolution
    (center-cropped to remove edge artefacts, as per loader default).
  - Physical calibration via EuclidLoader: Σ_b [kg/m²] and dx [m].
  - Five diagnostic panels: Σ_b | F'(Q) | Σ_eff | g_N field | Q field.
  - Console report includes physical regime diagnosis and ICM warning.
  - M/L override available for quick sensitivity checks.

Usage
-----
  # Default target (bullet_cluster, registry M/L):
  python3 tools/fits_to_mimetic.py

  # Specify target and optional M/L override:
  python3 tools/fits_to_mimetic.py --target abell_370 --ml 3.5

  # Point at a non-registry FITS (requires --target for redshift/zeropoint):
  python3 tools/fits_to_mimetic.py --fits path/to/file.fits --target el_gordo
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.mimetic_engine import MimeticEngine
from tools.euclid_loader import EuclidLoader, TARGET_REGISTRY

try:
    from core.physics_constants import A0, KPC_TO_M
except ImportError:
    A0       = 1.21e-10
    KPC_TO_M = 3.08567758e19

# Default FITS paths per target (relative to FITS_ROOT)
_DEFAULT_FITS = {
    "bullet_cluster" : "mastDownload/HST/j90701010/j90701010_drz.fits",
    "abell_370"      : "mastDownload/HST/jabu01030/jabu01030_drz.fits",
    "el_gordo"       : "mastDownload/HST/jbqz31010/jbqz31010_drz.fits",
    "hudf"           : "mastDownload/HST/j8wc7c010/j8wc7c010_drz.fits",
}

FITS_ROOT  = "data/raw_fits/"
OUTPUT_DIR = "survey_outputs/"


def render_diagnostic(
    target_key: str,
    sigma_b: np.ndarray,
    sigma_eff: np.ndarray,
    F_prime: np.ndarray,
    g_N: np.ndarray,
    Q: np.ndarray,
    telemetry: dict,
    provenance: dict,
    dx: float,
    output_dir: str,
) -> str:
    """
    Five-panel diagnostic report.

    Panels
    ------
    1. Σ_b          — baryonic input (log scale)
    2. F'(Q)        — MOND interpolation field
    3. Σ_eff        — effective lensing density (log scale)
    4. g_N / a₀     — Newtonian acceleration in units of a₀ (log scale)
                       Shows where the map is in MOND vs Newtonian regime
    5. Q field      — kinetic invariant (log scale)
                       Q ≪ 1 → deep MOND  |  Q ≫ 1 → Newtonian
    """
    fig = plt.figure(figsize=(26, 5), facecolor="#0a0a0a")
    gs  = gridspec.GridSpec(1, 5, figure=fig, wspace=0.35)

    def _log(arr):
        return np.log10(np.maximum(arr, 1e-40))

    # Panel 1: Σ_b
    ax0 = fig.add_subplot(gs[0])
    im0 = ax0.imshow(_log(sigma_b), cmap="hot", origin="lower")
    ax0.set_title(
        f"Σ_b  [log₁₀ kg/m²]\nM/L={provenance['ml_ratio_used']}  z={provenance['redshift_z']}",
        color="white", fontsize=8,
    )
    plt.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)

    # Panel 2: F'(Q)
    ax1 = fig.add_subplot(gs[1])
    im1 = ax1.imshow(F_prime, cmap="viridis", origin="lower", vmin=0.0, vmax=1.0)
    ax1.set_title(
        f"F'(Q)  MOND field\n"
        f"mean={telemetry['F_prime_mean']:.3f}  min={telemetry['F_prime_min']:.4f}",
        color="white", fontsize=8,
    )
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

    # Panel 3: Σ_eff
    ax2 = fig.add_subplot(gs[2])
    im2 = ax2.imshow(_log(sigma_eff), cmap="coolwarm", origin="lower")
    ax2.set_title(
        f"Σ_eff  [log₁₀ kg/m²]\n"
        f"DM: {telemetry['dm_ratio_mean']:.2f}x (mean)  "
        f"{telemetry['dm_ratio_peak']:.2f}x (peak)",
        color="white", fontsize=8,
    )
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

    # Panel 4: g_N / a₀ (regime map)
    ax3 = fig.add_subplot(gs[3])
    g_over_a0 = g_N / A0
    im3 = ax3.imshow(
        np.log10(np.maximum(g_over_a0, 1e-6)),
        cmap="plasma", origin="lower",
    )
    ax3.set_title(
        f"log₁₀(g_N / a₀)\n"
        f"<0 = MOND  |  >0 = Newtonian  |  max={np.max(g_over_a0):.1f}",
        color="white", fontsize=8,
    )
    cb3 = plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)
    cb3.ax.axhline(y=0, color="white", lw=1.5, ls="--")   # marks g_N = a₀

    # Panel 5: Q field
    ax4 = fig.add_subplot(gs[4])
    im4 = ax4.imshow(
        np.log10(np.maximum(Q, 1e-10)),
        cmap="magma", origin="lower",
    )
    ax4.set_title(
        f"log₁₀(Q)  kinetic invariant\n"
        f"Q≪1=deep MOND  Q≫1=Newtonian  max={telemetry['Q_max']:.2e}",
        color="white", fontsize=8,
    )
    plt.colorbar(im4, ax=ax4, fraction=0.046, pad=0.04)

    for ax in [ax0, ax1, ax2, ax3, ax4]:
        ax.axis("off")

    fig.suptitle(
        f"CODE-GEO V4.2  Diagnostic  |  {provenance['target_name']}  |  "
        f"dx={provenance['dx_pc']:.2f} pc/pix  |  "
        f"g_N_max/a₀={telemetry['g_N_max_ms2']/A0:.1f}  |  "
        f"λ={telemetry['lam']}  β={telemetry['beta']}",
        color="white", fontsize=9, fontweight="bold",
    )

    output_path = os.path.join(
        output_dir, f"DIAGNOSTIC_{target_key}_V42.png"
    )
    plt.savefig(output_path, facecolor="#0a0a0a", dpi=150)
    plt.close()
    return output_path


def run_inspection(
    fits_path: str,
    target_key: str,
    ml_ratio: float | None = None,
):
    """
    Run the full diagnostic pipeline on a single FITS file.
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    engine = MimeticEngine()
    loader = EuclidLoader(data_path=FITS_ROOT)

    target_name = TARGET_REGISTRY[target_key]["name"]

    print("=" * 60)
    print(f"CODE-GEO V4.2 — FITS Inspection Tool")
    print(f"Target : {target_name}")
    print(f"FITS   : {fits_path}")
    print(f"M/L    : {ml_ratio if ml_ratio else 'registry default'}")
    print("=" * 60)

    # Check file exists
    full_path = os.path.join(FITS_ROOT, fits_path)
    if not os.path.exists(full_path):
        print(
            f"\n[ERROR] FITS file not found: {full_path}\n"
            f"Run tools/fetch_eso_data.py {target_key} to download."
        )
        return

    # Load and calibrate
    sigma_b, dx, provenance = loader.load_and_calibrate(
        fits_path  = fits_path,
        target_key = target_key,
        ml_ratio   = ml_ratio,
    )

    # Engine — extract Q and g_N for diagnostic panels
    # We call the internal methods directly to get intermediate fields
    Q, g_N = engine._compute_q_field(sigma_b, dx)
    sigma_eff, F_prime, telemetry = engine.compute_effective_density(sigma_b, dx)

    # Console diagnostic report
    print(f"\n── Physical Calibration ────────────────────────────────")
    print(f"  dx               : {dx:.4e} m  ({provenance['dx_pc']:.3f} pc/pixel)")
    print(f"  Σ_b peak         : {provenance['sigma_b_peak']:.4e} kg/m²")
    print(f"  Σ_b mean         : {provenance['sigma_b_mean']:.4e} kg/m²")

    print(f"\n── Q-Field & Acceleration ──────────────────────────────")
    print(f"  g_N max          : {telemetry['g_N_max_ms2']:.4e} m/s²")
    print(f"  g_N max / a₀     : {telemetry['g_N_max_ms2']/A0:.2f}")
    print(f"  Q max            : {telemetry['Q_max']:.4e}")
    print(f"  Q mean           : {telemetry['Q_mean']:.4e}")
    regime = "Newtonian core" if telemetry['g_N_max_ms2']/A0 > 10 else "MOND transition"
    print(f"  Regime at peak   : {regime}")

    print(f"\n── MOND Interpolation ──────────────────────────────────")
    print(f"  F'(Q) mean       : {telemetry['F_prime_mean']:.4f}")
    print(f"  F'(Q) min        : {telemetry['F_prime_min']:.4f}  (strongest MOND enhancement)")
    print(f"  DM ratio mean    : {telemetry['dm_ratio_mean']:.2f}x  ← PREDICTED by F'(Q)")
    print(f"  DM ratio peak    : {telemetry['dm_ratio_peak']:.2f}x")

    print(f"\n── Astrophysical Assumptions ───────────────────────────")
    print(f"  M/L used         : {provenance['ml_ratio_used']} M☉/L☉  ({provenance['ml_source']})")
    print(f"  Filter zeropoint : {provenance['zeropoint_AB']} AB")
    print(f"  ICM gas included : {provenance['gas_fraction'] is not None}")
    if provenance["gas_fraction"] is None:
        print(f"  ⚠  ICM WARNING   : {provenance['notes'][:80]}...")

    EuclidLoader.print_provenance(provenance)

    # Render diagnostic image
    out = render_diagnostic(
        target_key, sigma_b, sigma_eff, F_prime, g_N, Q,
        telemetry, provenance, dx, OUTPUT_DIR,
    )
    print(f"[✓] Diagnostic image → {out}")

    # M/L sensitivity hint
    print(
        f"\n  Tip: run with --ml 1.5 through --ml 5.0 to sweep M/L systematic.\n"
        f"  DM ratio will scale approximately linearly with M/L in this regime."
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CODE-GEO V4.2 — Single-target FITS inspection tool"
    )
    parser.add_argument(
        "--target",
        type=str,
        default="bullet_cluster",
        choices=list(TARGET_REGISTRY.keys()),
        help="Target key from TARGET_REGISTRY (default: bullet_cluster)",
    )
    parser.add_argument(
        "--fits",
        type=str,
        default=None,
        help=(
            "FITS path relative to data/raw_fits/. "
            "If not provided, uses the default path for --target."
        ),
    )
    parser.add_argument(
        "--ml",
        type=float,
        default=None,
        help="M/L ratio override in M☉/L☉. Default: registry value for target.",
    )
    args = parser.parse_args()

    fits_path = args.fits or _DEFAULT_FITS.get(args.target)
    if fits_path is None:
        print(f"[ERROR] No FITS path for target '{args.target}'. Use --fits.")
        sys.exit(1)

    run_inspection(
        fits_path  = fits_path,
        target_key = args.target,
        ml_ratio   = args.ml,
    )
