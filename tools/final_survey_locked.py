"""
tools/final_survey_locked.py
=============================
CODE-GEO V4.2 — Locked Survey Runner (Standard Three-Panel Output)
===================================================================

Runs the calibrated V4.2 pipeline across all survey targets with a clean
three-panel output matching the standard SENTINEL_REPORT format used by
run_full_survey.py. This script is the "production locked" variant —
no experimental rendering, no radial profiles, just the canonical output.

Use this when you want output images that are directly comparable
across pipeline versions or parameter sweeps.

Changelog vs original
---------------------
Original:
  - Loaded raw FITS with no physical calibration.
  - Computed Q via arbitrary gradient normalisation (peak = 2.0 / 10.0).
  - Called old engine interface: compute_effective_density(rho_baryon_phys, Q_field).
  - Results were artefacts of the hardcoded 12.0 factor in the engine.

V4.2:
  - EuclidLoader.load_and_calibrate() provides Σ_b [kg/m²] and dx [m].
  - Engine derives Q from thin-disk Poisson — no free normalisation.
  - DM ratio is a physics prediction from F'(Q).
  - Output labelled clearly with physical parameters (M/L, z, g_N/a₀).
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.mimetic_engine import MimeticEngine
from tools.euclid_loader import EuclidLoader, TARGET_REGISTRY

SURVEY_MANIFEST = {
    "bullet_cluster" : "mastDownload/HST/j90701010/j90701010_drz.fits",
    "abell_370"      : "mastDownload/HST/jabu01030/jabu01030_drz.fits",
    "el_gordo"       : "mastDownload/HST/jbqz31010/jbqz31010_drz.fits",
    "hudf"           : "mastDownload/HST/j8wc7c010/j8wc7c010_drz.fits",
}

FITS_ROOT  = "data/raw_fits/"
OUTPUT_DIR = "survey_outputs/"


def render_locked_report(
    target_key: str,
    sigma_b: np.ndarray,
    sigma_eff: np.ndarray,
    F_prime: np.ndarray,
    telemetry: dict,
    provenance: dict,
    output_dir: str,
) -> str:
    """
    Standard three-panel report: Σ_b | F'(Q) | Σ_eff.
    Output filename matches SENTINEL_REPORT_{target}.png convention.
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), facecolor="#0a0a0a")

    def _log(arr):
        return np.log10(arr + 1e-40)

    im0 = axes[0].imshow(_log(sigma_b), cmap="hot", origin="lower")
    axes[0].set_title(
        f"Baryonic Σ_b  [log₁₀ kg/m²]\n"
        f"M/L = {provenance['ml_ratio_used']} M☉/L☉  |  z = {provenance['redshift_z']}",
        color="white", fontsize=9,
    )
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(F_prime, cmap="viridis", origin="lower", vmin=0.0, vmax=1.0)
    axes[1].set_title(
        f"F'(Q) — MOND interpolation field\n"
        f"mean = {telemetry['F_prime_mean']:.3f}  "
        f"min = {telemetry['F_prime_min']:.3f}",
        color="white", fontsize=9,
    )
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    im2 = axes[2].imshow(_log(sigma_eff), cmap="coolwarm", origin="lower")
    axes[2].set_title(
        f"Effective Σ_eff  [log₁₀ kg/m²]\n"
        f"DM ratio  mean: {telemetry['dm_ratio_mean']:.2f}x  "
        f"peak: {telemetry['dm_ratio_peak']:.2f}x",
        color="white", fontsize=9,
    )
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    for ax in axes:
        ax.axis("off")

    fig.patch.set_facecolor("#0a0a0a")
    fig.suptitle(
        f"CODE-GEO V4.2  |  {provenance['target_name']}  |  "
        f"g_N_max/a₀ = {telemetry['g_N_max_ms2'] / 1.21e-10:.1f}  |  "
        f"λ={telemetry['lam']}  β={telemetry['beta']}",
        color="white", fontsize=10, fontweight="bold",
    )
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"SENTINEL_REPORT_{target_key}.png")
    plt.savefig(output_path, facecolor="#0a0a0a", dpi=150)
    plt.close()
    return output_path


def run_locked_survey(targets: list = None):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    engine = MimeticEngine()
    loader = EuclidLoader(data_path=FITS_ROOT)

    active = targets or list(SURVEY_MANIFEST.keys())

    print("🛰️  CODE-GEO V4.2 — LOCKED SURVEY (Standard Output)")
    print("=" * 55)

    results = []
    for target_key in active:
        fits_path = SURVEY_MANIFEST.get(target_key)
        if fits_path is None:
            print(f"[SKIP] '{target_key}' not in manifest.")
            continue

        full_path = os.path.join(FITS_ROOT, fits_path)
        if not os.path.exists(full_path):
            print(f"[SKIP] {full_path} not found.")
            continue

        print(f"\n🔭  {TARGET_REGISTRY[target_key]['name']}")

        try:
            sigma_b, dx, provenance = loader.load_and_calibrate(
                fits_path=fits_path, target_key=target_key
            )
            sigma_eff, F_prime, telemetry = engine.compute_effective_density(
                sigma_b, dx
            )
            out = render_locked_report(
                target_key, sigma_b, sigma_eff, F_prime,
                telemetry, provenance, OUTPUT_DIR,
            )
            print(f"  [✓] {out}")
            print(
                f"  DM ratio: mean={telemetry['dm_ratio_mean']:.2f}x  "
                f"peak={telemetry['dm_ratio_peak']:.2f}x  "
                f"g_N/a₀={telemetry['g_N_max_ms2']/1.21e-10:.1f}"
            )
            results.append((target_key, telemetry, provenance))

        except Exception as exc:
            print(f"  [ERROR] {target_key}: {exc}")

    print(f"\n🏁 LOCKED SURVEY COMPLETE — {len(results)}/{len(active)} targets.")

    if results:
        print(f"\n  {'Target':<25} {'DM mean':>9} {'DM peak':>9} {'g_N/a₀':>8}")
        print(f"  {'-'*53}")
        for tk, tel, _ in results:
            print(
                f"  {tk:<25} "
                f"{tel['dm_ratio_mean']:>9.2f}x "
                f"{tel['dm_ratio_peak']:>9.2f}x "
                f"{tel['g_N_max_ms2']/1.21e-10:>8.1f}"
            )


if __name__ == "__main__":
    run_locked_survey()
