"""
tools/final_survey_halo_enhanced.py
====================================
CODE-GEO V4.2 — Enhanced Survey Runner (High-Contrast Rendering)
=================================================================

Runs the calibrated V4.2 pipeline across all survey targets with
an enhanced rendering style emphasising the spatial structure of the
MOND interpolation field F'(Q).

Changelog vs original
---------------------
Original:
  - Applied gaussian_filter(rho_eff, sigma=15.0) as a post-processing step
    to produce "ghost halos" around galaxies. This was a cosmetic operation —
    not physics. It visually mimicked a halo by blurring the output density,
    but the blur was unrelated to any prediction of the theory.
  - Called old engine interface with arbitrary Q normalisation.
  - No physical calibration of input flux.

V4.2:
  - Gaussian post-processing removed entirely. The F'(Q) field already encodes
    the physically correct spatial distribution of gravitational enhancement.
    Where F'(Q) is small (deep MOND), Σ_eff/Σ_b is large — that IS the halo
    signal, derived from the theory rather than imposed by a blur kernel.
  - Uses EuclidLoader.load_and_calibrate() for physical Σ_b [kg/m²].
  - Four-panel rendering: Σ_b | F'(Q) | Σ_eff | radial profile.
  - Radial profile panel shows the DM enhancement vs radius — the key
    observable signature of the Mimetic-Conformal framework.

Note on "halo" signal
---------------------
In the V4.2 framework, the gravitational enhancement is strongest where
g_N is smallest (outer regions, low surface density). This is the MOND
regime (Q ≪ 1, F'(Q) ≈ √Q ≪ 1). The enhancement is therefore naturally
extended beyond the baryonic core — an emergent spatial halo from the
physics, not a post-processing blur.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.mimetic_engine import MimeticEngine
from tools.euclid_loader import EuclidLoader, TARGET_REGISTRY

# Survey manifest — same as run_full_survey.py
SURVEY_MANIFEST = {
    "bullet_cluster" : "mastDownload/HST/j90701010/j90701010_drz.fits",
    "abell_370"      : "mastDownload/HST/jabu01030/jabu01030_drz.fits",
    "el_gordo"       : "mastDownload/HST/jbqz31010/jbqz31010_drz.fits",
    "hudf"           : "mastDownload/HST/j8wc7c010/j8wc7c010_drz.fits",
}

FITS_ROOT  = "data/raw_fits/"
OUTPUT_DIR = "survey_outputs/"


def _log_display(arr: np.ndarray) -> np.ndarray:
    return np.log10(arr + 1e-40)


def radial_profile(arr: np.ndarray, centre: tuple = None) -> tuple:
    """
    Compute the azimuthally averaged radial profile of a 2D array.

    Returns
    -------
    r_bins   : 1-D array  pixel radius
    profile  : 1-D array  mean value at each radius
    """
    ny, nx = arr.shape
    cy, cx = centre or (ny // 2, nx // 2)
    y_idx, x_idx = np.indices((ny, nx))
    r = np.sqrt((x_idx - cx)**2 + (y_idx - cy)**2).astype(int)
    r_max = min(cy, cx, ny - cy, nx - cx)
    r_bins  = np.arange(0, r_max)
    profile = np.array([arr[r == ri].mean() if np.any(r == ri) else 0.0
                        for ri in r_bins])
    return r_bins, profile


def render_enhanced_report(
    target_key: str,
    sigma_b: np.ndarray,
    sigma_eff: np.ndarray,
    F_prime: np.ndarray,
    telemetry: dict,
    provenance: dict,
    dx: float,
    output_dir: str,
) -> str:
    """
    Four-panel enhanced report: Σ_b | F'(Q) | Σ_eff | radial DM profile.
    """
    try:
        from core.physics_constants import KPC_TO_M
    except ImportError:
        KPC_TO_M = 3.08567758e19

    fig = plt.figure(figsize=(22, 6), facecolor="#0a0a0a")
    gs  = gridspec.GridSpec(1, 4, figure=fig, wspace=0.35)

    # Panel 1: Baryonic input
    ax0 = fig.add_subplot(gs[0])
    im0 = ax0.imshow(_log_display(sigma_b), cmap="hot", origin="lower")
    ax0.set_title(
        f"Σ_b  [log₁₀ kg/m²]\nM/L = {provenance['ml_ratio_used']} M☉/L☉",
        color="white", fontsize=9,
    )
    plt.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)

    # Panel 2: MOND interpolation field
    ax1 = fig.add_subplot(gs[1])
    im1 = ax1.imshow(F_prime, cmap="viridis", origin="lower", vmin=0.0, vmax=1.0)
    ax1.set_title(
        f"F'(Q) — MOND field\n"
        f"mean={telemetry['F_prime_mean']:.3f}  "
        f"min={telemetry['F_prime_min']:.3f}",
        color="white", fontsize=9,
    )
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

    # Panel 3: Effective lensing density
    ax2 = fig.add_subplot(gs[2])
    im2 = ax2.imshow(_log_display(sigma_eff), cmap="coolwarm", origin="lower")
    ax2.set_title(
        f"Σ_eff  [log₁₀ kg/m²]\n"
        f"DM ratio: {telemetry['dm_ratio_mean']:.2f}x (mean)  "
        f"{telemetry['dm_ratio_peak']:.2f}x (peak)",
        color="white", fontsize=9,
    )
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

    # Panel 4: Radial DM enhancement profile
    ax3 = fig.add_subplot(gs[3])
    r_bins_b,   prof_b   = radial_profile(sigma_b)
    r_bins_eff, prof_eff = radial_profile(sigma_eff)

    # Convert pixel radius to kpc
    r_kpc = r_bins_b * dx / KPC_TO_M

    ax3.plot(r_kpc, prof_eff,  color="cyan",   lw=2,
             label="Σ_eff (mimetic)")
    ax3.plot(r_kpc, prof_b,    color="orange", lw=2, ls="--",
             label="Σ_b (baryons)")

    # Enhancement ratio profile
    ax3_r = ax3.twinx()
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_prof = np.where(prof_b > 0, prof_eff / prof_b, 1.0)
    ax3_r.plot(r_kpc, ratio_prof, color="lime", lw=1.5, ls=":",
               label="DM ratio (right)")
    ax3_r.set_ylabel("DM ratio", color="lime", fontsize=8)
    ax3_r.tick_params(axis="y", colors="lime", labelsize=7)
    ax3_r.set_ylim(bottom=0)

    ax3.set_yscale("log")
    ax3.set_xlabel("Radius [kpc]", color="white", fontsize=8)
    ax3.set_ylabel("Surface density [kg/m²]", color="white", fontsize=8)
    ax3.set_title("Radial Profile\nΣ_eff vs Σ_b", color="white", fontsize=9)
    ax3.tick_params(colors="white", labelsize=7)
    ax3.set_facecolor("#0a0a0a")
    ax3.spines[["bottom","top","left","right"]].set_color("#555555")
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_r.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2,
               fontsize=7, facecolor="#1a1a1a", labelcolor="white")

    for ax in [ax0, ax1, ax2]:
        ax.axis("off")

    fig.suptitle(
        f"CODE-GEO V4.2 — Enhanced Survey  |  {provenance['target_name']}  |  "
        f"z={provenance['redshift_z']}  |  "
        f"g_N_max/a₀={telemetry['g_N_max_ms2']/1.21e-10:.1f}",
        color="white", fontsize=10, fontweight="bold",
    )

    output_path = os.path.join(
        output_dir, f"SENTINEL_HALO_{target_key}_V42.png"
    )
    plt.savefig(output_path, facecolor="#0a0a0a", dpi=150)
    plt.close()
    return output_path


def run_halo_survey(targets: list = None):
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    engine = MimeticEngine()
    loader = EuclidLoader(data_path=FITS_ROOT)

    active = targets or list(SURVEY_MANIFEST.keys())

    print("🛰️  CODE-GEO V4.2 — ENHANCED HALO SURVEY")
    print(f"   Physical halo from F'(Q) — no cosmetic blur")
    print("=" * 55)

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
            out = render_enhanced_report(
                target_key, sigma_b, sigma_eff, F_prime,
                telemetry, provenance, dx, OUTPUT_DIR,
            )
            print(f"  [✓] {out}")

        except Exception as exc:
            print(f"  [ERROR] {target_key}: {exc}")

    print("\n🏁 ENHANCED SURVEY COMPLETE.")


if __name__ == "__main__":
    run_halo_survey()