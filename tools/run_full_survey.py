"""
tools/run_full_survey.py
========================
CODE-GEO V4.2 — Calibrated Full Survey Runner
==============================================

Executes the Mimetic-Conformal V4.2 pipeline across all four primary survey
targets, producing physically calibrated output maps and a real SUMMARY_DATA.csv
where every number is derived from the physics, not imposed.

Changelog vs original
---------------------
Original:
  - Loaded raw FITS flux directly with no physical calibration.
  - Q-field normalised to an arbitrary peak of 2.0.
  - Density scaled by hardcoded 8e-23 with no physical derivation.
  - Called old two-argument engine interface: compute_effective_density(rho, Q).
  - SUMMARY_DATA.csv contained results artefacted by the 12.0 hardcoded factor.
  - Silent failures swallowed real error information.

V4.2:
  - Uses EuclidLoader.load_and_calibrate() → (sigma_b, dx, provenance).
  - Uses MimeticEngine.compute_effective_density(sigma_b, dx) — V4.2 interface.
  - DM ratio is a physical prediction, not an imposed constant.
  - SUMMARY_DATA.csv records real telemetry: Q ranges, g_N, F'(Q), DM ratios.
  - Provenance dict logged per target — every assumption auditable.
  - Per-target error handling distinguishes missing files from physics failures.
  - Sensitivity sweep (M/L variation) available via --ml-sweep CLI flag.

Usage
-----
  # Standard survey (registry M/L defaults):
  python3 tools/run_full_survey.py

  # M/L sensitivity sweep for systematic error budget:
  python3 tools/run_full_survey.py --ml-sweep

  # Single target:
  python3 tools/run_full_survey.py --target bullet_cluster
"""

import os
import sys
import csv
import argparse
import warnings
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.mimetic_engine import MimeticEngine
from tools.euclid_loader import EuclidLoader, TARGET_REGISTRY

# ---------------------------------------------------------------------------
# Survey configuration
# ---------------------------------------------------------------------------

# Maps target_key → FITS path relative to the loader's data_path.
# Update paths here if your local HST download directory differs.
SURVEY_MANIFEST = {
    "bullet_cluster" : "mastDownload/HST/j90701010/j90701010_drz.fits",
    "abell_370"      : "mastDownload/HST/jabu01030/jabu01030_drz.fits",
    "el_gordo"       : "mastDownload/HST/jbqz31010/jbqz31010_drz.fits",
    "hudf"           : "mastDownload/HST/j8wc7c010/j8wc7c010_drz.fits",
}

# Root directory for FITS data (loader prepends this to SURVEY_MANIFEST paths)
FITS_ROOT = "data/raw_fits/"

# Output directory for report images and CSV
OUTPUT_DIR = "survey_outputs/"

# M/L values for sensitivity sweep (--ml-sweep flag)
ML_SWEEP_VALUES = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]

# CSV output column order
CSV_FIELDS = [
    "timestamp",
    "target_key",
    "target_name",
    "redshift_z",
    "ml_ratio_used",
    "ml_source",
    "dx_m",
    "dx_pc",
    "Q_mean",
    "Q_max",
    "g_N_max_ms2",
    "g_N_max_over_a0",
    "F_prime_mean",
    "F_prime_min",
    "dm_ratio_mean",
    "dm_ratio_peak",
    "sigma_b_peak_kg_m2",
    "sigma_eff_peak_kg_m2",
    "lam",
    "beta",
    "a0",
    "gas_icm_included",
    "notes",
]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def render_report(
    target_key: str,
    sigma_b: np.ndarray,
    sigma_eff: np.ndarray,
    F_prime: np.ndarray,
    telemetry: dict,
    provenance: dict,
    output_dir: str,
) -> str:
    """
    Render a three-panel report image for one survey target.

    Panels
    ------
    Left  : log₁₀ Σ_b   — baryonic surface mass density input
    Centre: F'(Q) field  — MOND interpolation field (0=deep MOND, 1=Newtonian)
    Right : log₁₀ Σ_eff  — effective lensing surface density (physics output)

    Returns
    -------
    output_path : str  path to the saved PNG
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), facecolor="#0a0a0a")

    # Panel 1: Baryonic input
    im0 = axes[0].imshow(
        np.log10(sigma_b + 1e-40), cmap="hot", origin="lower"
    )
    axes[0].set_title(
        f"Baryonic Σ_b [log₁₀ kg/m²]\n"
        f"M/L = {provenance['ml_ratio_used']} M☉/L☉  |  "
        f"z = {provenance['redshift_z']}",
        color="white", fontsize=9,
    )
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    # Panel 2: MOND interpolation field
    im1 = axes[1].imshow(
        F_prime, cmap="viridis", origin="lower", vmin=0.0, vmax=1.0
    )
    axes[1].set_title(
        f"F'(Q) — MOND interpolation field\n"
        f"0 = deep MOND  |  1 = Newtonian  |  "
        f"mean = {telemetry['F_prime_mean']:.3f}",
        color="white", fontsize=9,
    )
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    # Panel 3: Effective lensing density
    im2 = axes[2].imshow(
        np.log10(sigma_eff + 1e-40), cmap="coolwarm", origin="lower"
    )
    axes[2].set_title(
        f"Effective Σ_eff [log₁₀ kg/m²]\n"
        f"DM ratio (mean): {telemetry['dm_ratio_mean']:.2f}x  |  "
        f"(peak): {telemetry['dm_ratio_peak']:.2f}x",
        color="white", fontsize=9,
    )
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    for ax in axes:
        ax.axis("off")
        ax.title.set_color("white")

    plt.suptitle(
        f"CODE-GEO V4.2  |  {provenance['target_name']}  |  "
        f"g_N_max = {telemetry['g_N_max_ms2']:.2e} m/s²  |  "
        f"Q_max = {telemetry['Q_max']:.2e}",
        color="white", fontsize=10, fontweight="bold",
    )
    fig.patch.set_facecolor("#0a0a0a")
    plt.tight_layout()

    output_path = os.path.join(output_dir, f"SENTINEL_REPORT_{target_key}.png")
    plt.savefig(output_path, facecolor="#0a0a0a", dpi=150)
    plt.close()
    return output_path


# ---------------------------------------------------------------------------
# Single target pipeline
# ---------------------------------------------------------------------------

def run_target(
    target_key: str,
    fits_path: str,
    engine: MimeticEngine,
    loader: EuclidLoader,
    output_dir: str,
    ml_ratio: float | None = None,
) -> dict | None:
    """
    Run the full V4.2 pipeline for one survey target.

    Parameters
    ----------
    target_key : str   registry key (e.g. 'bullet_cluster')
    fits_path  : str   path to FITS file relative to loader.data_path
    engine     : MimeticEngine instance
    loader     : EuclidLoader instance
    output_dir : str   directory for output images
    ml_ratio   : float | None   M/L override (None = registry default)

    Returns
    -------
    row : dict of CSV-ready result fields, or None on failure
    """
    target_name = TARGET_REGISTRY[target_key]["name"]
    ml_label    = f"M/L={ml_ratio}" if ml_ratio is not None else "M/L=registry"

    print(f"\n{'='*60}")
    print(f"  TARGET: {target_name}  [{ml_label}]")
    print(f"{'='*60}")

    # ── 1. Check FITS file exists ─────────────────────────────────────
    full_path = os.path.join(loader.data_path, fits_path)
    if not os.path.exists(full_path):
        print(
            f"  [SKIP] FITS file not found: {full_path}\n"
            f"  Run tools/fetch_eso_data.py to download."
        )
        return None

    # ── 2. Load and calibrate ─────────────────────────────────────────
    try:
        sigma_b, dx, provenance = loader.load_and_calibrate(
            fits_path  = fits_path,
            target_key = target_key,
            ml_ratio   = ml_ratio,
        )
    except Exception as exc:
        print(f"  [ERROR] Calibration failed for {target_key}: {exc}")
        return None

    # ── 3. Engine: effective density ──────────────────────────────────
    try:
        sigma_eff, F_prime, telemetry = engine.compute_effective_density(
            sigma_b, dx
        )
    except Exception as exc:
        print(f"  [ERROR] Engine failed for {target_key}: {exc}")
        return None

    # ── 4. Render report image ────────────────────────────────────────
    try:
        img_path = render_report(
            target_key, sigma_b, sigma_eff, F_prime,
            telemetry, provenance, output_dir,
        )
        print(f"  [✓] Report image → {img_path}")
    except Exception as exc:
        print(f"  [WARN] Render failed (non-fatal): {exc}")
        img_path = "render_failed"

    # ── 5. Print provenance ───────────────────────────────────────────
    EuclidLoader.print_provenance(provenance)

    # ── 6. Assemble CSV row ───────────────────────────────────────────
    row = {
        "timestamp"            : datetime.utcnow().isoformat(),
        "target_key"           : target_key,
        "target_name"          : target_name,
        "redshift_z"           : provenance["redshift_z"],
        "ml_ratio_used"        : provenance["ml_ratio_used"],
        "ml_source"            : provenance["ml_source"],
        "dx_m"                 : f"{provenance['dx_m']:.6e}",
        "dx_pc"                : f"{provenance['dx_pc']:.4f}",
        "Q_mean"               : f"{telemetry['Q_mean']:.6e}",
        "Q_max"                : f"{telemetry['Q_max']:.6e}",
        "g_N_max_ms2"          : f"{telemetry['g_N_max_ms2']:.6e}",
        "g_N_max_over_a0"      : f"{telemetry['g_N_max_ms2'] / engine.a0:.4f}",
        "F_prime_mean"         : f"{telemetry['F_prime_mean']:.6f}",
        "F_prime_min"          : f"{telemetry['F_prime_min']:.6f}",
        "dm_ratio_mean"        : f"{telemetry['dm_ratio_mean']:.4f}",
        "dm_ratio_peak"        : f"{telemetry['dm_ratio_peak']:.4f}",
        "sigma_b_peak_kg_m2"   : f"{provenance['sigma_b_peak']:.6e}",
        "sigma_eff_peak_kg_m2" : f"{np.max(sigma_eff):.6e}",
        "lam"                  : telemetry["lam"],
        "beta"                 : telemetry["beta"],
        "a0"                   : f"{telemetry['a0']:.6e}",
        "gas_icm_included"     : provenance["gas_fraction"] is not None,
        "notes"                : provenance["notes"],
    }

    print(
        f"\n  ── Result summary ──────────────────────────────────────\n"
        f"  g_N_max / a₀  : {float(row['g_N_max_over_a0']):.2f}  "
        f"({'Newtonian' if float(row['g_N_max_over_a0']) > 10 else 'MOND-transition'} regime at peak)\n"
        f"  Q_max         : {telemetry['Q_max']:.4e}\n"
        f"  F'(Q) mean    : {telemetry['F_prime_mean']:.4f}\n"
        f"  DM ratio mean : {telemetry['dm_ratio_mean']:.2f}x\n"
        f"  DM ratio peak : {telemetry['dm_ratio_peak']:.2f}x\n"
        f"  ───────────────────────────────────────────────────────"
    )

    return row


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

def write_csv(rows: list[dict], output_dir: str, suffix: str = "") -> str:
    """Write survey results to SUMMARY_DATA.csv."""
    csv_path = os.path.join(output_dir, f"SUMMARY_DATA{suffix}.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n[Survey] SUMMARY_DATA → {csv_path}")
    return csv_path


# ---------------------------------------------------------------------------
# Main survey runner
# ---------------------------------------------------------------------------

def run_survey(targets: list[str] | None = None, ml_sweep: bool = False):
    """
    Run the full survey pipeline.

    Parameters
    ----------
    targets  : list of target_key strings to process, or None for all
    ml_sweep : if True, run each target at all ML_SWEEP_VALUES M/L ratios
               and write a separate sensitivity CSV
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    engine = MimeticEngine()
    loader = EuclidLoader(data_path=FITS_ROOT)

    active_targets = targets or list(SURVEY_MANIFEST.keys())

    print(f"\n🛰️  CODE-GEO V4.2 — CALIBRATED SURVEY")
    print(f"    Targets    : {active_targets}")
    print(f"    Engine     : λ={engine.lam} | β={engine.beta} | a₀={engine.a0:.3e}")
    print(f"    Output dir : {OUTPUT_DIR}")
    if ml_sweep:
        print(f"    M/L sweep  : {ML_SWEEP_VALUES}")

    # ── Standard survey (registry defaults) ──────────────────────────
    standard_rows = []
    for target_key in active_targets:
        if target_key not in SURVEY_MANIFEST:
            print(f"[WARN] '{target_key}' not in SURVEY_MANIFEST — skipping.")
            continue
        row = run_target(
            target_key = target_key,
            fits_path  = SURVEY_MANIFEST[target_key],
            engine     = engine,
            loader     = loader,
            output_dir = OUTPUT_DIR,
            ml_ratio   = None,   # use registry default
        )
        if row is not None:
            standard_rows.append(row)

    if standard_rows:
        write_csv(standard_rows, OUTPUT_DIR)

    # ── M/L sensitivity sweep ─────────────────────────────────────────
    if ml_sweep:
        print(f"\n\n{'='*60}")
        print(f"  M/L SENSITIVITY SWEEP")
        print(f"{'='*60}")

        sweep_rows = []
        for target_key in active_targets:
            if target_key not in SURVEY_MANIFEST:
                continue
            for ml in ML_SWEEP_VALUES:
                row = run_target(
                    target_key = target_key,
                    fits_path  = SURVEY_MANIFEST[target_key],
                    engine     = engine,
                    loader     = loader,
                    output_dir = OUTPUT_DIR,
                    ml_ratio   = ml,
                )
                if row is not None:
                    sweep_rows.append(row)

        if sweep_rows:
            write_csv(sweep_rows, OUTPUT_DIR, suffix="_ML_SWEEP")
            _print_sweep_summary(sweep_rows)

    # ── Final summary ─────────────────────────────────────────────────
    print(f"\n\n🏁  SURVEY COMPLETE")
    print(f"    Processed : {len(standard_rows)} / {len(active_targets)} targets")
    print(f"    Output    : {OUTPUT_DIR}")

    if standard_rows:
        print(f"\n    ── DM Ratio Summary (registry M/L) ─────────────────")
        print(f"    {'Target':<25} {'DM (mean)':>10} {'DM (peak)':>10} {'g_N/a₀':>8}")
        print(f"    {'-'*55}")
        for row in standard_rows:
            print(
                f"    {row['target_key']:<25} "
                f"{float(row['dm_ratio_mean']):>10.2f}x "
                f"{float(row['dm_ratio_peak']):>10.2f}x "
                f"{float(row['g_N_max_over_a0']):>8.2f}"
            )
        print(f"    {'-'*55}")
        print(
            "    Note: DM ratio is a physics PREDICTION from F'(Q).\n"
            "    It is NOT imposed. Compare against weak-lensing observations\n"
            "    to evaluate the Mimetic-Conformal V4.2 framework."
        )


def _print_sweep_summary(sweep_rows: list[dict]):
    """Print a compact M/L sensitivity table grouped by target."""
    print(f"\n    ── M/L Sensitivity Summary ──────────────────────────────")
    print(f"    {'Target':<20} {'M/L':>6} {'DM (mean)':>10} {'DM (peak)':>10}")
    print(f"    {'-'*50}")
    current_target = None
    for row in sweep_rows:
        if row["target_key"] != current_target:
            if current_target is not None:
                print(f"    {'-'*50}")
            current_target = row["target_key"]
        print(
            f"    {row['target_key']:<20} "
            f"{float(row['ml_ratio_used']):>6.1f} "
            f"{float(row['dm_ratio_mean']):>10.2f}x "
            f"{float(row['dm_ratio_peak']):>10.2f}x"
        )
    print(f"    {'-'*50}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="CODE-GEO V4.2 Calibrated Survey Runner"
    )
    parser.add_argument(
        "--target",
        type=str,
        default=None,
        help=(
            "Run a single target by key "
            f"(choices: {list(SURVEY_MANIFEST.keys())}). "
            "Default: all targets."
        ),
    )
    parser.add_argument(
        "--ml-sweep",
        action="store_true",
        help=(
            f"Run M/L sensitivity sweep at values {ML_SWEEP_VALUES} "
            "for systematic error budget. Writes SUMMARY_DATA_ML_SWEEP.csv."
        ),
    )
    args = parser.parse_args()

    targets = [args.target] if args.target else None
    run_survey(targets=targets, ml_sweep=args.ml_sweep)