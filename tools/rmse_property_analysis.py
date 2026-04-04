"""
tools/rmse_property_analysis.py
================================
CODE-GEO V4.2 — RMSE vs Galaxy Property Correlation Analysis
=============================================================

Answers the key referee question: are the 56 poor-fit galaxies (RMSE ≥ 20 km/s)
a random sample, or do they cluster in specific regions of galaxy property space?

Properties extracted from the SPARC rotmod files and correlated with RMSE:

  1. V_flat         [km/s]    — galaxy mass proxy (last measured V_obs)
  2. g_bar/a₀ max  [–]       — dynamical regime at peak acceleration
  3. Bulge fraction [–]       — f_bul = Υ_bul V_bul² / V_bar² at outer radius
  4. Gas fraction   [–]       — f_gas = V_gas² / V_bar² at outer radius
  5. log SB_disk    [L⊙/pc²]  — central disk surface brightness (innermost point)
  6. MOND fraction  [–]       — fraction of radii with g_bar/a₀ < 1
  7. Profile shape  [–]       — (V_last - V_peak) / V_peak (declining < 0, flat ≈ 0)
  8. n_points       [–]       — number of rotation curve data points

Spearman rank correlation ρ and p-value computed for each property.
High |ρ| with low p identifies properties that predict fit quality.

Outputs
-------
  survey_outputs/RMSE_PROPERTY_ANALYSIS.png  — 8-panel figure
  survey_outputs/RMSE_PROPERTY_ANALYSIS.csv  — full property table
  Console: ranked correlation table for paper Table 2

Usage
-----
  python3 tools/rmse_property_analysis.py

References
----------
  Lelli et al. (2016) AJ 152 157   — SPARC database
  Li et al. (2018) A&A 615 A3      — MOND fits to SPARC (comparison baseline)
"""

import os
import sys
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.physics_constants import A0, KPC_TO_M, KM_S_TO_M_S

SPARC_DIR    = "data/sparc/"
RESULTS_CSV  = "survey_outputs/SPARC_RESULTS_V42.csv"
OUTPUT_DIR   = "survey_outputs/"
ML_DISK      = 0.50
ML_BUL       = 0.70
POOR_FIT_KMS = 20.0


# ---------------------------------------------------------------------------
# Load SPARC_RESULTS_V42.csv
# ---------------------------------------------------------------------------

def load_results(path: str) -> dict[str, dict]:
    """Load per-galaxy RMSE results from sparc_refinery output."""
    results = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            results[row["galaxy"]] = {
                "rmse_kms"         : float(row["rmse_kms"]),
                "v_flat_kms"       : float(row["v_flat_kms"]),
                "g_bar_max_over_a0": float(row["g_bar_max_over_a0"]),
                "n_points"         : int(row["n_points"]),
            }
    print(f"[Analysis] Loaded {len(results)} galaxies from {path}")
    return results


# ---------------------------------------------------------------------------
# Extract properties from rotmod files
# ---------------------------------------------------------------------------

def extract_properties(dat_path: str) -> dict | None:
    """
    Parse a SPARC rotmod file and extract galaxy properties.

    SPARC rotmod columns:
        Rad[kpc]  Vobs  errV  Vgas  Vdisk  Vbul  SBdisk  SBbul
    All velocities in km/s; SB in L_sun/pc².
    """
    rows = []
    try:
        with open(dat_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    rows.append([float(p) for p in parts[:8]])
                except ValueError:
                    continue
    except OSError:
        return None

    if len(rows) < 3:
        return None

    data  = np.array(rows)
    r_kpc = data[:, 0]
    v_obs = data[:, 1]
    err_v = data[:, 2]
    v_gas = data[:, 3]
    v_disk= data[:, 4]
    v_bul = data[:, 5]
    sb_disk = data[:, 6] if data.shape[1] > 6 else np.zeros(len(r_kpc))
    sb_bul  = data[:, 7] if data.shape[1] > 7 else np.zeros(len(r_kpc))

    # Quality cut — positive errV
    valid = (r_kpc > 0) & (v_obs > 0) & (err_v > 0)
    if valid.sum() < 3:
        return None

    r_kpc = r_kpc[valid]; v_obs = v_obs[valid]; err_v = err_v[valid]
    v_gas = v_gas[valid]; v_disk= v_disk[valid]; v_bul = v_bul[valid]
    sb_disk = sb_disk[valid]; sb_bul = sb_bul[valid]

    # ── Baryonic velocity ─────────────────────────────────────────────────
    # V_bar² = V_gas² + Υ_d V_disk² + Υ_b V_bul²  (sign convention: SPARC)
    v_bar_sq = (
        np.sign(v_gas)  * v_gas**2
        + ML_DISK * np.sign(v_disk) * v_disk**2
        + ML_BUL  * np.sign(v_bul)  * v_bul**2
    )
    v_bar_sq = np.maximum(v_bar_sq, 0.0)
    v_bar    = np.sqrt(v_bar_sq)

    # ── Property 1: Bulge fraction at outer radius ────────────────────────
    # Use mean of outer 30% of radii to reduce noise
    n_outer  = max(1, len(r_kpc) // 3)
    v_bar_outer = v_bar[-n_outer:]
    v_bul_outer = v_bul[-n_outer:]
    bul_contribution = ML_BUL * np.mean(v_bul_outer**2)
    bar_total        = np.mean(v_bar_sq[-n_outer:])
    bulge_fraction   = bul_contribution / (bar_total + 1e-10)
    bulge_fraction   = float(np.clip(bulge_fraction, 0.0, 1.0))

    # ── Property 2: Gas fraction at outer radius ──────────────────────────
    gas_contribution = np.mean(v_gas[-n_outer:]**2)
    gas_fraction     = gas_contribution / (bar_total + 1e-10)
    gas_fraction     = float(np.clip(gas_fraction, 0.0, 1.0))

    # ── Property 3: Central disk surface brightness ───────────────────────
    # Use median of inner 3 points to reduce noise from beam smearing
    n_inner    = min(3, len(sb_disk))
    sb_central = float(np.median(sb_disk[:n_inner]))
    log_sb     = float(np.log10(sb_central + 1e-3))   # [log L_sun/pc²]

    # ── Property 4: MOND fraction ─────────────────────────────────────────
    # Fraction of radii where g_bar/a₀ < 1 (in MOND regime)
    r_m  = r_kpc * KPC_TO_M
    g_bar_si = v_bar_sq * (KM_S_TO_M_S**2) / (r_m + 1e-30)
    mond_frac = float(np.mean(g_bar_si < A0))

    # ── Property 5: Profile shape (declining?) ────────────────────────────
    # (V_last - V_peak) / V_peak: negative → declining, positive → rising
    v_peak  = float(np.max(v_obs))
    v_last  = float(v_obs[-1])
    profile_shape = (v_last - v_peak) / (v_peak + 1e-10)

    # ── Property 6: Inner rise rate ───────────────────────────────────────
    # Normalised slope of inner curve: steep rise → compact mass distribution
    if len(v_obs) >= 4:
        inner_slope = (v_obs[3] - v_obs[0]) / (v_obs[-1] + 1e-10)
    else:
        inner_slope = 0.0

    return {
        "bulge_fraction" : bulge_fraction,
        "gas_fraction"   : gas_fraction,
        "log_sb_central" : log_sb,
        "mond_fraction"  : mond_frac,
        "profile_shape"  : float(profile_shape),
        "inner_slope"    : float(inner_slope),
        "v_peak_kms"     : float(v_peak),
    }


# ---------------------------------------------------------------------------
# Build full property table
# ---------------------------------------------------------------------------

def build_property_table(results: dict, sparc_dir: str) -> list[dict]:
    """
    Cross-match RMSE results with rotmod-derived properties.
    Returns list of dicts with all properties for correlation analysis.
    """
    dat_files = sorted(glob.glob(os.path.join(sparc_dir, "**/*.dat"), recursive=True) +
                       glob.glob(os.path.join(sparc_dir, "*.dat")))

    # Build name → path map
    name_to_path = {}
    for path in dat_files:
        name = os.path.splitext(os.path.basename(path))[0]
        for suffix in ["_rotmod", "_Rotmod"]:
            name = name.replace(suffix, "")
        name_to_path[name] = path

    rows = []
    matched = not_found = prop_failed = 0

    for gal_name, result in results.items():
        path = name_to_path.get(gal_name)
        if path is None:
            not_found += 1
            continue

        props = extract_properties(path)
        if props is None:
            prop_failed += 1
            continue

        rows.append({
            "galaxy"           : gal_name,
            "rmse_kms"         : result["rmse_kms"],
            "v_flat_kms"       : result["v_flat_kms"],
            "g_bar_max_over_a0": result["g_bar_max_over_a0"],
            "n_points"         : result["n_points"],
            **props,
        })
        matched += 1

    print(
        f"[Analysis] Matched: {matched}  |  "
        f"Not found: {not_found}  |  "
        f"Parse failed: {prop_failed}"
    )
    return rows


# ---------------------------------------------------------------------------
# Spearman correlation table
# ---------------------------------------------------------------------------

def compute_correlations(rows: list[dict]) -> list[dict]:
    """Compute Spearman ρ and p-value for RMSE vs each property."""
    rmse = np.array([r["rmse_kms"] for r in rows])

    PROPERTIES = {
        "v_flat_kms"        : "V_flat [km/s]",
        "g_bar_max_over_a0" : "g_bar/a₀ (max)",
        "bulge_fraction"    : "Bulge fraction",
        "gas_fraction"      : "Gas fraction",
        "log_sb_central"    : "log SB_disk (central)",
        "mond_fraction"     : "MOND fraction",
        "profile_shape"     : "Profile shape",
        "n_points"          : "N data points",
    }

    correlations = []
    for key, label in PROPERTIES.items():
        x = np.array([r[key] for r in rows])
        rho, p = stats.spearmanr(x, rmse)
        correlations.append({
            "property": key,
            "label"   : label,
            "rho"     : float(rho),
            "p_value" : float(p),
            "abs_rho" : abs(float(rho)),
        })

    # Sort by |ρ| descending
    correlations.sort(key=lambda x: x["abs_rho"], reverse=True)

    print(f"\n[Analysis] Spearman correlations with RMSE (n={len(rows)} galaxies):")
    print(f"  {'Property':<30} {'ρ':>8} {'p':>10} {'|ρ|':>6}")
    print(f"  {'-'*58}")
    for c in correlations:
        sig = "***" if c["p_value"] < 0.001 else "**" if c["p_value"] < 0.01 else "*" if c["p_value"] < 0.05 else ""
        print(f"  {c['label']:<30} {c['rho']:>8.4f} {c['p_value']:>10.2e} {sig}")

    return correlations


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def plot_correlation_panels(rows: list[dict], correlations: list[dict], output_dir: str):
    """
    8-panel figure: RMSE vs each property, colour-coded by dynamical regime.
    """
    rmse    = np.array([r["rmse_kms"]          for r in rows])
    g_ratio = np.array([r["g_bar_max_over_a0"]  for r in rows])

    # Regime colour coding
    regime_colour = np.where(
        g_ratio < 0.3,  0,   # deep MOND — blue
        np.where(g_ratio < 5.0, 1, 2)  # transition — orange; Newtonian — red
    )
    colours = ["#4fc3f7", "#ffb74d", "#ef5350"]
    labels  = ["Deep MOND (g/a₀<0.3)", "Transition (0.3–5)", "Newtonian (g/a₀≥5)"]

    PANELS = [
        ("v_flat_kms",         "V_flat [km/s]",           False),
        ("g_bar_max_over_a0",  "g_bar/a₀ (max)",          True),
        ("bulge_fraction",     "Bulge fraction",           False),
        ("gas_fraction",       "Gas fraction",             False),
        ("log_sb_central",     "log₁₀ SBdisk (L☉/pc²)",  False),
        ("mond_fraction",      "MOND fraction of radii",   False),
        ("profile_shape",      "Profile shape index",      False),
        ("n_points",           "N data points",            False),
    ]

    fig = plt.figure(figsize=(20, 14), facecolor="#0a0a0a")
    gs  = gridspec.GridSpec(2, 4, figure=fig, hspace=0.42, wspace=0.32)

    # Lookup correlation for each panel
    rho_lookup = {c["property"]: (c["rho"], c["p_value"]) for c in correlations}

    for idx, (key, xlabel, logx) in enumerate(PANELS):
        row_idx = idx // 4
        col_idx = idx % 4
        ax = fig.add_subplot(gs[row_idx, col_idx])
        ax.set_facecolor("#111111")

        x = np.array([r[key] for r in rows])

        # Plot each regime
        for regime_id, colour, rlabel in zip([0, 1, 2], colours, labels):
            mask = regime_colour == regime_id
            if mask.sum() == 0:
                continue
            ax.scatter(
                x[mask], rmse[mask],
                c=colour, s=18, alpha=0.75, linewidths=0,
                label=f"{rlabel} (n={mask.sum()})",
            )

        ax.axhline(POOR_FIT_KMS, color="white", lw=1.0, ls="--", alpha=0.5)

        if logx:
            ax.set_xscale("log")

        rho, p = rho_lookup.get(key, (0.0, 1.0))
        sig = " ***" if p < 0.001 else " **" if p < 0.01 else " *" if p < 0.05 else ""

        ax.set_xlabel(xlabel, color="white", fontsize=9)
        ax.set_ylabel("RMSE [km/s]" if col_idx == 0 else "", color="white", fontsize=9)
        ax.set_title(
            f"ρ = {rho:+.3f}{sig}  (p={p:.2e})",
            color="white", fontsize=9,
        )
        ax.tick_params(colors="white", labelsize=8)
        ax.spines[["bottom","top","left","right"]].set_color("#444")
        ax.set_ylim(bottom=0)

        if idx == 0:
            leg = ax.legend(
                fontsize=6.5, facecolor="#1a1a1a", labelcolor="white",
                loc="upper left", markerscale=1.5,
            )

    fig.suptitle(
        "CODE-GEO V4.2 — RMSE vs Galaxy Properties  |  171 SPARC Galaxies\n"
        "Dashed line: RMSE = 20 km/s (poor-fit threshold)  |  "
        "*** p<0.001  ** p<0.01  * p<0.05",
        color="white", fontsize=11, fontweight="bold",
    )

    output_path = os.path.join(output_dir, "RMSE_PROPERTY_ANALYSIS.png")
    plt.savefig(output_path, facecolor="#0a0a0a", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n[Analysis] Figure → {output_path}")
    return output_path


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

def write_property_csv(rows: list[dict], output_dir: str):
    """Write full property table for reproducibility."""
    path = os.path.join(output_dir, "RMSE_PROPERTY_ANALYSIS.csv")
    fields = [
        "galaxy", "rmse_kms", "v_flat_kms", "g_bar_max_over_a0", "n_points",
        "bulge_fraction", "gas_fraction", "log_sb_central",
        "mond_fraction", "profile_shape", "inner_slope",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[Analysis] Property CSV → {path}")


# ---------------------------------------------------------------------------
# Paper-ready interpretation
# ---------------------------------------------------------------------------

def print_paper_summary(rows: list[dict], correlations: list[dict]):
    """Print the key numbers for Section 4.6 of the preprint."""
    rmse = np.array([r["rmse_kms"] for r in rows])

    poor  = rmse >= POOR_FIT_KMS
    n_poor = poor.sum()
    n_tot  = len(rows)

    print(f"\n{'='*60}")
    print(f"  Paper Summary — Section 4.6")
    print(f"{'='*60}")
    print(f"  Total galaxies : {n_tot}")
    print(f"  Poor fits      : {n_poor} ({100*n_poor/n_tot:.0f}%)")
    print(f"  Good fits      : {n_tot - n_poor} ({100*(n_tot-n_poor)/n_tot:.0f}%)")

    print(f"\n  Top 3 predictors of RMSE (Spearman |ρ|):")
    for c in correlations[:3]:
        print(
            f"    {c['label']:<30} ρ={c['rho']:+.3f}  "
            f"p={c['p_value']:.2e}"
        )

    # Mean properties of good vs poor fits
    good_mask = rmse < POOR_FIT_KMS
    poor_mask = ~good_mask

    print(f"\n  Mean properties: good fits vs poor fits")
    print(f"  {'Property':<25} {'Good (<20)':>12} {'Poor (≥20)':>12}")
    print(f"  {'-'*50}")
    for key, label, _ in [
        ("v_flat_kms",        "V_flat [km/s]",   False),
        ("bulge_fraction",    "Bulge fraction",   False),
        ("gas_fraction",      "Gas fraction",     False),
        ("g_bar_max_over_a0", "g_bar/a₀ max",    False),
        ("mond_fraction",     "MOND fraction",    False),
    ]:
        x = np.array([r[key] for r in rows])
        good_mean = np.mean(x[good_mask])
        poor_mean = np.mean(x[poor_mask])
        print(f"  {label:<25} {good_mean:>12.3f} {poor_mean:>12.3f}")

    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_analysis():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("🔬  CODE-GEO V4.2 — RMSE Property Correlation Analysis")
    print("=" * 55)

    # Load RMSE results
    if not os.path.exists(RESULTS_CSV):
        print(f"[ERROR] {RESULTS_CSV} not found. Run sparc_refinery_v4.py first.")
        return

    results = load_results(RESULTS_CSV)

    # Build property table
    rows = build_property_table(results, SPARC_DIR)
    if len(rows) < 10:
        print(f"[ERROR] Only {len(rows)} galaxies matched. Check SPARC data path.")
        return

    # Compute correlations
    correlations = compute_correlations(rows)

    # Write CSV
    write_property_csv(rows, OUTPUT_DIR)

    # Figure
    plot_correlation_panels(rows, correlations, OUTPUT_DIR)

    # Paper summary
    print_paper_summary(rows, correlations)


if __name__ == "__main__":
    run_analysis()