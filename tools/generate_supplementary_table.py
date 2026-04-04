"""
tools/generate_supplementary_table.py
======================================
CODE-GEO V4.2 — Supplementary Table 1 Generator
=================================================

Merges SPARC_RESULTS_V42.csv and RMSE_PROPERTY_ANALYSIS.csv into a
single, formatted supplementary table for the SPARC paper.

Outputs
-------
  docs/TABLE_S1.md       — Markdown version (for PREPRINT.md appendix)
  docs/TABLE_S1.tex      — LaTeX version (for journal submission)
  docs/TABLE_S1_full.csv — Merged full table (for reproducibility)

Columns in output
-----------------
  Galaxy name
  N_pts       — number of rotation curve data points
  RMSE        — km/s, unweighted, at λ_opt=0
  V_flat      — km/s, observed velocity at last data point
  g_bar/a₀    — peak baryonic acceleration in units of a₀
  f_bul       — bulge fraction at outer radius
  f_gas       — gas fraction at outer radius
  log SB      — log₁₀ central disk surface brightness [L⊙/pc²]
  f_MOND      — fraction of radii in MOND regime (g_bar/a₀ < 1)
  Flag        — ✓ excellent (<10), — acceptable (10–20), ✗ poor (≥20)

Sorting: by RMSE ascending (best fits first).

Usage
-----
  python3 tools/generate_supplementary_table.py
"""

import os
import sys
import csv
import math

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RESULTS_CSV   = "survey_outputs/SPARC_RESULTS_V42.csv"
PROPERTY_CSV  = "survey_outputs/RMSE_PROPERTY_ANALYSIS.csv"
OUTPUT_DIR    = "docs/"

POOR_FIT  = 20.0
GOOD_FIT  = 10.0


# ---------------------------------------------------------------------------
# Load and merge
# ---------------------------------------------------------------------------

def load_and_merge() -> list[dict]:
    """Load both CSVs and merge on galaxy name."""

    # Load RMSE results
    results = {}
    with open(RESULTS_CSV, newline="") as f:
        for row in csv.DictReader(f):
            results[row["galaxy"]] = row

    # Load property analysis
    props = {}
    if os.path.exists(PROPERTY_CSV):
        with open(PROPERTY_CSV, newline="") as f:
            for row in csv.DictReader(f):
                props[row["galaxy"]] = row

    # Merge
    merged = []
    for name, res in results.items():
        row = {
            "galaxy"           : name,
            "n_pts"            : int(res["n_points"]),
            "rmse_kms"         : float(res["rmse_kms"]),
            "v_flat_kms"       : float(res["v_flat_kms"]),
            "g_bar_max_over_a0": float(res["g_bar_max_over_a0"]),
        }
        if name in props:
            p = props[name]
            row["bulge_fraction" ] = float(p.get("bulge_fraction",  "nan"))
            row["gas_fraction"   ] = float(p.get("gas_fraction",    "nan"))
            row["log_sb_central" ] = float(p.get("log_sb_central",  "nan"))
            row["mond_fraction"  ] = float(p.get("mond_fraction",   "nan"))
        else:
            row["bulge_fraction" ] = float("nan")
            row["gas_fraction"   ] = float("nan")
            row["log_sb_central" ] = float("nan")
            row["mond_fraction"  ] = float("nan")

        merged.append(row)

    # Sort by RMSE ascending
    merged.sort(key=lambda r: r["rmse_kms"])
    return merged


def quality_flag(rmse: float) -> str:
    if rmse < GOOD_FIT:
        return "✓"
    elif rmse < POOR_FIT:
        return "—"
    else:
        return "✗"


def quality_flag_latex(rmse: float) -> str:
    if rmse < GOOD_FIT:
        return r"$\checkmark$"
    elif rmse < POOR_FIT:
        return r"$-$"
    else:
        return r"$\times$"


def fmt(val: float, decimals: int = 2, width: int = 0) -> str:
    """Format float; return '—' for NaN."""
    if math.isnan(val):
        return "—"
    s = f"{val:.{decimals}f}"
    return s.rjust(width) if width else s


# ---------------------------------------------------------------------------
# Markdown table
# ---------------------------------------------------------------------------

def write_markdown(rows: list[dict], path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)

    n_excellent = sum(1 for r in rows if r["rmse_kms"] < GOOD_FIT)
    n_good      = sum(1 for r in rows if GOOD_FIT <= r["rmse_kms"] < POOR_FIT)
    n_poor      = sum(1 for r in rows if r["rmse_kms"] >= POOR_FIT)

    with open(path, "w") as f:
        f.write("# Supplementary Table 1\n")
        f.write("## Per-Galaxy Rotation Curve Fit Results — CODE-GEO V4.2\n\n")
        f.write(
            "**Framework:** Mimetic-Conformal V4.2 | "
            "**Free function:** $F_\\text{MOND}(\\mathcal{Q}) = "
            "\\sqrt{\\mathcal{Q}(1+\\mathcal{Q})} - \\text{arcsinh}(\\sqrt{\\mathcal{Q}})$ | "
            "**Parameters:** $\\lambda = 0$ (optimised), "
            "$a_0 = 1.21 \\times 10^{-10}$ m s$^{-2}$ (pinned), "
            "$\\Upsilon_\\text{disk} = 0.50$, $\\Upsilon_\\text{bul} = 0.70$ "
            "(McGaugh & Schombert 2015)\n\n"
        )
        f.write(
            f"**Sample:** 171 galaxies from the SPARC database (Lelli et al. 2016). "
            f"4 galaxies excluded: fewer than 5 data points with positive velocity error.\n\n"
        )
        f.write(
            f"**Fit quality:** "
            f"✓ Excellent RMSE < 10 km s$^{{-1}}$: {n_excellent} galaxies | "
            f"— Acceptable 10–20 km s$^{{-1}}$: {n_good} galaxies | "
            f"✗ Poor ≥ 20 km s$^{{-1}}$: {n_poor} galaxies\n\n"
        )
        f.write(
            "Galaxies sorted by RMSE ascending. "
            "$f_\\text{bul}$ and $f_\\text{gas}$ computed at the outer 30% of radii. "
            "$\\log_{10}\\text{SB}_\\text{disk}$ is the central disk surface brightness "
            "[L$_\\odot$ pc$^{-2}$] from the median of the three innermost points. "
            "$f_\\text{MOND}$ is the fraction of radii where $g_\\text{bar}/a_0 < 1$.\n\n"
        )

        # Table header
        f.write("| Galaxy | $N$ | RMSE | $V_\\text{flat}$ | $g_\\text{bar}/a_0$ | "
                "$f_\\text{bul}$ | $f_\\text{gas}$ | $\\log\\text{SB}$ | $f_\\text{MOND}$ | Flag |\n")
        f.write("|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|\n")
        f.write("| | | [km/s] | [km/s] | (max) | | | [L⊙/pc²] | | |\n")
        f.write("|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|\n")

        for r in rows:
            f.write(
                f"| {r['galaxy']:<20} "
                f"| {r['n_pts']:>3} "
                f"| {fmt(r['rmse_kms'], 1):>6} "
                f"| {fmt(r['v_flat_kms'], 1):>7} "
                f"| {fmt(r['g_bar_max_over_a0'], 2):>7} "
                f"| {fmt(r['bulge_fraction'], 3):>7} "
                f"| {fmt(r['gas_fraction'], 3):>7} "
                f"| {fmt(r['log_sb_central'], 2):>6} "
                f"| {fmt(r['mond_fraction'], 3):>7} "
                f"| {quality_flag(r['rmse_kms']):^5} |\n"
            )

        # Summary statistics
        rmse_vals = [r["rmse_kms"] for r in rows]
        import statistics
        f.write(f"\n**Summary statistics (n = {len(rows)}):** "
                f"Median RMSE = {statistics.median(rmse_vals):.2f} km/s | "
                f"Mean RMSE = {sum(rmse_vals)/len(rmse_vals):.2f} km/s | "
                f"Min = {min(rmse_vals):.2f} km/s ({rows[0]['galaxy']}) | "
                f"Max = {max(rmse_vals):.2f} km/s ({rows[-1]['galaxy']})\n")

    print(f"[Table] Markdown → {path}  ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# LaTeX table
# ---------------------------------------------------------------------------

def write_latex(rows: list[dict], path: str):
    """
    Generate a journal-ready LaTeX longtable.
    Uses longtable for multi-page output.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)

    n_excellent = sum(1 for r in rows if r["rmse_kms"] < GOOD_FIT)
    n_good      = sum(1 for r in rows if GOOD_FIT <= r["rmse_kms"] < POOR_FIT)
    n_poor      = sum(1 for r in rows if r["rmse_kms"] >= POOR_FIT)

    with open(path, "w") as f:
        f.write("% Supplementary Table 1 — CODE-GEO V4.2 SPARC Results\n")
        f.write("% Generated by tools/generate_supplementary_table.py\n")
        f.write("% Requires: \\usepackage{longtable, booktabs, amssymb}\n\n")

        f.write("\\begin{longtable}{lrrrrrrrrr}\n")
        f.write("\\caption{Per-galaxy rotation curve fit results for the Mimetic-Conformal V4.2 "
                "framework ($\\lambda = 0$, $a_0 = 1.21 \\times 10^{-10}$ m\\,s$^{-2}$ pinned, "
                "$\\Upsilon_{\\rm disk} = 0.50$, $\\Upsilon_{\\rm bul} = 0.70$) "
                "applied to the SPARC sample of 171 disk galaxies "
                "(Lelli et al.~\\citeyear{Lelli2016}). "
                f"Fit quality: $\\checkmark$ excellent ($<10$ km\\,s$^{{-1}}$, {n_excellent} galaxies); "
                f"$-$ acceptable ($10$--$20$ km\\,s$^{{-1}}$, {n_good} galaxies); "
                f"$\\times$ poor ($\\geq 20$ km\\,s$^{{-1}}$, {n_poor} galaxies). "
                "Sorted by RMSE ascending.} "
                "\\label{tab:sparc_results}\\\\\n")

        # Header (first page)
        f.write("\\toprule\n")
        f.write("Galaxy & $N$ & RMSE & $V_{\\rm flat}$ & $g_{\\rm bar}/a_0$ & "
                "$f_{\\rm bul}$ & $f_{\\rm gas}$ & $\\log\\,{\\rm SB}$ & "
                "$f_{\\rm MOND}$ & Flag \\\\\n")
        f.write(" & & [km\\,s$^{-1}$] & [km\\,s$^{-1}$] & (max) & "
                " & & [$L_\\odot\\,{\\rm pc}^{-2}$] & & \\\\\n")
        f.write("\\midrule\n")
        f.write("\\endfirsthead\n\n")

        # Header (subsequent pages)
        f.write("\\multicolumn{10}{c}{\\tablename\\ \\thetable\\ -- continued from previous page}\\\\\n")
        f.write("\\toprule\n")
        f.write("Galaxy & $N$ & RMSE & $V_{\\rm flat}$ & $g_{\\rm bar}/a_0$ & "
                "$f_{\\rm bul}$ & $f_{\\rm gas}$ & $\\log\\,{\\rm SB}$ & "
                "$f_{\\rm MOND}$ & Flag \\\\\n")
        f.write("\\midrule\n")
        f.write("\\endhead\n\n")

        # Footer
        f.write("\\midrule\n")
        f.write("\\multicolumn{10}{r}{Continued on next page}\\\\\n")
        f.write("\\endfoot\n\n")
        f.write("\\bottomrule\n")
        f.write("\\endlastfoot\n\n")

        # Data rows — add quality separator rules
        prev_quality = None
        for r in rows:
            curr_quality = quality_flag(r["rmse_kms"])
            if prev_quality is not None and curr_quality != prev_quality:
                f.write("\\midrule\n")
            prev_quality = curr_quality

            # Escape underscores in galaxy names for LaTeX
            gal_name = r["galaxy"].replace("_", "\\_")

            f.write(
                f"{gal_name} & "
                f"{r['n_pts']} & "
                f"{fmt(r['rmse_kms'], 1)} & "
                f"{fmt(r['v_flat_kms'], 1)} & "
                f"{fmt(r['g_bar_max_over_a0'], 2)} & "
                f"{fmt(r['bulge_fraction'], 3)} & "
                f"{fmt(r['gas_fraction'], 3)} & "
                f"{fmt(r['log_sb_central'], 2)} & "
                f"{fmt(r['mond_fraction'], 3)} & "
                f"{quality_flag_latex(r['rmse_kms'])} \\\\\n"
            )

        f.write("\\end{longtable}\n")

    print(f"[Table] LaTeX    → {path}  ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# Full merged CSV
# ---------------------------------------------------------------------------

def write_full_csv(rows: list[dict], path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fields = [
        "galaxy", "n_pts", "rmse_kms", "v_flat_kms", "g_bar_max_over_a0",
        "bulge_fraction", "gas_fraction", "log_sb_central", "mond_fraction",
        "quality_flag",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            r["quality_flag"] = quality_flag(r["rmse_kms"])
            writer.writerow({k: (fmt(r[k], 4) if isinstance(r.get(k), float)
                                 else r.get(k, "")) for k in fields})
    print(f"[Table] Full CSV → {path}  ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# Console summary for paper text
# ---------------------------------------------------------------------------

def print_summary(rows: list[dict]):
    import statistics
    rmse_vals = [r["rmse_kms"] for r in rows]

    print(f"\n{'='*60}")
    print(f"  Supplementary Table 1 — Summary")
    print(f"{'='*60}")
    print(f"  n          : {len(rows)}")
    print(f"  Median RMSE: {statistics.median(rmse_vals):.2f} km/s")
    print(f"  Mean RMSE  : {sum(rmse_vals)/len(rmse_vals):.2f} km/s")
    print(f"  Best fit   : {rows[0]['galaxy']} ({rows[0]['rmse_kms']:.2f} km/s)")
    print(f"  Worst fit  : {rows[-1]['galaxy']} ({rows[-1]['rmse_kms']:.2f} km/s)")

    # Quality breakdown with thresholds
    thresholds = [5, 10, 15, 20, 30, 50]
    print(f"\n  Cumulative fraction by RMSE threshold:")
    for t in thresholds:
        n = sum(1 for r in rows if r["rmse_kms"] < t)
        print(f"    RMSE < {t:>3} km/s : {n:>3}/{len(rows)} ({100*n/len(rows):.0f}%)")

    # Best and worst 5
    print(f"\n  5 best fits:")
    for r in rows[:5]:
        print(f"    {r['galaxy']:<20} {r['rmse_kms']:>6.2f} km/s")
    print(f"\n  5 worst fits:")
    for r in rows[-5:]:
        print(f"    {r['galaxy']:<20} {r['rmse_kms']:>6.2f} km/s")
    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run():
    print("📋  CODE-GEO V4.2 — Supplementary Table 1 Generator")
    print("=" * 55)

    # Check inputs
    for path in [RESULTS_CSV, PROPERTY_CSV]:
        if not os.path.exists(path):
            print(f"[ERROR] Missing: {path}")
            print(f"  Run sparc_refinery_v4.py and rmse_property_analysis.py first.")
            return

    rows = load_and_merge()
    print(f"[Table] {len(rows)} galaxies loaded and merged.")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    write_markdown(rows, os.path.join(OUTPUT_DIR, "TABLE_S1.md"))
    write_latex(rows,    os.path.join(OUTPUT_DIR, "TABLE_S1.tex"))
    write_full_csv(rows, os.path.join(OUTPUT_DIR, "TABLE_S1_full.csv"))

    print_summary(rows)

    print(f"\n  Files written to {OUTPUT_DIR}")
    print(f"  TABLE_S1.md  — paste into PREPRINT.md appendix")
    print(f"  TABLE_S1.tex — include in LaTeX manuscript with \\input{{TABLE_S1}}")
    print(f"  TABLE_S1_full.csv — archive with paper data release")


if __name__ == "__main__":
    run()