#!/usr/bin/env python3
"""
cleanup.py — Euclid Sentinel repo cleanup script
=================================================

Run from the project root:
    python3 cleanup.py

What this does
--------------
1. Deletes stray files accidentally created in the project root
   (np, os, sys, 0.01, — created by shell import accidents)
2. Moves old Phase VIII PNG outputs to archive/phase_08/
   (generated with hardcoded 12.0 DM factor — scientifically invalid)
3. Archives old SUMMARY_DATA.csv
4. Reports remaining items that need manual attention

What this does NOT do
---------------------
- Does not touch data/ (SPARC .dat files, HST FITS) — too large to move
- Does not modify any .py source files
- Does not delete anything without printing what it's deleting

Run with --dry-run to preview without making changes.
"""

import os
import sys
import shutil
import argparse
from pathlib import Path

ROOT = Path(__file__).parent


def ensure_archive():
    archive = ROOT / "archive" / "phase_08"
    archive.mkdir(parents=True, exist_ok=True)
    return archive


def cleanup(dry_run: bool = False):
    tag = "[DRY RUN] " if dry_run else ""
    archive = ensure_archive()
    deleted  = []
    archived = []
    skipped  = []

    # ── 1. Stray root-level files ────────────────────────────────────────
    STRAY_FILES = ["np", "os", "sys", "0.01,"]
    print("\n── Stray files ─────────────────────────────────────────────────")
    for name in STRAY_FILES:
        path = ROOT / name
        if path.exists() and path.is_file():
            print(f"  {tag}DELETE  {path}")
            if not dry_run:
                path.unlink()
            deleted.append(str(path))
        else:
            print(f"  SKIP    {path}  (not found)")

    # ── 2. Old Phase VIII PNG outputs (root level) ───────────────────────
    PHASE08_PNGS = [
        "SENTINEL_REPORT_Abell_370.png",
        "SENTINEL_REPORT_Bullet_Cluster.png",
        "SENTINEL_REPORT_El_Gordo.png",
        "SENTINEL_REPORT_HUDF_DeepField.png",
        "Sentinel_FINAL_REPORT.png",
    ]
    print("\n── Phase VIII artefacts → archive/phase_08/ ────────────────────")
    for name in PHASE08_PNGS:
        src = ROOT / name
        dst = archive / name
        if src.exists():
            print(f"  {tag}ARCHIVE {src.name}")
            if not dry_run:
                shutil.move(str(src), str(dst))
            archived.append(name)
        else:
            print(f"  SKIP    {name}  (not found)")

    # ── 3. Old SUMMARY_DATA.csv (root level) ─────────────────────────────
    print("\n── Old SUMMARY_DATA.csv ─────────────────────────────────────────")
    old_csv = ROOT / "SUMMARY_DATA.csv"
    if old_csv.exists():
        dst = archive / "SUMMARY_DATA_phase08.csv"
        print(f"  {tag}ARCHIVE {old_csv.name}  → archive/phase_08/SUMMARY_DATA_phase08.csv")
        print(f"           (contains hardcoded 12.0× DM ratio artefacts)")
        if not dry_run:
            shutil.move(str(old_csv), str(dst))
        archived.append("SUMMARY_DATA.csv")
    else:
        print(f"  SKIP    SUMMARY_DATA.csv  (not found)")

    # ── 4. sentinel_engine.log ────────────────────────────────────────────
    print("\n── Runtime logs ─────────────────────────────────────────────────")
    log = ROOT / "sentinel_engine.log"
    if log.exists():
        print(f"  {tag}DELETE  sentinel_engine.log  (runtime log, gitignored)")
        if not dry_run:
            log.unlink()
        deleted.append("sentinel_engine.log")
    else:
        print(f"  SKIP    sentinel_engine.log  (not found)")

    # ── 5. Items needing manual review ────────────────────────────────────
    print("\n── Manual review needed ─────────────────────────────────────────")

    # core/physics directory (unexpected)
    phys_dir = ROOT / "core" / "physics"
    if phys_dir.exists():
        print(f"  REVIEW  core/physics/  — unexpected directory, check contents")
        skipped.append("core/physics/")

    # tools/data_processing directory
    dp_dir = ROOT / "tools" / "data_processing"
    if dp_dir.exists():
        contents = list(dp_dir.iterdir())
        print(f"  REVIEW  tools/data_processing/  — {len(contents)} items")
        skipped.append("tools/data_processing/")

    # simulations/results directory
    sim_res = ROOT / "simulations" / "results"
    if sim_res.exists():
        contents = list(sim_res.iterdir())
        if not contents:
            print(f"  {tag}DELETE  simulations/results/  — empty directory")
            if not dry_run:
                sim_res.rmdir()
        else:
            print(f"  REVIEW  simulations/results/  — {len(contents)} items")
            skipped.append("simulations/results/")

    # ── 6. Empty source files ─────────────────────────────────────────────
    print("\n── Empty source files (check stubs applied) ─────────────────────")
    CHECK_STUBS = [
        "core/conformal_constraints.py",
        "core/solvers.py",
        "tools/lensing_audit.py",
        "tools/lensing_sim_v31.py",
        "tools/euclid_mock_gen.py",
        "tools/sentinel_survey.py",
    ]
    for rel in CHECK_STUBS:
        p = ROOT / rel
        if p.exists():
            size = p.stat().st_size
            status = "stubbed" if size > 100 else "EMPTY — add stub comment"
            print(f"  {status:<30s}  {rel}")
        else:
            print(f"  MISSING  {rel}")

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n{'='*55}")
    print(f"  Deleted  : {len(deleted)} files")
    print(f"  Archived : {len(archived)} files → archive/phase_08/")
    print(f"  Review   : {len(skipped)} items need manual inspection")
    if dry_run:
        print(f"\n  DRY RUN — no changes made. Remove --dry-run to apply.")
    else:
        print(f"\n  Cleanup complete.")
    print(f"{'='*55}")

    print("\nNext steps:")
    print("  1. Copy .gitignore to project root")
    print("  2. git add .gitignore requirements.txt MANIFEST.json docs/")
    print("  3. git rm --cached data/sparc/*.dat data/sparc/*.zip  (if tracked)")
    print("  4. git rm --cached survey_outputs/  (if tracked)")
    print("  5. git commit -m 'chore: repo cleanup, V4.2 gitignore, updated manifest'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Euclid Sentinel repo cleanup")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview changes without applying them",
    )
    args = parser.parse_args()
    cleanup(dry_run=args.dry_run)
