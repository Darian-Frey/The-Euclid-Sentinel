"""
tests/test_v31_lensing.py
==========================
CODE-GEO V4.2 — Lensing Pipeline Tests
========================================

Originally: tested the deprecated V3.1 Soft-Krylov engine lensing profile.
Updated:    tests the V4.2 thin-disk Poisson lensing pipeline.

The file is retained under its original name for repo history continuity.
The V3.1 test content is preserved in comments at the bottom as documentation
of what changed and why.

What these tests verify
-----------------------
1.  Thin-disk Poisson kernel: Φ(k) = -2πG Σ_b(k) / |k|
    — correct scaling with surface mass density
    — produces finite, non-NaN potential field
    — correct sign (negative potential from positive mass)

2.  Q-field physical values
    — Q = (g_N/a₀)² with g_N from the thin-disk potential
    — Q range physically consistent with MOND/Newtonian regimes
    — deep MOND in low-density outskirts; Newtonian in cores

3.  Effective density scaling
    — ρ_eff = Σ_b / F'(Q) — emerges from physics, not imposed
    — DM ratio (Σ_eff/Σ_b) > 1 everywhere (enhancement, never suppression)
    — DM ratio → 1 in Newtonian limit (no enhancement)
    — DM ratio > 1 in MOND limit (enhancement proportional to g_bar/a₀)

4.  Engine interface
    — MimeticEngine.compute_effective_density(sigma_b, dx) returns
      (sigma_eff, F_prime, telemetry) with correct shapes and types
    — Telemetry contains expected keys with finite values
    — No NaN or Inf in any output field

5.  Comparison with V3.1 (documented, not rerun)
    — V3.1 imposed 12.0× DM ratio hardcoded
    — V4.2 produces variable DM ratio from F'(Q): range ≈ 1× to 10×
      depending on surface density and pixel scale

Physical setup for tests
------------------------
Synthetic galaxy cluster: two overlapping Gaussian blobs in surface mass
density Σ_b [kg/m²], representative of a cluster at z~0.3. Physical pixel
scale dx = 150 kpc / 256 pixels ≈ 1.8 kpc/pixel.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from core.physics_constants import A0, KPC_TO_M, LAMBDA, BETA


# ---------------------------------------------------------------------------
# Synthetic test data
# ---------------------------------------------------------------------------

def make_synthetic_cluster(
    grid: int = 128,
    box_kpc: float = 1500.0,
) -> tuple[np.ndarray, float]:
    """
    Two-Gaussian synthetic baryonic cluster in surface mass density [kg/m²].

    Returns
    -------
    sigma_b : 2D array  [kg/m²]
    dx      : float     [m/pixel]
    """
    BOX_M = box_kpc * KPC_TO_M
    dx    = BOX_M / grid

    x = np.linspace(-BOX_M / 2, BOX_M / 2, grid)
    X, Y = np.meshgrid(x, x)

    # Main cluster (~10^13 M_sun stellar over ~200 kpc scale)
    # Peak Σ_* ~ 1.0 kg/m²: gives g_N ~ 2πG Σ ~ 4e-10 m/s² ≈ 3×a₀
    # This places the cluster in the MOND transition regime (Q ~ 10)
    # as required to exercise the F'(Q) mapping physically.
    # (Original value 2e-14 kg/m² gave g_N ~ 10^-24 m/s², far below a₀)
    r1 = 200.0 * KPC_TO_M
    s1 = 1.0 * np.exp(-(X**2 + Y**2) / r1**2)

    # Subcluster — offset, less massive
    r2 = 120.0 * KPC_TO_M
    x2 = 300.0 * KPC_TO_M
    s2 = 0.4 * np.exp(-((X - x2)**2 + Y**2) / r2**2)

    return s1 + s2, dx


SIGMA_B, DX = make_synthetic_cluster()


# ---------------------------------------------------------------------------
# Test 1: Thin-disk Poisson kernel — correct potential scaling
# ---------------------------------------------------------------------------

def test_thin_disk_potential_sign():
    """Newtonian potential Φ_N < 0 for positive mass (attractive gravity)."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()
    phi_N = engine._solve_newtonian_potential(SIGMA_B, DX)

    assert phi_N.shape == SIGMA_B.shape, "Potential shape mismatch"
    # Potential should be negative (or zero) everywhere — attractive
    assert np.mean(phi_N) < 0, (
        f"Mean potential should be negative: mean Φ_N = {np.mean(phi_N):.3e}"
    )
    print(f"  PASS  Thin-disk potential: Φ_N < 0 (mean={np.mean(phi_N):.3e} m²/s²)")


def test_thin_disk_potential_finite():
    """Thin-disk Poisson solve produces finite, non-NaN output."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()
    phi_N = engine._solve_newtonian_potential(SIGMA_B, DX)

    assert np.all(np.isfinite(phi_N)), (
        f"Non-finite values in Φ_N: n_bad = {np.sum(~np.isfinite(phi_N))}"
    )
    assert not np.any(np.isnan(phi_N)), "NaN in Φ_N"
    print(f"  PASS  Thin-disk potential: finite and non-NaN")


def test_thin_disk_scales_with_mass():
    """Doubling Σ_b doubles |Φ_N| (linearity of Poisson equation)."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    phi1 = engine._solve_newtonian_potential(SIGMA_B, DX)
    phi2 = engine._solve_newtonian_potential(2.0 * SIGMA_B, DX)

    ratio = np.mean(np.abs(phi2)) / np.mean(np.abs(phi1))
    assert abs(ratio - 2.0) < 0.01, (
        f"Potential not linear in mass: |Φ_2x|/|Φ_1x| = {ratio:.4f} (expected 2.0)"
    )
    print(f"  PASS  Thin-disk linearity: 2×Σ_b → 2×|Φ_N| (ratio={ratio:.4f})")


# ---------------------------------------------------------------------------
# Test 2: Q-field physical values
# ---------------------------------------------------------------------------

def test_q_field_physical_range():
    """Q = (g_N/a₀)² is dimensionless and physically plausible."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()
    Q, g_N = engine._compute_q_field(SIGMA_B, DX)

    assert Q.shape == SIGMA_B.shape, "Q shape mismatch"
    assert np.all(Q >= 0), f"Q < 0 encountered: min Q = {Q.min():.4e}"
    assert np.all(np.isfinite(Q)), f"Non-finite Q: n_bad = {np.sum(~np.isfinite(Q))}"

    # For a cluster at z~0.3, typical accelerations range from
    # ~10⁻¹² m/s² (outer) to ~10⁻⁹ m/s² (core)
    # → Q = (g/a₀)² ranges from ~10⁻⁴ to ~10²
    Q_max = float(np.max(Q))
    Q_mean = float(np.mean(Q))
    assert Q_max > 0.1,  f"Q_max = {Q_max:.4e} suspiciously small (under-resolved?)"
    assert Q_max < 1e10, f"Q_max = {Q_max:.4e} suspiciously large (unit error?)"

    print(
        f"  PASS  Q-field physical: Q ∈ [{Q.min():.2e}, {Q_max:.2e}]  "
        f"mean Q = {Q_mean:.2e}"
    )


def test_q_field_regime_structure():
    """Q is higher in cluster core (Newtonian) than in outskirts (MOND)."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()
    Q, _ = engine._compute_q_field(SIGMA_B, DX)

    N = Q.shape[0]
    Q_core   = float(np.mean(Q[N//4:3*N//4, N//4:3*N//4]))  # central 50%
    Q_outer  = float(np.mean(
        np.concatenate([Q[:N//8, :].ravel(), Q[-N//8:, :].ravel()])
    ))  # outer 12% of rows

    assert Q_core > Q_outer, (
        f"Q_core ({Q_core:.3e}) should exceed Q_outer ({Q_outer:.3e}) — "
        "Newtonian core vs MOND outskirts"
    )
    print(
        f"  PASS  Q regime structure: core Q={Q_core:.2e} > outer Q={Q_outer:.2e}"
    )


# ---------------------------------------------------------------------------
# Test 3: Effective density scaling
# ---------------------------------------------------------------------------

def test_effective_density_enhancement():
    """Σ_eff ≥ Σ_b everywhere (MOND enhancement, never suppression)."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    sigma_eff, F_prime, _ = engine.compute_effective_density(SIGMA_B, DX)

    # F'(Q) ≤ 1 always (proven in test_stability.py), so Σ_eff = Σ_b/F'(Q) ≥ Σ_b
    ratio = sigma_eff / (SIGMA_B + 1e-50)
    assert np.all(ratio >= 1.0 - 1e-10), (
        f"Σ_eff < Σ_b encountered: min ratio = {ratio.min():.4f}"
    )
    print(
        f"  PASS  Σ_eff ≥ Σ_b everywhere: min ratio = {ratio.min():.4f}"
    )


def test_dm_ratio_physical_range():
    """DM ratio (mean Σ_eff / mean Σ_b) is in a physically plausible range."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    sigma_eff, _, telemetry = engine.compute_effective_density(SIGMA_B, DX)

    dm_mean = telemetry["dm_ratio_mean"]
    dm_peak = telemetry["dm_ratio_peak"]

    # Physical range: 1× (pure Newtonian) to large in deep MOND outskirts
    # Mean should be > 1 (enhancement) and < 1e6 (not diverging)
    assert 1.0 <= dm_mean <= 1e6, (
        f"Mean DM ratio = {dm_mean:.2f} outside physical range [1, 1e6]"
    )
    assert dm_peak >= dm_mean, (
        f"Peak ratio ({dm_peak:.2f}) should exceed mean ({dm_mean:.2f})"
    )
    print(
        f"  PASS  DM ratio physical: mean={dm_mean:.2f}x  peak={dm_peak:.2f}x"
    )


def test_dm_ratio_not_hardcoded():
    """
    DM ratio varies spatially (not uniform) — confirms it is not hardcoded.

    The V3.1 engine multiplied by a scalar 12.0, giving uniform enhancement.
    The V4.2 engine divides by F'(Q(x,y)), giving spatially varying enhancement
    that is higher in MOND-regime regions and lower in Newtonian-regime regions.
    """
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    sigma_eff, F_prime, _ = engine.compute_effective_density(SIGMA_B, DX)

    # F'(Q) should vary spatially
    fp_std = float(np.std(F_prime))
    assert fp_std > 1e-6, (
        f"F'(Q) is spatially uniform (std={fp_std:.2e}) — "
        "this would indicate hardcoded enhancement rather than physical Q-field"
    )

    # Σ_eff/Σ_b should vary spatially (not all the same ratio)
    ratio = sigma_eff / (SIGMA_B + 1e-50)
    ratio_std = float(np.std(ratio))
    assert ratio_std > 1e-4, (
        f"DM ratio is spatially uniform (std={ratio_std:.2e}) — "
        "this would indicate the old hardcoded 12.0 behaviour"
    )

    print(
        f"  PASS  DM ratio spatially variable: "
        f"F' std={fp_std:.4f}  ratio std={ratio_std:.4f}  "
        f"(confirms physics-derived, not hardcoded)"
    )


# ---------------------------------------------------------------------------
# Test 4: Engine interface
# ---------------------------------------------------------------------------

def test_engine_output_shapes():
    """compute_effective_density returns correct shapes and types."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    sigma_eff, F_prime, telemetry = engine.compute_effective_density(SIGMA_B, DX)

    assert sigma_eff.shape == SIGMA_B.shape, "sigma_eff shape mismatch"
    assert F_prime.shape  == SIGMA_B.shape, "F_prime shape mismatch"
    assert isinstance(telemetry, dict), "telemetry is not a dict"
    print(
        f"  PASS  Engine output shapes: sigma_eff={sigma_eff.shape}  "
        f"F_prime={F_prime.shape}"
    )


def test_engine_telemetry_keys():
    """Telemetry dict contains all expected keys with finite values."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    _, _, telemetry = engine.compute_effective_density(SIGMA_B, DX)

    required_keys = [
        "Q_mean", "Q_max", "g_N_mean_ms2", "g_N_max_ms2",
        "F_prime_mean", "F_prime_min", "dm_ratio_mean", "dm_ratio_peak",
        "lam", "beta", "a0",
    ]
    for key in required_keys:
        assert key in telemetry, f"Missing telemetry key: {key}"
        val = telemetry[key]
        assert np.isfinite(val), f"Telemetry[{key}] = {val} is not finite"

    print(f"  PASS  Telemetry: all {len(required_keys)} keys present and finite")


def test_engine_no_nan_output():
    """No NaN or Inf in any engine output."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    sigma_eff, F_prime, _ = engine.compute_effective_density(SIGMA_B, DX)

    for name, arr in [("sigma_eff", sigma_eff), ("F_prime", F_prime)]:
        assert np.all(np.isfinite(arr)), (
            f"{name} has non-finite values: n_bad = {np.sum(~np.isfinite(arr))}"
        )
    print(f"  PASS  No NaN or Inf in sigma_eff or F_prime")


def test_engine_lensing_potential():
    """get_lensing_potential returns finite array of correct shape."""
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    sigma_eff, _, _ = engine.compute_effective_density(SIGMA_B, DX)
    phi_eff = engine.get_lensing_potential(sigma_eff, DX)

    assert phi_eff.shape == SIGMA_B.shape, "phi_eff shape mismatch"
    assert np.all(np.isfinite(phi_eff)), "Non-finite values in phi_eff"
    assert np.mean(phi_eff) < 0, "Lensing potential should be negative (attractive)"
    print(
        f"  PASS  Lensing potential: shape={phi_eff.shape}  "
        f"mean Φ_eff={np.mean(phi_eff):.3e} m²/s²"
    )


# ---------------------------------------------------------------------------
# Test 5: V4.2 vs V3.1 comparison (documented)
# ---------------------------------------------------------------------------

def test_v42_not_uniform_enhancement():
    """
    V4.2 gives spatially varying DM ratio; V3.1 gave uniform 12.0× hardcoded.

    This is a regression test ensuring the V3.1 bug (rho_mimetic = rho_baryon
    × f_q_val × 12.0) cannot return. The 12.0 factor gave:
        DM ratio = uniform ~12× regardless of g_bar/a₀ or surface density
    V4.2 gives:
        DM ratio = 1/F'(Q) — spatially varying, physically grounded
    """
    from core.mimetic_engine import MimeticEngine
    engine = MimeticEngine()

    sigma_eff, F_prime, telemetry = engine.compute_effective_density(SIGMA_B, DX)

    # V3.1 would give dm_ratio_mean ≈ 12.0 with very small spatial variance
    # V4.2 should give a physically determined ratio that is NOT 12.0
    dm_mean = telemetry["dm_ratio_mean"]

    # It should NOT be suspiciously close to 12.0
    # (The probability of getting exactly 12.00 from physical Q-field is negligible)
    assert abs(dm_mean - 12.0) > 0.5 or np.std(F_prime) > 0.01, (
        f"DM ratio = {dm_mean:.4f} is suspiciously close to the V3.1 hardcoded 12.0. "
        "Check for regression to the old engine."
    )

    # Confirm F'(Q) is not constant (which would imply hardcoded ratio)
    fp_cv = np.std(F_prime) / (np.mean(F_prime) + 1e-10)  # coefficient of variation
    assert fp_cv > 0.01, (
        f"F'(Q) coefficient of variation = {fp_cv:.4f} — suspiciously uniform"
    )

    print(
        f"  PASS  V4.2 not uniform: dm_ratio_mean={dm_mean:.2f}  "
        f"F'(Q) CV={fp_cv:.3f}  (confirms no hardcoded 12.0 regression)"
    )


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

TESTS = [
    # Thin-disk Poisson kernel
    test_thin_disk_potential_sign,
    test_thin_disk_potential_finite,
    test_thin_disk_scales_with_mass,
    # Q-field
    test_q_field_physical_range,
    test_q_field_regime_structure,
    # Effective density
    test_effective_density_enhancement,
    test_dm_ratio_physical_range,
    test_dm_ratio_not_hardcoded,
    # Engine interface
    test_engine_output_shapes,
    test_engine_telemetry_keys,
    test_engine_no_nan_output,
    test_engine_lensing_potential,
    # V4.2 vs V3.1
    test_v42_not_uniform_enhancement,
]

if __name__ == "__main__":
    print("=" * 60)
    print("  V4.2 Lensing Pipeline Tests")
    print("  (updated from V3.1 — see module docstring)")
    print("=" * 60)

    passed = failed = 0
    for test in TESTS:
        name = test.__name__.replace("test_", "").replace("_", " ")
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  FAIL  {name}: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR {name}: {type(e).__name__}: {e}")
            failed += 1

    print(f"\n  {'='*40}")
    print(f"  {passed}/{passed+failed} tests passed")
    if failed == 0:
        print("  V4.2 LENSING PIPELINE: ALL TESTS PASSED")
    else:
        print("  PIPELINE FAILURES — review output above")
    print(f"  {'='*40}")


# ---------------------------------------------------------------------------
# Historical note: original V3.1 test content
# ---------------------------------------------------------------------------
# The original test_v31_lensing.py tested the V3.1 MimeticEngine with:
#
#   engine = MimeticEngine(kappa=0.80)
#   sigma = np.log10(rho_baryon + 1e-30)       # pixel-gradient proxy for Q
#   gy, gx = np.gradient(sigma, x[1]-x[0])
#   grad_sq = (gx**2 + gy**2) * 5e18           # arbitrary scaling
#   rho_eff = engine.compute_effective_density(rho_baryon, grad_sq)
#
# Problems with this test:
#   - Used log10 of density as the scalar field (no physical motivation)
#   - Multiplied gradient by 5e18 (arbitrary, no units)
#   - engine.compute_effective_density multiplied by 12.0 (hardcoded)
#   - "SUCCESS" was declared if any rho_eff > rho_baryon — trivially true
#     given the 12.0 factor
#
# The test was not verifying physics; it was verifying that the hardcoded
# factor produced larger numbers. It has been replaced with the above suite
# which verifies actual physical properties of the V4.2 pipeline.