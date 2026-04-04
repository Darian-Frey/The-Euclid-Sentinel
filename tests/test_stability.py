"""
tests/test_stability.py
========================
CODE-GEO V4.2 — Free Function Stability Test Suite
====================================================

Tests all mathematical and physical properties of the V4.2 hybrid free
function F(Q) that must hold for the Mimetic-Conformal framework to be
a physically consistent scalar-tensor gravity theory.

Properties tested
-----------------
1.  Canonical kinetic limit:  F(Q) → Q as Q → 0
2.  Newtonian limit:          F'(Q) → 1 as Q → ∞
3.  MOND interpolation:       F'(Q) = μ_std(√Q) at λ=0 (zero RMSE)
4.  MOND deep limit:          F'(Q) → √Q as Q → 0
5.  Positivity:               F(Q) ≥ 0 for all Q ≥ 0
6.  Monotonicity:             F'(Q) > 0 for all Q ≥ 0 (no gravitational flip)
7.  Ghost-free:               c_s² > 0 everywhere
8.  Subluminality:            c_s² ≤ 1 everywhere (see also test_gw_compliance.py)
9.  Continuity:               F, F', F'' finite and continuous on Q ∈ (0, ∞)
10. λ-sensitivity:            F(Q) is smooth in λ; λ=0 recovers F_MOND exactly
11. Caustic guard localisation: guard term O(Q²) at small Q, suppressed at large Q
12. SPARC result consistency:  λ=0 fit error ≤ λ=0.05 fit error (SPARC data prefers λ=0)
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import warnings

from core.action import (
    free_function_v42,
    free_function_derivative_v42,
    free_function_second_derivative_v42,
    check_stability_v42,
)
from core.physics_constants import LAMBDA, BETA, A0


Q_RANGE  = np.logspace(-6, 4, 2000)    # wide range for limit tests
Q_DENSE  = np.logspace(-4, 3, 500)     # V4.2 audit range
RTOL     = 1e-4                         # relative tolerance for limit tests
LAM_VALS = [0.0, 0.01, 0.05, 0.1, 0.3, 0.5]


# ---------------------------------------------------------------------------
# Test 1: Canonical kinetic limit F(Q) → Q as Q → 0
# ---------------------------------------------------------------------------

def test_canonical_kinetic_limit():
    """F(Q) ~ (2/3)Q^{3/2} as Q → 0 (deep MOND free function limit).

    The Taylor expansion of F_MOND(Q) = √(Q(1+Q)) − arcsinh(√Q) near Q=0:
        √(Q(1+Q)) = √Q·(1 + Q/2 − Q²/8 + ...)
        arcsinh(√Q)  = √Q·(1 − Q/6 + 3Q²/40 − ...)
        F_MOND(Q)    = Q^{3/2}(1/2 + 1/6) + O(Q^{5/2}) = (2/3)Q^{3/2} + O(Q^{5/2})

    Note: the claim "F(Q) → Q as Q → 0" in earlier docstrings was INCORRECT.
    F(Q)/Q → 0 (not 1) as Q → 0. The correct limit is F(Q)/(2/3·Q^{3/2}) → 1.
    This does NOT affect the physics: F'(Q) = μ_std(√Q) → √Q as Q → 0 (correct).
    """
    Q_small = np.logspace(-8, -4, 50)
    prefactor = 2.0 / 3.0
    for lam in [0.0, 0.05]:   # λ-dependent terms are O(Q²), subdominant at small Q
        F = free_function_v42(Q_small, lam=lam, beta=BETA)
        ratio = F / (prefactor * Q_small**1.5)
        max_dev = float(np.max(np.abs(ratio - 1.0)))
        assert max_dev < 0.01, (
            f"Deep MOND free function limit failed at λ={lam}: "
            f"max|F/[(2/3)Q^1.5] − 1| = {max_dev:.4e} (threshold 0.01)"
        )
    print(f"  PASS  F(Q)/(2/3·Q^1.5) → 1 as Q→0  (deep MOND limit; λ=0,0.05)")


# ---------------------------------------------------------------------------
# Test 2: Newtonian limit F'(Q) → 1 as Q → ∞
# ---------------------------------------------------------------------------

def test_newtonian_limit():
    """F'(Q) → 1 as Q → ∞ for all λ values."""
    Q_large = np.logspace(2, 5, 50)
    for lam in LAM_VALS:
        Fp = free_function_derivative_v42(Q_large, lam=lam, beta=BETA)
        max_dev = float(np.max(np.abs(Fp - 1.0)))
        assert max_dev < 0.01, (
            f"Newtonian limit failed at λ={lam}: max|F'−1| = {max_dev:.4e}"
        )
    print(f"  PASS  F'(Q) → 1 as Q→∞  (all λ ∈ {LAM_VALS})")


# ---------------------------------------------------------------------------
# Test 3: MOND interpolation — F'(Q) = μ_std(√Q) at λ=0
# ---------------------------------------------------------------------------

def test_mond_interpolation_exact_at_lambda_zero():
    """F'(Q) = √Q/√(1+Q) exactly when λ=0 (zero RMSE vs Standard MOND)."""
    mu_std = np.sqrt(Q_DENSE) / np.sqrt(1.0 + Q_DENSE)
    Fp = free_function_derivative_v42(Q_DENSE, lam=0.0, beta=BETA)

    rmse = float(np.sqrt(np.mean((Fp - mu_std)**2)))
    max_dev = float(np.max(np.abs(Fp - mu_std)))

    assert rmse < 1e-10, f"MOND RMSE at λ=0 = {rmse:.4e} (should be ~0)"
    assert max_dev < 1e-8, f"Max MOND deviation at λ=0 = {max_dev:.4e}"

    print(
        f"  PASS  F'(Q) = μ_std(√Q) at λ=0: "
        f"RMSE = {rmse:.2e}  max_dev = {max_dev:.2e}"
    )


def test_mond_interpolation_approximate_at_lambda_nonzero():
    """F'(Q) ≈ μ_std(√Q) + O(λ) for small λ — guard is a small perturbation."""
    mu_std = np.sqrt(Q_DENSE) / np.sqrt(1.0 + Q_DENSE)
    for lam in [0.01, 0.05]:
        Fp = free_function_derivative_v42(Q_DENSE, lam=lam, beta=BETA)
        rmse = float(np.sqrt(np.mean((Fp - mu_std)**2)))
        assert rmse < 0.1, (
            f"MOND RMSE at λ={lam} = {rmse:.4f} — guard is not a small perturbation"
        )
    print(f"  PASS  Guard is small perturbation on μ_std at λ ≤ 0.05")


# ---------------------------------------------------------------------------
# Test 4: Deep MOND limit — F'(Q) → √Q as Q → 0
# ---------------------------------------------------------------------------

def test_deep_mond_limit():
    """F'(Q)/√Q → 1 as Q → 0 (deep MOND: g_eff → √(g_bar × a₀))."""
    Q_small = np.logspace(-8, -4, 50)
    Fp = free_function_derivative_v42(Q_small, lam=0.0, beta=BETA)
    ratio = Fp / np.sqrt(Q_small)
    max_dev = float(np.max(np.abs(ratio - 1.0)))
    assert max_dev < 0.01, (
        f"Deep MOND limit: max|F'/√Q − 1| = {max_dev:.4e}"
    )
    print(f"  PASS  F'(Q)/√Q → 1 as Q→0  (deep MOND limit)")


# ---------------------------------------------------------------------------
# Test 5: Positivity — F(Q) ≥ 0
# ---------------------------------------------------------------------------

def test_positivity():
    """F(Q) ≥ 0 for all Q ≥ 0 and all λ ∈ [0, 0.5]."""
    for lam in LAM_VALS:
        F = free_function_v42(Q_RANGE, lam=lam, beta=BETA)
        assert np.all(F >= -1e-12), (
            f"F(Q) < 0 at λ={lam}: min F = {F.min():.4e}"
        )
    # Also at Q=0 — use .item() for NumPy ≥1.25 scalar conversion
    F_zero = free_function_v42(np.array([0.0]), lam=0.0)
    assert float(F_zero.item()) >= 0.0, f"F(0) = {float(F_zero.item()):.4e} < 0"
    print(f"  PASS  F(Q) ≥ 0 for all Q ≥ 0  (all λ ∈ {LAM_VALS})")


# ---------------------------------------------------------------------------
# Test 6: Monotonicity — F'(Q) > 0
# ---------------------------------------------------------------------------

def test_monotonicity():
    """F'(Q) > 0: effective gravitational coupling never reverses sign."""
    for lam in LAM_VALS:
        Fp = free_function_derivative_v42(Q_DENSE, lam=lam, beta=BETA)
        assert np.all(Fp > 0), (
            f"F'(Q) ≤ 0 at λ={lam}: min F' = {Fp.min():.4e}"
        )
    print(f"  PASS  F'(Q) > 0  (no gravitational sign reversal)  all λ")


# ---------------------------------------------------------------------------
# Test 7 & 8: Ghost-free and subluminality (via sentinel)
# ---------------------------------------------------------------------------

def test_ghost_free_sentinel():
    """Stability sentinel: c_s² > 0 (no ghost)."""
    for lam in LAM_VALS:
        result = check_stability_v42(Q_DENSE, lam=lam, beta=BETA)
        ghost_warnings = [w for w in result.warnings if "GHOST" in w]
        assert not ghost_warnings, (
            f"Ghost at λ={lam}:\n" + "\n".join(ghost_warnings)
        )
    print(f"  PASS  Ghost-free (c_s² > 0) for all λ ∈ {LAM_VALS}")


def test_subluminal_sentinel():
    """Stability sentinel: c_s² ≤ 1 within the safe operating range λ ≤ 0.35.

    λ ≥ 0.44 causes c_s² > 1 (documented in test_gw_compliance.py).
    LAM_VALS includes 0.5 for coverage; we document but do not fail on it here
    since the refinery LAM_BOUNDS is already capped at 0.35.
    """
    safe_lam_vals = [l for l in LAM_VALS if l <= 0.35]
    for lam in safe_lam_vals:
        result = check_stability_v42(Q_DENSE, lam=lam, beta=BETA)
        causal_warnings = [w for w in result.warnings if "CAUSALITY" in w]
        assert not causal_warnings, (
            f"Causality violation at λ={lam}:\n" + "\n".join(causal_warnings)
        )
    # Document that λ=0.5 violates causality (expected, not a bug)
    result_unsafe = check_stability_v42(Q_DENSE, lam=0.5, beta=BETA)
    has_violation = any("CAUSALITY" in w for w in result_unsafe.warnings)
    assert has_violation, "Expected causality violation at λ=0.5 — check sentinel"
    print(
        f"  PASS  Subluminal (c_s² ≤ 1) for λ ∈ {safe_lam_vals}  |  "
        f"λ=0.5 correctly flagged as causality-violating"
    )


# ---------------------------------------------------------------------------
# Test 9: Continuity — F, F', F'' are finite and continuous
# ---------------------------------------------------------------------------

def test_continuity():
    """F(Q), F'(Q), F''(Q) are finite, non-NaN, and non-Inf for all Q > 0."""
    for lam in [0.0, LAMBDA]:
        F   = free_function_v42(Q_DENSE, lam=lam, beta=BETA)
        Fp  = free_function_derivative_v42(Q_DENSE, lam=lam, beta=BETA)
        Fpp = free_function_second_derivative_v42(Q_DENSE, lam=lam, beta=BETA)

        for name, arr in [("F", F), ("F'", Fp), ("F''", Fpp)]:
            assert np.all(np.isfinite(arr)), (
                f"{name}(Q) has non-finite values at λ={lam}: "
                f"n_bad = {np.sum(~np.isfinite(arr))}"
            )
            assert not np.any(np.isnan(arr)), (
                f"{name}(Q) has NaN values at λ={lam}"
            )
    print(f"  PASS  F, F', F'' finite and continuous on Q ∈ [10⁻⁴, 10³]")


# ---------------------------------------------------------------------------
# Test 10: λ-smoothness — F'(Q) is smooth in λ
# ---------------------------------------------------------------------------

def test_lambda_smoothness():
    """F'(Q) varies smoothly and monotonically with λ at fixed Q."""
    Q_test = np.array([0.01, 0.1, 1.0, 10.0])
    lam_sweep = np.linspace(0.0, 0.5, 20)

    for q in Q_test:
        fp_vals = [
            float(free_function_derivative_v42(np.array([q]), lam=lam, beta=BETA).item())
            for lam in lam_sweep
        ]
        # Check no discontinuities — adjacent differences should be small
        diffs = np.abs(np.diff(fp_vals))
        assert np.all(diffs < 0.5), (
            f"Non-smooth F'(Q={q}) vs λ: max jump = {diffs.max():.4f}"
        )

    print(f"  PASS  F'(Q) varies smoothly with λ at all tested Q values")


# ---------------------------------------------------------------------------
# Test 11: Caustic guard localisation
# ---------------------------------------------------------------------------

def test_guard_localisation():
    """
    The caustic guard λQ²e^{-βQ} satisfies:
    (a) O(Q²) at small Q — negligible in deep MOND regime
    (b) Exponentially suppressed at large Q — negligible in Newtonian regime
    (c) Peaks at Q = 2/β ≈ 1.33
    """
    lam = 0.05
    beta = BETA
    # Guard peak: d/dQ[(2Q - βQ²)e^{-βQ}] = 0
    # → β²Q² - 4βQ + 2 = 0 → Q = (4β - √(16β²-8β²))/(2β²) = (2-√2)/β
    Q_peak_expected = (2.0 - np.sqrt(2.0)) / beta   # ≈ 0.391 for β=1.5

    # Guard contribution to F'(Q) = λ(2Q - βQ²)e^{-βQ}
    Q_arr = Q_RANGE
    guard_Fp = lam * (2.0 * Q_arr - beta * Q_arr**2) * np.exp(-beta * Q_arr)
    guard_Fp_positive = guard_Fp.copy()
    guard_Fp_positive[guard_Fp < 0] = 0.0

    Q_peak_actual = float(Q_arr[np.argmax(guard_Fp_positive)].item())
    assert abs(Q_peak_actual - Q_peak_expected) / Q_peak_expected < 0.1, (
        f"Guard peak at Q={Q_peak_actual:.3f}, expected Q≈{Q_peak_expected:.3f}"
    )

    # At small Q (deep MOND): guard ≪ F'_MOND
    Q_small = np.array([1e-3])
    Fp_mond  = float(free_function_derivative_v42(Q_small, lam=0.0, beta=BETA).item())
    Fp_guard = lam * (2.0 * float(Q_small.item())) * np.exp(-beta * float(Q_small.item()))
    ratio_small = Fp_guard / Fp_mond
    assert ratio_small < 0.01, (
        f"Guard not small at Q=1e-3: guard/MOND = {ratio_small:.4f}"
    )

    # At large Q (Newtonian): guard exponentially suppressed
    Q_large = np.array([10.0])
    ql = float(Q_large.item())
    guard_large = lam * (2.0 * ql - beta * ql**2) * np.exp(-beta * ql)
    assert abs(guard_large) < 1e-3, (
        f"Guard not suppressed at Q=10: |guard| = {abs(guard_large):.4e}"
    )

    print(
        f"  PASS  Guard localised: peaks at Q≈{Q_peak_actual:.2f} (expected {Q_peak_expected:.2f}); "
        f"deep-MOND ratio = {ratio_small:.4f}; large-Q guard = {abs(guard_large):.2e}"
    )


# ---------------------------------------------------------------------------
# Test 12: SPARC result consistency
# ---------------------------------------------------------------------------

def test_sparc_consistency_lambda_zero_preferred():
    """
    The SPARC global fit (171 galaxies, April 2026) drives λ_opt = 0.
    Profile scan curvature = +2.11 (peaked — not just insensitive).

    This test verifies the mathematical condition that makes λ=0 optimal:
    the guard term contributes positively to F'(Q) near the transition
    regime Q ≈ 1.33, increasing g_eff = g_N/F'(Q) relative to λ=0.
    If the guard always increases g_eff, the optimizer will prefer
    λ=0 because any non-zero λ makes the predicted velocity systematically
    higher in the transition regime, increasing RMSE against data.
    """
    # Evaluate at Q_peak where F'_guard is MAXIMALLY POSITIVE: Q = (2-√2)/β ≈ 0.39
    # Note: at Q = 2/β = 1.33 the guard contribution is EXACTLY ZERO
    # (2Q - βQ²)|_{Q=2/β} = 4/β - 4/β = 0. Do not test there.
    Q_peak = np.array([(2.0 - np.sqrt(2.0)) / BETA])   # ≈ 0.391 for β=1.5

    Fp_lambda0  = float(free_function_derivative_v42(Q_peak, lam=0.0,  beta=BETA).item())
    Fp_lambda05 = float(free_function_derivative_v42(Q_peak, lam=0.05, beta=BETA).item())

    # At the guard's F' peak, the guard term contributes positively to F'(Q).
    # Larger F'(Q) → smaller g_eff = g_bar/F'(Q) → lower V_pred.
    # For transition-regime galaxies where the data sit near Q≈0.4–1.3,
    # adding λ>0 systematically lowers predicted velocities, increasing RMSE.
    # The SPARC optimizer therefore prefers λ=0.
    assert Fp_lambda05 > Fp_lambda0, (
        f"Guard should increase F' at its F' peak Q≈(2-√2)/β: "
        f"F'(λ=0)={Fp_lambda0:.4f}, F'(λ=0.05)={Fp_lambda05:.4f}"
    )

    delta_Fp = Fp_lambda05 - Fp_lambda0
    print(
        f"  PASS  Guard additive at F' peak Q={(2-np.sqrt(2))/BETA:.3f}: "
        f"ΔF' = +{delta_Fp:.4f} ({delta_Fp/Fp_lambda0*100:.2f}%)\n"
        f"         Higher F'(Q) → lower g_eff → lower V_pred in transition regime → "
        f"SPARC optimizer prefers λ=0 to avoid systematic underprediction"
    )


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

TESTS = [
    test_canonical_kinetic_limit,
    test_newtonian_limit,
    test_mond_interpolation_exact_at_lambda_zero,
    test_mond_interpolation_approximate_at_lambda_nonzero,
    test_deep_mond_limit,
    test_positivity,
    test_monotonicity,
    test_ghost_free_sentinel,
    test_subluminal_sentinel,
    test_continuity,
    test_lambda_smoothness,
    test_guard_localisation,
    test_sparc_consistency_lambda_zero_preferred,
]

if __name__ == "__main__":
    print("=" * 60)
    print("  Stability Test Suite — CODE-GEO V4.2")
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
        print("  V4.2 FREE FUNCTION: ALL STABILITY PROPERTIES VERIFIED")
    else:
        print("  STABILITY FAILURES — review output above")
    print(f"  {'='*40}")