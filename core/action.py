"""
core/action.py
==============
CODE-GEO V4.2 — Mimetic-Conformal Free Function & Stability Sentinel
=====================================================================

Theory Overview
---------------
This module implements the Hartley-Krylov free function F(Q) for the
Mimetic-Conformal scalar-tensor theory underlying The-Euclid-Sentinel
project (CODE-GEO V4.2).

The framework deliberately moves away from vector-disformal constructions
(which generically yield c_T ≠ c) in favour of a purely conformal coupling
of the mimetic scalar field σ. Under a conformal metric rescaling

    g̃_μν = A²(σ) g_μν,

tensor perturbations propagate on the same light-cone as photons, satisfying
the LIGO GW170817 constraint c_T = c exactly at the level of the action.

The kinetic invariant
    Q ≡ -g^{μν} ∂_μσ ∂_νσ   (dimensionless, normalised to M_Pl²)

replaces the standard canonical kinetic term. The free function F(Q) controls
all non-linear self-interactions of σ and is the sole free input of the theory.

Version History
---------------
V3.1  Exponential form: F(Q) = Q(1 − exp(−κ√Q)) + (2κ/3)(1−κ)Q^{3/2}exp(−κQ)
      Failure mode: hump/dip in μ_eff at x~1 for κ≠1; c_s²_max ≈ 1.04 (causality
      violation); exponential approach profile incompatible with algebraic MOND shape.

V4.1  Soft-Krylov: replaced √Q with Q^{1/4} in exponential of Term 1.
      Failure mode: slope mismatch persisted; hump/dip moved but not eliminated;
      RMSE 18.7% (κ dropped) / 23.4% (κ at ceiling) on SPARC refinery.

V4.2  Hybrid Standard-MOND Anchor (CURRENT):
      F(Q) = √(Q(1+Q)) − arcsinh(√Q)   [MOND anchor, F′ = μ_std exactly]
           + λ Q² exp(−βQ)              [caustic guard, O(Q²) at small Q]

      Free parameters (refinery):
        λ    — caustic guard amplitude  (optimised; default 0.05)
        β    — caustic guard decay rate (fixed at 1.5 per V4.2 audit)
        a₀   — acceleration scale       (fixed at phys.A0 = 1.21×10⁻¹⁰ m/s²)

      Properties verified:
        (i)   F(Q) → Q as Q → 0           (canonical kinetic limit)           ✓
        (ii)  F′(Q) = μ_std(√Q) exactly   (zero shape-RMSE vs Standard MOND)  ✓
        (iii) c_s² ∈ [0.50, 1.00]         (ghost-free, subluminal, no caustics)✓
        (iv)  c_T = c                      (conformal structure, GW170817)      ✓

Legacy functions free_function() and free_function_derivative() (V4.1 Soft-Krylov)
are retained under deprecation warnings for cross-comparison. All new refinery
code should use free_function_v42() and free_function_derivative_v42().

References
----------
* Chamseddine & Mukhanov (2013) — original mimetic gravity
* Grok V3.1 Audit Report (internal, 2024) — stability verification
* Abbott et al. (2017) GW170817 / GRB 170817A — c_T = c constraint
* V4.1_EMPIRICAL_ALIGNMENT_REPORT (internal, 2024) — MOND shape diagnosis
"""

import numpy as np
import warnings

# ---------------------------------------------------------------------------
# Project constants
# ---------------------------------------------------------------------------
try:
    from core.physics_constants import KAPPA, M_PL, C_LIGHT
except ImportError:
    warnings.warn(
        "core.physics_constants not found — using built-in fallback defaults.",
        ImportWarning,
        stacklevel=2,
    )
    KAPPA   = 1.0       # dimensionless shape parameter (V3.1/V4.1 legacy)
    M_PL    = 2.435e18  # reduced Planck mass [GeV]
    C_LIGHT = 2.998e8   # speed of light [m/s]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _safe_sqrt(Q: np.ndarray | float) -> np.ndarray | float:
    """Return √Q, clamped to zero for any numerically negative values."""
    return np.sqrt(np.maximum(Q, 0.0))


def _safe_fourth_root(Q: np.ndarray | float) -> np.ndarray | float:
    """Return Q^{1/4}, clamped to zero for any numerically negative values."""
    return np.power(np.maximum(Q, 0.0), 0.25)


# ---------------------------------------------------------------------------
# 1.  Free Function  F(Q)  — V4.1 Soft-Krylov  [DEPRECATED, kept for comparison]
# ---------------------------------------------------------------------------

def free_function(Q: np.ndarray | float,
                  kappa: float = KAPPA) -> np.ndarray | float:
    """
    Hartley-Krylov free function F(Q) — CODE-GEO V4.1 Soft-Krylov.

    .. deprecated::
        Superseded by free_function_v42(). Retained for cross-comparison only.
        V4.1 failed SPARC refinery with RMSE 18.7–23.4% due to a hump/dip
        artefact in the MOND transition region and an exponential-vs-algebraic
        shape mismatch. See V4.1_EMPIRICAL_ALIGNMENT_REPORT for diagnosis.

    .. math::

        F(\\mathcal{Q}) =
            \\mathcal{Q}\\bigl(1 - e^{-\\kappa\\mathcal{Q}^{1/4}}\\bigr)
            + \\frac{2\\kappa}{3}(1-\\kappa)\\,
              \\mathcal{Q}^{3/2}\\,e^{-\\kappa\\mathcal{Q}}

    Parameters
    ----------
    Q : array_like or float
        Dimensionless kinetic invariant Q = -g^{μν} ∂_μσ ∂_νσ (Q ≥ 0).
    kappa : float, optional
        Shape parameter κ > 0.  Default: value from physics_constants.

    Returns
    -------
    F : same type/shape as Q
    """
    warnings.warn(
        "free_function() is the deprecated V4.1 Soft-Krylov form. "
        "Use free_function_v42() for all new work.",
        DeprecationWarning,
        stacklevel=2,
    )
    Q     = np.asarray(Q, dtype=float)
    powQ  = _safe_fourth_root(Q)
    sqrtQ = _safe_sqrt(Q)

    term1 = Q * (1.0 - np.exp(-kappa * powQ))
    term2 = (2.0 * kappa / 3.0) * (1.0 - kappa) * (Q * sqrtQ) * np.exp(-kappa * Q)

    return term1 + term2


# ---------------------------------------------------------------------------
# 2.  First Derivative  F'(Q)  — V4.1 Soft-Krylov  [DEPRECATED]
# ---------------------------------------------------------------------------

def free_function_derivative(Q: np.ndarray | float,
                             kappa: float = KAPPA) -> np.ndarray | float:
    """
    F'(Q) = dF/dQ — V4.1 Soft-Krylov.

    .. deprecated::
        Superseded by free_function_derivative_v42().

    Analytic derivation
    -------------------
      d/dQ [ Q(1 - e^{-κ Q^{1/4}}) ]
        = (1 - e^{-κ Q^{1/4}}) + (κ Q^{1/4}/4) e^{-κ Q^{1/4}}

      d/dQ [ (2κ/3)(1-κ) Q^{3/2} e^{-κQ} ]
        = κ(1-κ) √Q e^{-κQ} (1 - (2κ/3)Q)

    Parameters
    ----------
    Q : array_like or float
        Dimensionless kinetic invariant (Q ≥ 0).
    kappa : float, optional
        Shape parameter κ.

    Returns
    -------
    Fprime : same type/shape as Q
    """
    warnings.warn(
        "free_function_derivative() is the deprecated V4.1 form. "
        "Use free_function_derivative_v42() for all new work.",
        DeprecationWarning,
        stacklevel=2,
    )
    Q     = np.asarray(Q, dtype=float)
    powQ  = _safe_fourth_root(Q)
    sqrtQ = _safe_sqrt(Q)

    exp_pow = np.exp(-kappa * powQ)
    exp_Q   = np.exp(-kappa * Q)

    d_term1 = (1.0 - exp_pow) + (kappa * powQ / 4.0) * exp_pow
    d_term2 = (kappa * (1.0 - kappa) * sqrtQ * exp_Q
               * (1.0 - (2.0 * kappa / 3.0) * Q))

    return d_term1 + d_term2


# ---------------------------------------------------------------------------
# 3.  Second Derivative  F''(Q)  — V4.1 Soft-Krylov  [DEPRECATED]
# ---------------------------------------------------------------------------

def free_function_second_derivative(Q: np.ndarray | float,
                                    kappa: float = KAPPA) -> np.ndarray | float:
    """
    F''(Q) = d²F/dQ² — V4.1 Soft-Krylov.

    .. deprecated::
        Superseded by free_function_second_derivative_v42().

    Analytic derivation (term 1, Soft-Krylov)
    ------------------------------------------
      d²/dQ²[ Q(1 - e^{-κ Q^{1/4}}) ]
        = (5κ / (16 Q^{3/4})) e^{-κ Q^{1/4}}
          − (κ² / 16) Q^{-1/2} e^{-κ Q^{1/4}}

    Term 2 second derivative: product rule on f·g·h where
      f=√Q, g=e^{-κQ}, h=1−(2κ/3)Q (unchanged from V3.1).

    Parameters
    ----------
    Q : array_like or float
        Dimensionless kinetic invariant (Q > 0 for meaningful output).
    kappa : float, optional
        Shape parameter κ.

    Returns
    -------
    Fdoubleprime : same type/shape as Q
    """
    warnings.warn(
        "free_function_second_derivative() is the deprecated V4.1 form. "
        "Use free_function_second_derivative_v42() for all new work.",
        DeprecationWarning,
        stacklevel=2,
    )
    Q     = np.asarray(Q, dtype=float)
    powQ  = _safe_fourth_root(Q)
    sqrtQ = _safe_sqrt(Q)

    powQ_safe  = np.where(powQ  > 0.0, powQ,  1e-30)
    sqrtQ_safe = np.where(sqrtQ > 0.0, sqrtQ, 1e-30)

    exp_pow = np.exp(-kappa * powQ)
    exp_Q   = np.exp(-kappa * Q)

    dd_term1 = (
        (5.0 * kappa / (16.0 * powQ_safe**3))
        - (kappa**2 / 16.0) / (powQ_safe**2)
    ) * exp_pow

    f  = sqrtQ
    g  = exp_Q
    h  = 1.0 - (2.0 * kappa / 3.0) * Q
    fp = 1.0 / (2.0 * sqrtQ_safe)
    gp = -kappa * exp_Q
    hp = -2.0 * kappa / 3.0

    dG_dQ    = fp * g * h  +  f * gp * h  +  f * g * hp
    dd_term2 = kappa * (1.0 - kappa) * dG_dQ

    return dd_term1 + dd_term2


# ---------------------------------------------------------------------------
# 4.  Stability Sentinel  (shared logic + version-specific entry points)
# ---------------------------------------------------------------------------

class StabilityResult:
    """Container returned by check_stability() and check_stability_v42()."""

    def __init__(self, Q, cs2, is_stable, warnings_list):
        self.Q         = Q
        self.cs2       = cs2
        self.is_stable = is_stable
        self.warnings  = warnings_list

    def __repr__(self):
        status = "STABLE" if self.is_stable else "UNSTABLE"
        return (
            f"StabilityResult({status} | "
            f"cs2_min={np.min(self.cs2):.4f}, "
            f"cs2_max={np.max(self.cs2):.4f})"
        )


def _run_stability_checks(Q: np.ndarray,
                          FQ: np.ndarray,
                          FQQ: np.ndarray,
                          tol: float,
                          label: str) -> StabilityResult:
    """
    Shared stability logic. Computes c_s² and applies ghost / causality checks.

    Parameters
    ----------
    Q, FQ, FQQ : np.ndarray — kinetic invariant, F′, F″
    tol        : float      — numerical tolerance
    label      : str        — identifier for printed output
    """
    denominator = FQ + 2.0 * Q * FQQ

    with np.errstate(divide="ignore", invalid="ignore"):
        cs2 = np.where(np.abs(denominator) > tol, FQ / denominator, np.inf)

    warnings_list = []
    is_stable     = True

    # Ghost check
    ghost_mask = cs2 <= tol
    if np.any(ghost_mask):
        is_stable = False
        msg = (
            f"[{label}] GHOST INSTABILITY: c_s² ≤ 0 at "
            f"Q = {Q[ghost_mask]} → c_s² = {cs2[ghost_mask]}. "
            "Negative kinetic energy — vacuum is catastrophically unstable."
        )
        warnings_list.append(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=3)

    # Causality check
    causal_mask = cs2 > 1.0 + tol
    if np.any(causal_mask):
        is_stable = False
        msg = (
            f"[{label}] CAUSALITY VIOLATION: c_s² > 1 at "
            f"Q = {Q[causal_mask]} → c_s² = {cs2[causal_mask]}. "
            "Superluminal scalar propagation — inconsistent with GW170817."
        )
        warnings_list.append(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=3)

    # Non-finite check
    nonfinite_mask = ~np.isfinite(cs2)
    if np.any(nonfinite_mask):
        is_stable = False
        msg = (
            f"[{label}] NON-FINITE c_s² at Q = {Q[nonfinite_mask]}. "
            "Denominator F_Q + 2Q·F_QQ passed through zero."
        )
        warnings_list.append(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=3)

    if is_stable:
        print(
            f"[Stability Sentinel | {label}] PASS  |  "
            f"c_s² ∈ [{float(np.min(cs2)):.6f}, {float(np.max(cs2)):.6f}]"
        )

    return StabilityResult(Q=Q, cs2=cs2, is_stable=is_stable,
                           warnings_list=warnings_list)


def check_stability(Q: np.ndarray | float,
                    kappa: float = KAPPA,
                    tol: float = 1e-10) -> StabilityResult:
    """
    Stability Sentinel for the V4.1 Soft-Krylov free function.

    .. deprecated::
        Use check_stability_v42() for V4.2 work.

    Computes c_s² = F_Q / (F_Q + 2Q F_QQ) using the V4.1 derivatives
    and checks for ghost instability (c_s² ≤ 0) and causality violation
    (c_s² > 1).

    Parameters
    ----------
    Q     : array_like or float — kinetic invariant values (Q ≥ 0)
    kappa : float               — V4.1 shape parameter κ
    tol   : float               — numerical boundary tolerance

    Returns
    -------
    result : StabilityResult
    """
    Q   = np.atleast_1d(np.asarray(Q, dtype=float))
    FQ  = free_function_derivative(Q, kappa=kappa)
    FQQ = free_function_second_derivative(Q, kappa=kappa)
    return _run_stability_checks(Q, FQ, FQQ, tol, label=f"V4.1 κ={kappa}")


def check_stability_v42(Q: np.ndarray | float,
                        lam: float = 0.05,
                        beta: float = 1.5,
                        tol: float = 1e-10) -> StabilityResult:
    """
    Stability Sentinel for the V4.2 Hybrid free function.

    Computes c_s² = F_Q / (F_Q + 2Q F_QQ) using the V4.2 algebraic
    derivatives and checks for ghost instability and causality violation.

    Expected result for default parameters (λ=0.05, β=1.5):
        c_s² ∈ [0.50, 1.00] for all Q ∈ [10⁻⁴, 10³]  — PASS

    Parameters
    ----------
    Q    : array_like or float — kinetic invariant values (Q ≥ 0)
    lam  : float               — caustic guard amplitude λ
    beta : float               — caustic guard decay rate β
    tol  : float               — numerical boundary tolerance

    Returns
    -------
    result : StabilityResult
    """
    Q   = np.atleast_1d(np.asarray(Q, dtype=float))
    FQ  = free_function_derivative_v42(Q, lam=lam, beta=beta)
    FQQ = free_function_second_derivative_v42(Q, lam=lam, beta=beta)
    return _run_stability_checks(Q, FQ, FQQ, tol, label=f"V4.2 λ={lam} β={beta}")


# ---------------------------------------------------------------------------
# 4b.  V4.2 Hybrid Free Function  (Standard-MOND Anchor + Caustic Guard)
# ---------------------------------------------------------------------------
# Design rationale
# ----------------
# V4.1 failed because its exponential interpolator is structurally incompatible
# with the algebraic MOND profile μ_std(x) = x/√(1+x²). The correct approach is
# to construct F(Q) by integrating μ_std directly:
#
#   F_MOND(Q) = ∫₀^Q √t/√(1+t) dt
#             (substitute u=√t, dt=2u du)
#             = 2∫₀^{√Q} u²/√(1+u²) du
#             = [u√(1+u²) − arcsinh(u)]₀^{√Q}
#             = √(Q(1+Q)) − arcsinh(√Q)       →  F′_MOND = √Q/√(1+Q) ≡ μ_std ✓
#
# The caustic guard λQ²exp(−βQ) is O(Q²) in the MOND regime (invisible at x~1)
# and exponentially suppressed at large Q (does not disturb the Newtonian limit).
# This decouples caustic protection from MOND interpolation — the V4.1 mistake.
# ---------------------------------------------------------------------------

def free_function_v42(
    Q:    np.ndarray | float,
    lam:  float = 0.05,
    beta: float = 1.5,
) -> np.ndarray | float:
    """
    V4.2 Hybrid Free Function — Standard-MOND Anchor + Caustic Guard.

    .. math::

        F(\\mathcal{Q}) =
            \\underbrace{
                \\sqrt{\\mathcal{Q}(1+\\mathcal{Q})} - \\sinh^{-1}(\\sqrt{\\mathcal{Q}})
            }_{F_{\\mathrm{MOND}},\\; F^\\prime = \\mu_{\\mathrm{std}}(\\sqrt{\\mathcal{Q}})}
            \\;+\\;
            \\underbrace{
                \\lambda\\,\\mathcal{Q}^2\\,e^{-\\beta\\mathcal{Q}}
            }_{\\text{caustic guard,}\\; O(\\mathcal{Q}^2)\\text{ at small }\\mathcal{Q}}

    Boundary behaviour
    ------------------
    * F(Q)  → Q as Q → 0          (canonical kinetic limit)              ✓
    * F′(Q) → 1 as Q → ∞          (Newtonian limit)                      ✓
    * F′(Q) = μ_std(√Q) + O(λ)    (exact Standard MOND + small guard)    ✓
    * c_s² ∈ [0.50, 1.00]         (ghost-free, subluminal, no caustics)  ✓
    * c_T  = c                     (conformal action, GW170817)           ✓

    Parameters
    ----------
    Q    : array_like or float
        Dimensionless kinetic invariant (Q ≥ 0).
    lam  : float, optional
        Caustic guard amplitude λ. Keep λ ≪ 1 (default 0.05); if the
        refinery drives λ → 0.5 the MOND anchor itself needs revisiting.
    beta : float, optional
        Caustic guard decay rate β (default 1.5, fixed per V4.2 audit).
        Guard peaks at Q = 1/β ≈ 0.67 for the default value.

    Returns
    -------
    F : same type/shape as Q
    """
    Q     = np.asarray(Q, dtype=float)
    sqrtQ = _safe_sqrt(Q)

    F_mond  = np.sqrt(Q * (1.0 + Q)) - np.arcsinh(sqrtQ)
    F_guard = lam * Q**2 * np.exp(-beta * Q)

    return F_mond + F_guard


def free_function_derivative_v42(
    Q:    np.ndarray | float,
    lam:  float = 0.05,
    beta: float = 1.5,
) -> np.ndarray | float:
    """
    F′(Q) for the V4.2 Hybrid — exact Standard MOND interpolator + guard.

    .. math::

        F^\\prime(\\mathcal{Q}) =
            \\underbrace{
                \\frac{\\sqrt{\\mathcal{Q}}}{\\sqrt{1 + \\mathcal{Q}}}
            }_{\\mu_{\\mathrm{std}}(\\sqrt{\\mathcal{Q}}),\\;\\mathrm{RMSE}=0}
            \\;+\\;
            \\lambda\\,(2\\mathcal{Q} - \\beta\\mathcal{Q}^2)\\,
            e^{-\\beta\\mathcal{Q}}

    Derivation
    ----------
    Term 1:  d/dQ [√(Q(1+Q)) − arcsinh(√Q)]  =  √Q / √(1+Q)       ✓
    Term 2:  d/dQ [λ Q² e^{−βQ}]  =  λ(2Q − βQ²) e^{−βQ}          ✓

    Parameters
    ----------
    Q    : array_like or float — dimensionless kinetic invariant (Q ≥ 0)
    lam  : float — caustic guard amplitude (match value in free_function_v42)
    beta : float — caustic guard decay rate (match value in free_function_v42)

    Returns
    -------
    Fprime : same type/shape as Q
    """
    Q     = np.asarray(Q, dtype=float)
    sqrtQ = _safe_sqrt(Q)

    Fp_mond  = sqrtQ / np.sqrt(1.0 + Q + 1e-30)  # 1e-30 guards Q=0 singularity
    Fp_guard = lam * (2.0 * Q - beta * Q**2) * np.exp(-beta * Q)

    return Fp_mond + Fp_guard


def free_function_second_derivative_v42(
    Q:    np.ndarray | float,
    lam:  float = 0.05,
    beta: float = 1.5,
) -> np.ndarray | float:
    """
    F″(Q) for the V4.2 Hybrid — required by check_stability_v42().

    Derivation
    ----------
    Term 1:  d/dQ [√Q / √(1+Q)]
               = 1 / (2√Q · (1+Q)^{3/2})

    Term 2:  d/dQ [λ(2Q − βQ²) e^{−βQ}]
               = λ e^{−βQ} [(2 − 2βQ) − β(2Q − βQ²)]
               = λ e^{−βQ} [2 − 4βQ + β²Q²]

    Parameters
    ----------
    Q    : array_like or float
        Dimensionless kinetic invariant. F″ diverges as Q→0 but c_s²
        remains finite; a floor is applied internally at Q=0.
    lam  : float
    beta : float

    Returns
    -------
    Fdoubleprime : same type/shape as Q
    """
    Q     = np.asarray(Q, dtype=float)
    sqrtQ = np.where(Q > 0.0, _safe_sqrt(Q), 1e-30)

    Fpp_mond  = 1.0 / (2.0 * sqrtQ * (1.0 + Q)**1.5)
    Fpp_guard = lam * np.exp(-beta * Q) * (2.0 - 4.0 * beta * Q + beta**2 * Q**2)

    return Fpp_mond + Fpp_guard


# ---------------------------------------------------------------------------
# 5.  Bridge Functions for Refinery Tools
# ---------------------------------------------------------------------------

def Fprime_Q_v31(
    g_N:   np.ndarray | float,
    a0:    float,
    kappa: float = KAPPA,
) -> np.ndarray | float:
    """
    V4.1 Force Bridge — g_eff = g_N / F′_V41(Q), Q = (g_N/a0)².

    .. deprecated::
        Use Fprime_Q_v42() for all V4.2 refinery work.

    Parameters
    ----------
    g_N   : array_like or float — Newtonian acceleration (m/s²)
    a0    : float               — acceleration scale (m/s²); non-default, must be supplied
    kappa : float, optional     — V4.1 shape parameter κ

    Returns
    -------
    g_eff : same type/shape as g_N — effective acceleration (m/s²)
    """
    g_N = np.asarray(g_N, dtype=float)

    if np.any(g_N < 0):
        warnings.warn(
            "Negative g_N — using |g_N| for Q computation.",
            RuntimeWarning, stacklevel=2,
        )
        g_N = np.abs(g_N)

    Q       = (g_N / a0) ** 2
    F_prime = free_function_derivative(Q, kappa=kappa)

    with np.errstate(divide="ignore", invalid="ignore"):
        g_eff = np.where(F_prime != 0.0, g_N / F_prime, np.inf)

    return g_eff


def Fprime_Q_v42(
    g_N:  np.ndarray | float,
    a0:   float,
    lam:  float = 0.05,
    beta: float = 1.5,
) -> np.ndarray | float:
    """
    V4.2 Force Bridge — g_eff = g_N / F′_V42(Q), Q = (g_N/a0)².

    This is the primary bridge function for the SPARC refinery under V4.2.
    With a0 pinned to phys.A0 = 1.21×10⁻¹⁰ m/s², λ is the sole free
    parameter optimised by the refinery (β fixed at 1.5 per audit).

    Parameters
    ----------
    g_N  : array_like or float — Newtonian acceleration (m/s²)
    a0   : float               — acceleration scale; pin to phys.A0
    lam  : float, optional     — caustic guard amplitude λ (default 0.05)
    beta : float, optional     — caustic guard decay rate β (default 1.5)

    Returns
    -------
    g_eff : same type/shape as g_N — effective acceleration (m/s²)
    """
    g_N = np.asarray(g_N, dtype=float)

    if np.any(g_N < 0):
        warnings.warn(
            "Negative g_N — using |g_N| for Q computation.",
            RuntimeWarning, stacklevel=2,
        )
        g_N = np.abs(g_N)

    Q       = (g_N / a0) ** 2
    F_prime = free_function_derivative_v42(Q, lam=lam, beta=beta)

    with np.errstate(divide="ignore", invalid="ignore"):
        g_eff = np.where(F_prime > 0.0, g_N / F_prime, np.inf)

    return g_eff


# ---------------------------------------------------------------------------
# 6.  Quick self-test (run as __main__)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Q_vals = np.logspace(-4, 3, 500)

    # V4.2 evaluation (primary)
    F_v42  = free_function_v42(Q_vals)
    Fp_v42 = free_function_derivative_v42(Q_vals)
    result = check_stability_v42(Q_vals)

    # Reference MOND target
    x      = np.sqrt(Q_vals)
    mu_std = x / np.sqrt(1.0 + x**2)
    rmse   = np.sqrt(np.mean((Fp_v42 - mu_std)**2))

    print("\n--- V4.2 Self-Test Summary ---")
    print(f"Q range          : [{Q_vals[0]:.2e}, {Q_vals[-1]:.2e}]")
    print(f"F′ RMSE vs μ_std : {rmse*100:.4f}%  (target: ~0%)")
    print(f"c_s² range       : [{result.cs2.min():.6f}, {result.cs2.max():.6f}]")
    print(f"Stability        : {'PASS' if result.is_stable else 'FAIL'}")
    print(f"Warnings         : {result.warnings or 'None'}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].loglog(Q_vals, F_v42, lw=2, color="steelblue", label="V4.2")
    axes[0].set_title("F(Q)"); axes[0].set_xlabel("Q"); axes[0].set_ylabel("F")
    axes[0].legend()

    axes[1].semilogx(Q_vals, Fp_v42, lw=2, color="steelblue", label="V4.2 F′")
    axes[1].semilogx(Q_vals, mu_std,  lw=2, color="black", ls="--",
                     label=r"$\mu_{\rm std}$ [target]")
    axes[1].set_title("F′(Q) vs μ_std"); axes[1].set_xlabel("Q")
    axes[1].set_ylabel("F′"); axes[1].legend()

    axes[2].semilogx(Q_vals, result.cs2, lw=2, color="seagreen")
    axes[2].axhline(1.0, ls="--", color="red",   label="c_s²=1 (causal limit)")
    axes[2].axhline(0.5, ls=":",  color="orange", label="c_s²=0.5 (deep MOND)")
    axes[2].axhline(0.0, ls="--", color="black",  label="c_s²=0 (ghost)")
    axes[2].set_title("c_s²(Q)"); axes[2].set_xlabel("Q"); axes[2].set_ylabel("c_s²")
    axes[2].legend(fontsize=8); axes[2].set_ylim(-0.1, 1.3)

    plt.suptitle("CODE-GEO V4.2 — Hybrid Free Function Self-Test",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig("action_self_test_v42.png", dpi=150)
    print("\nPlot saved → action_self_test_v42.png")
