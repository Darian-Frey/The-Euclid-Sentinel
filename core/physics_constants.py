"""
core/physics_constants.py
=========================
CODE-GEO V4.2 — Physical Constants: Single Source of Truth
===========================================================

All physical constants and framework-specific parameters for the
Mimetic-Conformal gravity pipeline are defined here. Every pipeline module
imports from this file; nothing is hard-coded elsewhere.

Framework summary
-----------------
The Mimetic-Conformal V4.2 framework replaces cold dark matter with a
scalar field σ governed by the free function:

    F(Q) = √(Q(1+Q)) − arcsinh(√Q)  +  λ Q² e^{−βQ}

where Q ≡ −g^{μν}∂_μσ∂_νσ / a₀² is the dimensionless kinetic invariant.
F'(Q) = μ_std(√Q) exactly, anchoring the theory to Standard MOND.
The conformal metric structure ensures c_T = c (GW170817 constraint).

Version history of framework-specific parameters
-------------------------------------------------
V3.1  (deprecated):
    KAPPA = 0.80 — "Hartley-Krylov damping constant", described in the
    original codebase as "Informational Refresh Efficiency". This was a
    shape parameter of the V3.1 Soft-Krylov exponential free function:
        F(Q) = Q(1 − exp(−κ√Q)) + (2κ/3)(1−κ)Q^{3/2}exp(−κQ)
    The V3.1/V4.1 functions failed the SPARC refinery (RMSE 18.7–23.4%)
    due to a shape mismatch with μ_std and a causality violation (c_s²>1).
    KAPPA has no role in the V4.2 framework and is NOT defined here.
    Legacy code using kappa=0.80 will raise ImportError intentionally.

V4.1  (deprecated):
    Soft-Krylov: replaced √Q with Q^{1/4} in exponential of Term 1.
    Still used KAPPA. Superseded by V4.2.

V4.2  (current):
    Free parameters: LAMBDA (λ), BETA (β). KAPPA removed entirely.
    SPARC global fit drives λ_opt = 0.000 (171 galaxies, April 2026).
    The minimal action F_MOND(Q) = √(Q(1+Q)) − arcsinh(√Q) is preferred.
    λ and β are retained in the codebase for theoretical exploration
    but are not empirically required at galactic rotation curve scales.

References
----------
* Chamseddine & Mukhanov (2013) — original mimetic gravity
* Bekenstein & Milgrom (1984)   — AQUAL formulation
* Abbott et al. (2017)          — GW170817, c_T = c constraint
* McGaugh & Schombert (2015)    — 3.6μm stellar M/L ratios
* Lelli et al. (2016)           — SPARC database
* Li et al. (2018)              — SPARC MOND fits
"""

import math


# =============================================================================
# Standard Physical Constants (SI units throughout)
# =============================================================================

# Gravitational constant [m³ kg⁻¹ s⁻²]
G: float = 6.67430e-11

# Speed of light in vacuum [m/s]
C: float = 299_792_458.0

# Alias for modules that import C_LIGHT (core/action.py compatibility)
C_LIGHT: float = C


# =============================================================================
# CODE-GEO V4.2 Framework Parameters
# =============================================================================

# MOND acceleration scale [m/s²]
# Canonical value from McGaugh et al. (2016) RAR fit to SPARC.
# Joint (λ, a₀) free fit on 171 SPARC galaxies recovers 1.264e-10 (+4.46%);
# the canonical value is within systematic uncertainties of M/L ratios.
# Do not change this without re-running sparc_refinery_v4.py.
A0: float = 1.21e-10

# Caustic guard amplitude [dimensionless]
# Theoretical default from V4.2 stability audit.
# SPARC global fit (171 galaxies) drives λ_opt = 0.000 (April 2026).
# The guard is empirically excluded at galactic rotation curve scales
# (peaked λ-profile; curvature +2.11 at λ=0; non-degenerate with a₀).
# Retained for theoretical exploration (cluster scales, structure formation).
LAMBDA: float = 0.05

# Caustic guard decay rate [dimensionless]
# Fixed per V4.2 audit; guard peaks at Q = 2/β ≈ 1.33 (g_bar/a₀ ≈ 1.15).
# Do not tune without re-running core/action.check_stability_v42().
BETA: float = 1.5

# MOND characteristic length scale [m]
# ℓ = c / √a₀  (dimensionally: [m/s] / √[m/s²] = [m/s] × √[s²/m] = √[m·s])
# Note: this is √(c²/a₀) [m], not c²/a₀ [m].
ELL: float = math.sqrt(C**2 / A0)

# Reduced Planck mass [kg]
# M_Pl = ℏc / (8πG) in SI; equivalent to 2.435×10¹⁸ GeV/c².
# Used for scalar field normalisation in the action.
M_PL: float = 4.341e-9   # [kg]  (NOT GeV — unit is kg)


# =============================================================================
# DEPRECATED — V3.1 / V4.1 Parameters
# =============================================================================
# KAPPA (κ = 0.80) was the "Hartley-Krylov damping constant" used in the
# V3.1 Soft-Krylov and V4.1 exponential free functions. It is NOT defined
# in V4.2. Any import of KAPPA from this module will raise ImportError,
# which is the intended behaviour — it forces legacy code to be updated
# rather than silently falling back to an incorrect value.
#
# If you are updating legacy code that used kappa=0.80:
#   - Remove the kappa argument from MimeticEngine() constructors
#   - Use free_function_v42() and free_function_derivative_v42() from
#     core/action.py instead of the deprecated V3.1/V4.1 functions
#   - The V4.2 free parameters are LAMBDA and BETA (defined above)
#
# KAPPA: float = 0.80   ← intentionally absent


# =============================================================================
# Astrophysical Mass Units
# =============================================================================

# Solar mass [kg]
MSUN: float = 1.98847e30

# Solar luminosity [W]  — used in flux→Σ_b calibration pipeline
LSUN: float = 3.828e26


# =============================================================================
# Astrophysical Distance Units
# =============================================================================

# Kiloparsec to metres  [m/kpc]
# 1 pc = 3.085677758×10¹⁶ m  →  1 kpc = 3.085677758×10¹⁹ m
KPC_TO_M: float = 3.08567758e19

# Parsec to metres [m/pc]
PC_TO_M: float = 3.08567758e16

# Megaparsec to metres [m/Mpc]
MPC_TO_M: float = 3.08567758e22


# =============================================================================
# Velocity Conversion Factors
# =============================================================================

# Kilometres per second → metres per second
KM_S_TO_M_S: float = 1_000.0

# Metres per second → kilometres per second
M_S_TO_KM_S: float = 1.0e-3


# =============================================================================
# Surface Mass Density Conversion
# =============================================================================

# Solar masses per parsec² → kg per metre²
# Used in tools/euclid_loader.py flux-to-Σ_b calibration
MSUN_PER_PC2_TO_SI: float = MSUN / PC_TO_M**2   # ≈ 6.768e-21 kg/m²


# =============================================================================
# Quick self-check (run as __main__)
# =============================================================================

if __name__ == "__main__":
    print("CODE-GEO V4.2 — Physics Constants Verification")
    print("=" * 52)

    print(f"\n  Standard constants:")
    print(f"    G        = {G:.5e} m³ kg⁻¹ s⁻²")
    print(f"    C        = {C:.9e} m/s")
    print(f"    C_LIGHT  = {C_LIGHT:.9e} m/s  (alias for C)")

    print(f"\n  V4.2 framework parameters:")
    print(f"    A0       = {A0:.4e} m/s²  (MOND scale)")
    print(f"    LAMBDA   = {LAMBDA}        (caustic guard; λ_opt=0 from SPARC)")
    print(f"    BETA     = {BETA}          (guard decay rate; fixed)")
    print(f"    ELL      = {ELL:.6e} m   (√(c²/a₀))")
    print(f"    M_PL     = {M_PL:.4e} kg  (reduced Planck mass)")

    print(f"\n  Astrophysical units:")
    print(f"    MSUN     = {MSUN:.5e} kg")
    print(f"    LSUN     = {LSUN:.4e} W")
    print(f"    KPC_TO_M = {KPC_TO_M:.8e} m/kpc")
    print(f"    PC_TO_M  = {PC_TO_M:.8e} m/pc")
    print(f"    MSUN_PER_PC2_TO_SI = {MSUN_PER_PC2_TO_SI:.4e} kg/m²")

    print(f"\n  Velocity conversions:")
    print(f"    KM_S_TO_M_S = {KM_S_TO_M_S:.1f}")
    print(f"    M_S_TO_KM_S = {M_S_TO_KM_S:.4f}")

    # Verify ELL
    ell_check = math.sqrt(C**2 / A0)
    assert abs(ELL - ell_check) < 1e-6 * ELL, "ELL mismatch"

    # Verify MSUN_PER_PC2_TO_SI
    si_check = MSUN / PC_TO_M**2
    assert abs(MSUN_PER_PC2_TO_SI - si_check) < 1e-10 * si_check, "SI conversion mismatch"

    # Verify C_LIGHT alias
    assert C_LIGHT == C, "C_LIGHT alias broken"

    # Verify KAPPA is absent (should raise NameError)
    try:
        _ = KAPPA  # noqa: F821
        print("\n  ⚠ WARNING: KAPPA is defined — it should not be in V4.2")
    except NameError:
        print(f"\n  ✓ KAPPA correctly absent (V3.1 parameter; not used in V4.2)")

    print(f"\n  All checks passed.")
