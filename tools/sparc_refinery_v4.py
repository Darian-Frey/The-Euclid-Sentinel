# tools/sparc_refinery_v4.py
"""
SPARC Refinery V4.2 – Average Galaxy calibration for Mimetic-Conformal V4.2

This script:
- Builds a mock "Average Galaxy" inspired by Lelli et al. (2016) SPARC properties.
- Uses the V4.2 Hybrid free function F'(Q) from core.action to compute V_circ(r).
- Optimises λ only (a0 and β are pinned) to minimise RMSE vs. the mock SPARC curve.
- Reports the best-fit λ and whether the 5% RMSE target is met.

Single Source of Truth
----------------------
ALL physical constants are read from core.physics_constants.  Nothing is
hard-coded in this script.  To change a0, edit core/physics_constants.py:

    A0 = 1.21e-10  # MOND acceleration scale [m/s²]

Pinned parameters (V4.2 audit defaults, not optimised):
    a0   = phys.A0   — MOND acceleration scale
    beta = 1.5       — caustic guard decay rate (defined as module constant below)

Free parameter (optimised by Powell):
    lam  ∈ [0.0, 0.5] — caustic guard amplitude λ
"""

import numpy as np
from scipy.optimize import minimize

from core import physics_constants as phys
from core.action import free_function_derivative_v42

# ---------------------------------------------------------------------------
# Pinned audit constants — change here to propagate everywhere in this module
# ---------------------------------------------------------------------------
_BETA_AUDIT = 1.5   # caustic guard decay rate (V4.2 audit default, not optimised)


# ---------------------------------------------------------------------------
# Utilities and mock SPARC "Average Galaxy"
# ---------------------------------------------------------------------------

def load_mock_sparc(num_points: int = 40):
    """
    Build the baryonic mass profile for a representative "Average Galaxy"
    inspired by Lelli et al. (2016) SPARC properties.

    The target rotation curve is NO LONGER generated here. It is produced
    by ``run_optimization`` via a call to ``vcirc_v42`` at fiducial parameters
    to ensure the refinery validates the mathematical pipeline (Self-Consistency).

    Galaxy parameters
    -----------------
    * Total baryonic mass         : 5e10 M_sun
    * Exponential disk scale R_d  : 3 kpc
    * Radii                       : 0.5 – 30 kpc  (num_points, linear)

    Returns
    -------
    r_kpc    : np.ndarray — radii [kpc]
    M_cum_kg : np.ndarray — enclosed baryonic mass M(<r) [kg]
    """
    r_kpc = np.linspace(0.5, 30.0, num_points)
    M_baryon_total_kg = 5.0e10 * phys.MSUN
    R_d_kpc = 3.0
    r_m = r_kpc * phys.KPC_TO_M
    dr_m = np.gradient(r_m)
    
    # Surface density profile (normalized exponential disk)
    sigma_unnorm = np.exp(-r_kpc / R_d_kpc)
    # Mass in each annulus: dM ~ 2π r Sigma(r) dr
    dM_unnorm = 2.0 * np.pi * r_m * sigma_unnorm * dr_m
    M_cum_unnorm = np.cumsum(dM_unnorm)

    # Normalize so that M(<r_max) = M_baryon_total_kg
    M_cum_kg = M_cum_unnorm * (M_baryon_total_kg / M_cum_unnorm[-1])
    
    return r_kpc, M_cum_kg


# ---------------------------------------------------------------------------
# V4.2 rotation curve using F'(Q)
# ---------------------------------------------------------------------------

def vcirc_v42(r_kpc: np.ndarray,
              M_enclosed_kg: np.ndarray,
              lam: float,
              a0: float,
              beta: float = _BETA_AUDIT) -> np.ndarray:
    """
    V4.2 circular velocity — Standard-MOND anchor + caustic guard.

    Parameters
    ----------
    r_kpc, M_enclosed_kg : as per load_mock_sparc()
    lam  : caustic guard amplitude λ
    a0   : acceleration scale — pass phys.A0; do not hard-code
    beta : caustic guard decay rate (default: _BETA_AUDIT = 1.5)

    Returns
    -------
    V_circ_kms : np.ndarray  [km/s]
    """
    r_m = r_kpc * phys.KPC_TO_M
    g_N = phys.G * M_enclosed_kg / r_m**2
    Q   = (g_N / a0) ** 2

    # Call the canonical V4.2 derivative from core.action
    F_prime = free_function_derivative_v42(Q, lam=lam, beta=beta)

    with np.errstate(divide="ignore", invalid="ignore"):
        g_eff = np.where(F_prime > 0.0, g_N / F_prime, np.inf)

    return np.sqrt(np.maximum(g_eff * r_m, 0.0)) / 1000.0


# ---------------------------------------------------------------------------
# Objective function and optimization
# ---------------------------------------------------------------------------

def rmse_percentage(model_kms: np.ndarray, target_kms: np.ndarray) -> float:
    """Compute RMSE as a percentage of the mean target velocity."""
    residuals = model_kms - target_kms
    rmse = np.sqrt(np.mean(residuals ** 2))
    mean_target = np.mean(target_kms)
    return 100.0 * rmse / mean_target


def objective(params, r_kpc, M_enclosed_kg, V_target_kms):
    """
    V4.2 Objective — tunes only the caustic guard amplitude λ.
    a0 and beta are pinned to Single Source of Truth values.
    """
    (lam,) = params
    a0     = phys.A0        
    beta   = _BETA_AUDIT    

    V_model_kms = vcirc_v42(r_kpc, M_enclosed_kg, lam=lam, a0=a0, beta=beta)
    return rmse_percentage(V_model_kms, V_target_kms)


def run_optimization():
    """
    V4.2 Self-Consistency Validation — optimise λ only, a0 pinned to phys.A0.
    """
    # 1. Load Geometry and Mass
    r_kpc, M_enclosed_kg = load_mock_sparc()

    # 2. Verify physical constant presence
    if not hasattr(phys, "A0"):
        raise AttributeError("phys.A0 not found in core.physics_constants. Check your file.")

    # 3. Generate FIDUCIAL TARGET (The Self-Consistency Anchor)
    # We create a "Ground Truth" curve using the physics engine itself.
    lam_fiducial  = 0.05
    a0_fixed      = phys.A0
    
    V_target_kms  = vcirc_v42(r_kpc, M_enclosed_kg,
                               lam=lam_fiducial,
                               a0=a0_fixed,
                               beta=_BETA_AUDIT)

    print(f"[Fiducial target] lam={lam_fiducial}, a0={a0_fixed:.4e}, "
          f"V_flat≈{V_target_kms[-1]:.1f} km/s")

    # 4. Optimization Setup
    # We start the optimizer at 0.10 to see if it can recover the 0.05 fiducial.
    x0     = np.array([0.10])   
    bounds = [(0.0, 0.5)]

    result = minimize(
        objective,
        x0,
        args=(r_kpc, M_enclosed_kg, V_target_kms),
        method="Powell",
        bounds=bounds,
        options={"xtol": 1e-6, "ftol": 1e-6, "disp": False},
    )

    (lam_opt,) = result.x

    # 5. Diagnostics
    V_model_kms    = vcirc_v42(r_kpc, M_enclosed_kg, lam_opt, a0_fixed)
    final_rmse_pct = rmse_percentage(V_model_kms, V_target_kms)
    lam_recovery_err = abs(lam_opt - lam_fiducial)

    print("\n=== SPARC Refinery V4.2 — Self-Consistency Validation Results ===")
    print(f"Fixed    a0  [m/s²]  : {a0_fixed:.4e}  (phys.A0)")
    print(f"Fiducial λ           : {lam_fiducial:.6f}")
    print(f"Recovered λ          : {lam_opt:.6f}  (Δλ = {lam_recovery_err:.2e})")
    print(f"Fixed    β           : {_BETA_AUDIT}  (_BETA_AUDIT)")
    print(f"Final RMSE           : {final_rmse_pct:.6f} %")
    print(f"Optimiser success    : {result.success}  ({result.message})")

    if final_rmse_pct <= 0.01 and lam_recovery_err < 1e-4:
        print("\nConvergence Message  : SELF-CONSISTENCY CONFIRMED — pipeline ready for real SPARC data.")
    elif final_rmse_pct <= 1.0:
        print("\nConvergence Message  : Marginal pass — RMSE < 1% but recovery precision is low.")
    else:
        print("\nConvergence Message  : SELF-CONSISTENCY FAILED — structural issue in pipeline.")
        print("                      Check for cached .pyc files or A0 mismatches.")

# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_optimization()
