"""
core/mimetic_engine.py
======================
CODE-GEO V4.2 — Mimetic-Conformal Engine (Physics-Grounded Rebuild)
====================================================================

Changelog vs V3.1.3
--------------------
V3.1.3  DEPRECATED ISSUES:
  - Used free_function_f_q() — the V3.1 Soft-Krylov form (wrong MOND shape,
    RMSE 18.7-23.4% on SPARC, causality violation c_s²_max ≈ 1.04).
  - Hardcoded 12.0 amplification factor in compute_effective_density().
    The "DM ratio" was therefore IMPOSED, not predicted. All Phase VIII
    survey results (~12x universality) were artefacts of this constant.
  - Q-field derived from arbitrary pixel-gradient rescaling normalised to
    peak Q=2.0. No physical derivation.
  - V4.2 functions in core.action were never called by any pipeline tool.

V4.2  THIS VERSION:
  - Imports free_function_derivative_v42() from core.action — the Standard-MOND
    anchor free function with verified c_s² ∈ [0.50, 1.00] and c_T = c.
  - Q is computed from first principles: Q = (g_N / a₀)² where g_N = |∇Φ_N|
    and Φ_N is solved from the baryonic density via FFT Poisson.
  - The effective density rho_eff = rho_baryon / F'(Q) — the DM ratio
    emerges from the physics.
  - DM ratio is now a prediction, not a parameter.

Open Issue (Resolved → Flagged Systematic)
------------------------------------------
RESOLVED [Q-DERIVATION]: Hub-and-spoke session (2026-04) concluded:
  ✓  Thin-disk FFT kernel  Φ(k) = -2πG Σ_b(k)/|k|  is the correct procedure
  ✓  Input must be Σ_b [kg/m²] — NOT ρ [kg/m³], NOT ∇flux, NOT ∇κ
  ✓  Q computed post-deprojection (thin-disk is implicit deprojection)
  ✓  ∇κ and ∇flux are third derivatives of the lensing potential — wrong

REMAINING SYSTEMATIC [CLUSTER-PROJECTION]:
  The thin-disk approximation underestimates line-of-sight mass in triaxial
  merging cluster systems (Bullet, Abell 370, El Gordo). This is a bounded
  systematic, not a free parameter. Quantification pending X-ray / SZ data.

Theory
------
Modified Poisson equation (quasi-Newtonian AQUAL limit):

    ∇·(F'(Q) ∇Φ) = 4πG ρ_baryon

where Q ≡ (|∇Φ_N| / a₀)² = (g_N / a₀)², and Φ_N is the thin-disk potential.

Effective surface density seen by lensing:

    Σ_eff = Σ_b / F'(Q)

Effective acceleration:

    g_eff = g_N / F'(Q)

In the deep-MOND limit (Q ≪ 1):  F'(Q) → √Q → g_eff ∝ √(g_N · a₀)   ✓
In the Newtonian limit  (Q ≫ 1):  F'(Q) → 1   → g_eff = g_N            ✓

References
----------
* Chamseddine & Mukhanov (2013)     — original mimetic gravity
* Bekenstein & Milgrom (1984)       — AQUAL formulation
* core.action V4.2 docstring        — free function derivation and audit
* GW170817 / Abbott et al. (2017)   — c_T = c constraint (satisfied by V4.2)
"""

import numpy as np
import logging
import warnings

# ---------------------------------------------------------------------------
# Internal imports — Single Source of Truth for all constants and functions
# ---------------------------------------------------------------------------
try:
    from core.physics_constants import G, C, A0, LAMBDA, BETA, ELL, KPC_TO_M
    from core.action import (
        free_function_v42,
        free_function_derivative_v42,
        check_stability_v42,
    )
except ImportError as exc:
    raise ImportError(
        "MimeticEngine V4.2 requires core.physics_constants and core.action. "
        f"Import failed: {exc}"
    ) from exc


# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------
logging.basicConfig(
    filename="sentinel_engine.log",
    level=logging.INFO,
    format="%(asctime)s [SENTINEL_V42] %(message)s",
    filemode="a",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# MimeticEngine V4.2
# ---------------------------------------------------------------------------

class MimeticEngine:
    """
    Mimetic-Conformal gravity engine, CODE-GEO V4.2.

    Parameters
    ----------
    lam  : float
        Caustic guard amplitude λ (default: physics_constants.LAMBDA = 0.05).
        This is the sole free parameter of the V4.2 free function. It is
        optimised by sparc_refinery_v4.py against rotation curve data.
    beta : float
        Caustic guard decay rate β (default: physics_constants.BETA = 1.5).
        Fixed per V4.2 audit — do not tune without re-running the stability
        sentinel in core.action.check_stability_v42().
    a0   : float
        MOND acceleration scale [m/s²] (default: physics_constants.A0).
        Pinned to 1.21e-10 m/s². Change in physics_constants.py only.

    Notes
    -----
    The engine does NOT accept kappa. kappa was the V3.1/V4.1 shape parameter
    and has no role in the V4.2 hybrid free function. If you have legacy code
    passing kappa=0.80, remove it — it will raise TypeError here intentionally.
    """

    def __init__(self, lam: float = LAMBDA, beta: float = BETA, a0: float = A0):
        self.lam  = lam
        self.beta = beta
        self.a0   = a0

        log.info(
            f"ENGINE_START_V42: lam={self.lam}, beta={self.beta}, "
            f"a0={self.a0:.4e} m/s²"
        )
        print(
            f"[MimeticEngine V4.2] Initialized | "
            f"λ={self.lam} | β={self.beta} | a₀={self.a0:.3e} m/s²"
        )

    # ------------------------------------------------------------------
    # 1. Q-field derivation
    # ------------------------------------------------------------------

    def _solve_newtonian_potential(
        self,
        sigma_b: np.ndarray,
        dx: float,
    ) -> np.ndarray:
        """
        Solve for the in-plane Newtonian potential of a razor-thin baryonic
        sheet via the 2D thin-disk Poisson kernel (spectral method).

        Physics
        -------
        For a razor-thin mass sheet with surface mass density Σ_b(x,y), the
        in-plane gravitational potential satisfies (in Fourier space):

            Φ(k) = -2π G Σ_b(k) / |k|

        This is the exact result for a test particle in the plane of an
        infinitesimally thin sheet. It differs from the 3D volumetric kernel
        (which uses k² and requires volumetric density ρ [kg/m³]) in two ways:

          (a) Denominator |k| vs k²  — one fewer power of k, reflecting the
              fact that the sheet source is 2D not 3D.
          (b) Prefactor 2π vs 4π   — the Green's function for the 2D-embedded
              3D Laplacian evaluated at z=0.

        Using the 3D kernel on a projected 2D map conflates Σ [kg/m²] with
        ρ [kg/m³] and introduces an implicit line-of-sight depth scale with
        no physical justification. The thin-disk kernel eliminates this.

        Systematic limitation for cluster targets
        -----------------------------------------
        The thin-disk approximation is exact for face-on disk galaxies. For
        merging clusters (Bullet, Abell 370, El Gordo) it underestimates the
        line-of-sight mass contribution in triaxial structures. This is a
        known, bounded systematic — not a free parameter. Quantifying this
        bias for cluster targets is flagged as a follow-on task.

        Parameters
        ----------
        sigma_b : 2-D array  [kg/m²]   baryonic surface mass density
                  Must be in physical units. Convert from FITS flux using:
                  Σ_b = (flux / flux_scale) × (M/L) × Σ_unit
                  where Σ_unit = M_sun / pc² × KPC_TO_M² / MSUN.
                  See tools/euclid_loader.py for the calibration pipeline.
        dx      : float      physical pixel scale [m]
                  Derive from FITS WCS: dx = D_A × pixel_angle_rad
                  where D_A is the angular diameter distance at target redshift.

        Returns
        -------
        phi_N : 2-D array  [m²/s²]  in-plane Newtonian gravitational potential
        """
        ny, nx = sigma_b.shape

        kx = np.fft.fftfreq(nx, d=dx) * 2.0 * np.pi
        ky = np.fft.fftfreq(ny, d=dx) * 2.0 * np.pi
        KX, KY = np.meshgrid(kx, ky)

        # |k| — thin-disk denominator (NOT k², which is the 3D volumetric kernel)
        k_abs = np.sqrt(KX**2 + KY**2)
        k_abs[0, 0] = 1.0e-30   # guard DC mode; mean potential is arbitrary

        sigma_k = np.fft.fft2(sigma_b)

        # Thin-disk Poisson kernel: Φ(k) = -2πG Σ(k) / |k|
        phi_k = -2.0 * np.pi * G * sigma_k / k_abs
        phi_N = np.real(np.fft.ifft2(phi_k))

        log.info(
            f"THIN_DISK_POISSON: grid={ny}x{nx} | dx={dx:.3e} m | "
            f"Σ_peak={np.max(sigma_b):.3e} kg/m² | "
            f"Φ_peak={np.max(phi_N):.3e} m²/s²"
        )
        return phi_N

    def _compute_q_field(
        self,
        sigma_b: np.ndarray,
        dx: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute Q = (g_N / a₀)² from the baryonic surface mass density.

        Procedure (resolved — hub-and-spoke session 2026-04)
        -----------------------------------------------------
        1. Thin-disk Poisson:  Φ(k) = -2πG Σ_b(k) / |k|  →  Φ_N
           (replaces 3D volumetric kernel ∇²Φ = 4πGρ, which mixed units)
        2. In-plane acceleration:  g_N = |∇Φ_N|            →  [m/s²]
        3. Kinetic invariant:      Q = (g_N / a₀)²         →  dimensionless

        Resolved decisions
        ------------------
        ✓  Use baryonic Σ_b (NOT κ gradient, NOT flux gradient).
           ∇κ and ∇flux are third derivatives of the lensing potential —
           wrong physical quantity for g_N (Reply 3, hub-and-spoke session).
        ✓  g_N derived via Poisson (non-local integral), not local gradient
           of surface density — local ∇Σ has wrong dimensions and scaling.
        ✓  Q computed post-deprojection (thin-disk = implicit deprojection
           for face-on maps; inclination correction applied upstream).

        Remaining open issue — cluster targets
        ---------------------------------------
        TODO [CLUSTER-PROJECTION]: The thin-disk kernel is exact for face-on
        disk galaxies. For merging clusters (Bullet, Abell 370, El Gordo) the
        line-of-sight mass distribution is triaxial and not well-approximated
        by a thin sheet. The systematic bias in g_N (and therefore Q) for
        cluster targets is unquantified. Candidate approaches:
          (a) X-ray temperature profile → hydrostatic g_N independent of Σ_b
          (b) SZ decrement → independent pressure/mass proxy
          (c) Weak lensing shear map → total (not baryonic) g_N upper bound
        This systematic does not affect SPARC galaxy-scale fits.

        Parameters
        ----------
        sigma_b : 2-D array  [kg/m²]   baryonic surface mass density
        dx      : float      physical pixel scale [m]

        Returns
        -------
        Q   : 2-D array  [dimensionless]  kinetic invariant field
        g_N : 2-D array  [m/s²]          Newtonian acceleration field
                                          (returned for diagnostics/telemetry)
        """
        phi_N = self._solve_newtonian_potential(sigma_b, dx)

        # In-plane Newtonian acceleration: g_N = |∇Φ_N|
        # np.gradient uses central differences — appropriate for smooth FFT output
        gx, gy = np.gradient(phi_N, dx)
        g_N = np.sqrt(gx**2 + gy**2)   # [m/s²]

        Q = (g_N / self.a0) ** 2       # dimensionless

        log.info(
            f"Q_FIELD: Q_mean={np.mean(Q):.4e} | Q_max={np.max(Q):.4e} | "
            f"g_N_mean={np.mean(g_N):.4e} m/s² | g_N_max={np.max(g_N):.4e} m/s²"
        )
        return Q, g_N

    # ------------------------------------------------------------------
    # 2. Effective density
    # ------------------------------------------------------------------

    def compute_effective_density(
        self,
        sigma_b: np.ndarray,
        dx: float,
    ) -> tuple[np.ndarray, np.ndarray, dict]:
        """
        Compute the effective lensing surface density Σ_eff = Σ_b / F'(Q).

        The DM-equivalent ratio emerges from the physics:

            DM ratio = <Σ_eff> / <Σ_b> = <1 / F'(Q)>

        In the deep-MOND limit (Q ≪ 1):  F'(Q) ≈ √Q ≪ 1  →  ratio ≫ 1
        In the Newtonian limit  (Q ≫ 1):  F'(Q) → 1         →  ratio → 1

        Parameters
        ----------
        sigma_b : 2-D array  [kg/m²]
            Baryonic surface mass density. Must be in physical units —
            convert from FITS flux via M/L ratio and distance before calling.
            Do NOT pass raw pixel counts or normalised flux.
        dx : float
            Physical pixel scale [m]. Extract from FITS WCS header:
            dx = D_A(z) × pixel_angular_size_rad.
            Use tools/euclid_loader.py for calibrated extraction.

        Returns
        -------
        sigma_eff : 2-D array  [kg/m²]   effective lensing surface density
        F_prime   : 2-D array  [–]       MOND interpolation field F'(Q)
        telemetry : dict                  diagnostic quantities for dashboard
        """
        # Step 1: Q-field and g_N from thin-disk Poisson
        Q, g_N = self._compute_q_field(sigma_b, dx)

        # Step 2: V4.2 MOND interpolation field F'(Q)
        F_prime = free_function_derivative_v42(Q, lam=self.lam, beta=self.beta)

        # Guard against F'(Q) ≤ 0 (should not occur for V4.2, but be safe)
        F_prime_safe = np.where(F_prime > 1.0e-10, F_prime, 1.0e-10)
        if np.any(F_prime <= 1.0e-10):
            warnings.warn(
                "F'(Q) ≤ 0 encountered — check Q-field range and V4.2 parameters.",
                RuntimeWarning,
                stacklevel=2,
            )

        # Step 3: Effective surface density — ratio emerges from physics
        sigma_eff = sigma_b / F_prime_safe

        # Step 4: Telemetry — g_N stats now included for physical auditing
        dm_ratio_mean = np.mean(sigma_eff) / (np.mean(sigma_b) + 1.0e-40)
        dm_ratio_peak = np.max(sigma_eff)  / (np.max(sigma_b)  + 1.0e-40)

        telemetry = {
            "Q_mean"        : float(np.mean(Q)),
            "Q_max"         : float(np.max(Q)),
            "g_N_mean_ms2"  : float(np.mean(g_N)),
            "g_N_max_ms2"   : float(np.max(g_N)),
            "F_prime_mean"  : float(np.mean(F_prime)),
            "F_prime_min"   : float(np.min(F_prime)),
            "dm_ratio_mean" : float(dm_ratio_mean),
            "dm_ratio_peak" : float(dm_ratio_peak),
            "lam"           : self.lam,
            "beta"          : self.beta,
            "a0"            : self.a0,
        }

        log.info(
            f"STABILITY_METRIC: Q_mean={telemetry['Q_mean']:.4e} | "
            f"Q_max={telemetry['Q_max']:.4e} | "
            f"g_N_max={telemetry['g_N_max_ms2']:.4e} m/s² | "
            f"F'_mean={telemetry['F_prime_mean']:.4f} | "
            f"Ratio_mean={telemetry['dm_ratio_mean']:.2f}x | "
            f"Ratio_peak={telemetry['dm_ratio_peak']:.2f}x"
        )

        print(
            f"  [Engine] g_N_max={telemetry['g_N_max_ms2']:.3e} m/s² | "
            f"Q_max={telemetry['Q_max']:.3e} | "
            f"F'(Q) ∈ [{telemetry['F_prime_min']:.4f}, 1.000] | "
            f"DM ratio (mean): {telemetry['dm_ratio_mean']:.2f}x | "
            f"(peak): {telemetry['dm_ratio_peak']:.2f}x"
        )

        return sigma_eff, F_prime, telemetry

    # ------------------------------------------------------------------
    # 3. Lensing potential (unchanged — FFT Poisson on effective density)
    # ------------------------------------------------------------------

    def get_lensing_potential(
        self,
        sigma_eff: np.ndarray,
        dx: float,
    ) -> np.ndarray:
        """
        Solve for the effective lensing potential Φ_eff via thin-disk Poisson.

        Called AFTER compute_effective_density(). The lensing potential is
        sourced by the mimetic-enhanced surface density Σ_eff = Σ_b / F'(Q).

        Parameters
        ----------
        sigma_eff : 2-D array  [kg/m²]  from compute_effective_density()
        dx        : float      pixel scale [m]

        Returns
        -------
        phi_eff : 2-D array  [m²/s²]
        """
        return self._solve_newtonian_potential(sigma_eff, dx)

    # ------------------------------------------------------------------
    # 4. Stability check (expose V4.2 sentinel to pipeline tools)
    # ------------------------------------------------------------------

    def run_stability_check(self, Q: np.ndarray):
        """
        Run the V4.2 Stability Sentinel over the Q-field.

        Checks c_s² ∈ [0, 1] (ghost-free and subluminal) and c_T = c
        (conformal structure). Wraps core.action.check_stability_v42().

        Parameters
        ----------
        Q : 2-D or 1-D array

        Returns
        -------
        result : StabilityResult
        """
        return check_stability_v42(Q.ravel(), lam=self.lam, beta=self.beta)


# ---------------------------------------------------------------------------
# Quick self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print("=" * 55)
    print("MimeticEngine V4.2 — Self-Test (Thin-Disk Kernel)")
    print("=" * 55)

    # Synthetic baryonic map: two Gaussian blobs in SURFACE MASS DENSITY [kg/m²]
    # Representative of a galaxy cluster at z~0.3, ~650 kpc field
    GRID = 256
    BOX  = 2.0e22   # ~650 kpc physical extent
    dx   = BOX / GRID

    x = np.linspace(-BOX / 2, BOX / 2, GRID)
    X, Y = np.meshgrid(x, x)

    # Surface mass density [kg/m²] — physically motivated for a cluster
    # Peak ~ 1e9 M_sun/kpc² ~ 2e-14 kg/m²
    sigma_a = 2.0e-14 * np.exp(-(X**2 + Y**2) / (0.08 * BOX) ** 2)
    sigma_b_map = 1.2e-14 * np.exp(-((X - 0.2 * BOX)**2 + Y**2) / (0.05 * BOX) ** 2)
    sigma_b = sigma_a + sigma_b_map

    engine = MimeticEngine()

    # Stability check over expected Q range
    Q_test = np.logspace(-4, 3, 500)
    stability = engine.run_stability_check(Q_test)
    print(f"\nStability sentinel: {stability}")

    # Effective surface density
    sigma_eff, F_prime, telemetry = engine.compute_effective_density(sigma_b, dx)

    print("\n--- Telemetry ---")
    for k, v in telemetry.items():
        print(f"  {k:>18s}: {v}")

    # Sanity check: g_N should be in the MOND regime for outer cluster regions
    # a0 = 1.21e-10 m/s²; clusters have g_N ~ 1e-10 to 1e-9 m/s² in core
    g_N_check = telemetry["g_N_max_ms2"]
    regime = "Newtonian" if g_N_check > 10 * engine.a0 else "MOND-transition"
    print(f"\n  g_N_max / a₀ = {g_N_check / engine.a0:.2f}  → {regime} regime at peak")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].imshow(np.log10(sigma_b + 1e-40), cmap="hot", origin="lower")
    axes[0].set_title("log₁₀ Σ_b  [kg/m²]")

    im1 = axes[1].imshow(F_prime, cmap="viridis", origin="lower", vmin=0.0, vmax=1.0)
    axes[1].set_title("F'(Q)  [MOND interpolation field]")
    plt.colorbar(im1, ax=axes[1])

    axes[2].imshow(np.log10(sigma_eff + 1e-40), cmap="coolwarm", origin="lower")
    axes[2].set_title("log₁₀ Σ_eff  [kg/m²]")

    plt.suptitle(
        f"MimeticEngine V4.2 — Thin-Disk Kernel  |  "
        f"DM ratio (mean): {telemetry['dm_ratio_mean']:.2f}x  |  "
        f"g_N_max: {telemetry['g_N_max_ms2']:.2e} m/s²",
        fontsize=10,
    )
    plt.tight_layout()
    plt.savefig("mimetic_engine_v42_selftest.png", dpi=150)
    print("\nPlot → mimetic_engine_v42_selftest.png")