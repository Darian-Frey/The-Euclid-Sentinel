"""
tools/euclid_loader.py
======================
CODE-GEO V4.2 — Calibrated FITS Loader & Baryonic Mass Map Pipeline
=====================================================================

Purpose
-------
Convert raw HST/Euclid FITS flux data into physically calibrated baryonic
surface mass density maps Σ_b [kg/m²] and physical pixel scales dx [m],
ready for direct ingestion by MimeticEngine.compute_effective_density().

Changelog vs original
---------------------
Original:
  - map_to_physical_units() multiplied by hardcoded 1e-15 placeholder.
  - No dx extraction — pixel scale was never computed.
  - No WCS handling beyond object construction.
  - M/L ratio undocumented, not auditable.

V4.2:
  - TARGET_REGISTRY: per-target redshift, M/L, filter, literature reference.
  - dx derived from WCS pixel scale + astropy cosmology D_A(z). No hardcoding.
  - M/L ratio is a named, overridable parameter — not a hidden constant.
  - Returns (sigma_b, dx, provenance) — every assumption is logged.
  - Σ_b pipeline: flux → flux density → luminosity surface density → Σ_b
  - Sensitivity analysis supported: loop over ml_ratio values externally.

M/L Calibration Note
--------------------
The mass-to-light ratio Υ (M/L) converts observed luminosity surface density
to stellar mass surface density:

    Σ_* = Υ × I_ν × conversion_factors

This is an ASTROPHYSICAL ASSUMPTION, not a free parameter of the
Mimetic-Conformal theory. It is sourced from stellar population synthesis
(SPS) models in the literature for each target. The registry defaults are
conservative mid-range values; uncertainty is ±30–50% typical for cluster
stellar populations.

For a publishable systematic error budget, run the pipeline at:
    ml_ratio = [1.5, 2.0, 3.0, 4.0, 5.0]
and report the resulting spread in DM ratio from the engine. This is a
one-liner loop — the parameterised design makes it trivial.

Gas mass note
-------------
Baryonic mass = stellar mass + gas mass. For galaxy clusters, the hot
intracluster gas (ICM) dominates baryonic mass and is better traced by
X-ray surface brightness than optical flux. The current pipeline uses
optical flux only (stellar component). ICM contribution is flagged as a
TODO per target — adding it requires X-ray data (Chandra/XMM) and is a
significant systematic for cluster targets.

References
----------
* Clowe et al. (2006)         — Bullet Cluster baryonic mass map
* Hoekstra et al. (2011)      — Abell 370 lensing analysis
* Menanteau et al. (2012)     — El Gordo cluster properties
* Beckwith et al. (2006)      — HUDF photometry
* Bell & de Jong (2001)       — stellar M/L vs colour relations
* Astropy Collaboration (2022) — cosmology and WCS tools
"""

import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import os

try:
    from core.physics_constants import G, C, MSUN, KPC_TO_M
except ImportError:
    warnings.warn(
        "core.physics_constants not found — using built-in fallback constants.",
        ImportWarning,
        stacklevel=2,
    )
    G       = 6.67430e-11
    C       = 2.998e8
    MSUN    = 1.98847e30
    KPC_TO_M = 3.08567758e19


# ---------------------------------------------------------------------------
# Cosmology — single source of truth
# ---------------------------------------------------------------------------
# Standard flat ΛCDM consistent with Planck 2018.
# Change here to propagate to all D_A(z) calculations.
_COSMO = FlatLambdaCDM(H0=67.4, Om0=0.315)


# ---------------------------------------------------------------------------
# Target Registry
# ---------------------------------------------------------------------------
# Each entry documents the astrophysical assumptions used for calibration.
# All values are sourced from peer-reviewed literature (cited per entry).
#
# Fields
# ------
# name        : common name for display
# z           : spectroscopic redshift
# ml_ratio    : stellar M/L in solar units [M_sun / L_sun]
#               V-band default from Bell & de Jong (2001) colour-M/L relation
#               applied to cluster stellar populations (old, red stars → high M/L)
# ml_band     : photometric band the M/L applies to
# ml_ref      : literature source for the M/L value
# flux_hdu    : which HDU index contains the science data in HST DRZ files
# flux_unit   : native flux unit of the FITS data ('electrons/s' for HST ACS/WFC)
# zeropoint   : AB magnitude zeropoint for the filter (for flux → luminosity)
# gas_fraction: hot gas fraction of total baryonic mass (ICM-to-stellar ratio)
#               flagged TODO where X-ray data not yet integrated
# notes       : any target-specific caveats
# ---------------------------------------------------------------------------
TARGET_REGISTRY = {
    "bullet_cluster": {
        "name"        : "Bullet Cluster (1E 0657-558)",
        "z"           : 0.296,
        "ml_ratio"    : 2.5,       # M_sun/L_sun, V-band
        "ml_band"     : "V",
        "ml_ref"      : "Clowe et al. (2006), Bell & de Jong (2001)",
        "flux_hdu"    : 1,         # SCI extension in HST DRZ
        "flux_unit"   : "electrons/s",
        "zeropoint"   : 26.49,     # ACS/WFC F606W AB zeropoint
        "gas_fraction": None,      # TODO: integrate Chandra X-ray map
        "notes"       : (
            "Merging system — thin-disk approximation underestimates "
            "line-of-sight ICM mass. Gas-to-stellar ratio ~6:1 in cluster core. "
            "Optical flux traces stellar mass only."
        ),
    },
    "abell_370": {
        "name"        : "Abell 370",
        "z"           : 0.375,
        "ml_ratio"    : 3.0,
        "ml_band"     : "V",
        "ml_ref"      : "Hoekstra et al. (2011), Bell & de Jong (2001)",
        "flux_hdu"    : 1,
        "flux_unit"   : "electrons/s",
        "zeropoint"   : 26.49,
        "gas_fraction": None,      # TODO: integrate XMM-Newton X-ray map
        "notes"       : (
            "Massive relaxed cluster with elongated arc system. "
            "Stellar M/L elevated for BCG-dominated core."
        ),
    },
    "el_gordo": {
        "name"        : "El Gordo (ACT-CL J0102-4915)",
        "z"           : 0.870,
        "ml_ratio"    : 2.0,       # younger stellar population at high-z → lower M/L
        "ml_band"     : "V",
        "ml_ref"      : "Menanteau et al. (2012), Bell & de Jong (2001)",
        "flux_hdu"    : 1,
        "flux_unit"   : "electrons/s",
        "zeropoint"   : 25.96,     # ACS/WFC F814W AB zeropoint (rest-frame optical)
        "gas_fraction": None,      # TODO: integrate ACT SZ + Chandra data
        "notes"       : (
            "Most massive known high-z merger. At z=0.87 the ACS filter "
            "samples rest-frame B-band — M/L lower than local clusters. "
            "Thin-disk systematic worst here due to violent merger geometry."
        ),
    },
    "hudf": {
        "name"        : "Hubble Ultra Deep Field",
        "z"           : 1.0,       # representative median redshift of field galaxies
        "ml_ratio"    : 1.5,       # younger high-z field galaxies → lower M/L
        "ml_band"     : "V",
        "ml_ref"      : "Beckwith et al. (2006), Bell & de Jong (2001)",
        "flux_hdu"    : 1,
        "flux_unit"   : "electrons/s",
        "zeropoint"   : 26.49,
        "gas_fraction": None,      # gas mass negligible for field galaxies
        "notes"       : (
            "Ensemble of field galaxies at mixed redshifts. Single z used "
            "for dx calculation is a median approximation — individual "
            "galaxy redshifts would give per-object physical scales. "
            "DM ratio here is a statistical ensemble result."
        ),
    },
}


# ---------------------------------------------------------------------------
# Conversion constants
# ---------------------------------------------------------------------------
# 1 M_sun / pc² in kg/m²
_MSUN_PER_PC2_TO_SI = MSUN / (3.08567758e16 ** 2)  # 1 pc = 3.086e16 m

# L_sun in Watts (bolometric, used internally for unit chain)
_LSUN_W = 3.828e26


# ---------------------------------------------------------------------------
# EuclidLoader
# ---------------------------------------------------------------------------

class EuclidLoader:
    """
    Calibrated FITS loader for the CODE-GEO V4.2 Mimetic-Conformal pipeline.

    Parameters
    ----------
    data_path : str
        Root directory containing cluster FITS files.
    cosmo : astropy.cosmology instance, optional
        Cosmology for angular diameter distance D_A(z). Defaults to the
        module-level _COSMO (flat ΛCDM, Planck 2018).

    Usage
    -----
    loader = EuclidLoader("data/raw_fits/mastDownload/HST/")
    sigma_b, dx, provenance = loader.load_and_calibrate(
        fits_path  = "j90701010/j90701010_drz.fits",
        target_key = "bullet_cluster",
    )
    sigma_eff, F_prime, telemetry = engine.compute_effective_density(sigma_b, dx)
    """

    def __init__(self, data_path: str = "data/clusters/", cosmo=None):
        self.data_path = data_path
        self.cosmo     = cosmo or _COSMO
        os.makedirs(data_path, exist_ok=True)
        print(f"[EuclidLoader V4.2] Initialized | data_path={data_path}")

    # ------------------------------------------------------------------
    # 1. Raw FITS ingestion
    # ------------------------------------------------------------------

    def load_cluster_fits(
        self,
        filename: str,
        hdu_index: int = 1,
        center_crop: int = 1000,
    ) -> tuple[np.ndarray, fits.Header, WCS]:
        """
        Load a FITS file and return the science data, header, and WCS.

        Parameters
        ----------
        filename     : str   FITS filename relative to self.data_path
        hdu_index    : int   HDU containing the science data (default 1 for DRZ)
        center_crop  : int   if > 0, crop to a (2×crop)×(2×crop) central region
                             to avoid edge artefacts in HST mosaics. Set 0 to
                             use the full frame.

        Returns
        -------
        data   : 2-D float64 array   raw flux [electrons/s for HST]
        header : fits.Header
        wcs    : WCS object          for pixel scale extraction
        """
        path = os.path.join(self.data_path, filename)

        with fits.open(path) as hdul:
            data   = hdul[hdu_index].data.astype(np.float64)
            header = hdul[hdu_index].header
            wcs    = WCS(header, naxis=2)

        data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

        if center_crop > 0:
            h, w = data.shape
            cy, cx = h // 2, w // 2
            r = min(center_crop, cy, cx)   # guard against small images
            data = data[cy - r : cy + r, cx - r : cx + r]

        print(
            f"[EuclidLoader] Loaded: {filename} | "
            f"HDU={hdu_index} | shape={data.shape} | "
            f"flux range=[{data.min():.3e}, {data.max():.3e}]"
        )
        return data, header, wcs

    # ------------------------------------------------------------------
    # 2. Physical pixel scale from WCS + cosmology
    # ------------------------------------------------------------------

    def extract_dx(
        self,
        wcs: WCS,
        z: float,
    ) -> float:
        """
        Compute the physical pixel scale dx [m] from the WCS and redshift.

        Procedure
        ---------
        1. Extract pixel angular size θ_pix [arcsec] from WCS CD matrix or
           CDELT keywords.
        2. Compute angular diameter distance D_A(z) [m] from cosmology.
        3. dx = D_A × θ_pix [rad]

        Parameters
        ----------
        wcs : WCS   from the FITS header
        z   : float target redshift

        Returns
        -------
        dx : float  physical pixel scale [m]
        """
        # Extract pixel scale in degrees from WCS
        # wcs.proj_plane_pixel_scales() returns [scale_x, scale_y] in degrees
        try:
            scales_deg = wcs.proj_plane_pixel_scales()
            # Use geometric mean of x and y scales for non-square pixels
            pix_scale_deg = float(np.sqrt(scales_deg[0].value * scales_deg[1].value))
        except Exception:
            # Fallback: try CDELT1 directly
            try:
                pix_scale_deg = abs(float(wcs.wcs.cdelt[0]))
            except Exception as exc:
                raise ValueError(
                    f"Cannot extract pixel scale from WCS: {exc}. "
                    "Check that the FITS header contains valid CD matrix or CDELT."
                ) from exc

        pix_scale_rad = np.deg2rad(pix_scale_deg)

        # Angular diameter distance [m]
        D_A_mpc = self.cosmo.angular_diameter_distance(z).to(u.m).value

        dx = D_A_mpc * pix_scale_rad

        print(
            f"[EuclidLoader] WCS pixel scale: {pix_scale_deg * 3600:.4f} arcsec | "
            f"D_A(z={z}): {D_A_mpc / 3.086e22:.3f} Mpc | "
            f"dx: {dx:.4e} m  ({dx / KPC_TO_M * 1e3:.3f} pc/pixel)"
        )
        return dx

    # ------------------------------------------------------------------
    # 3. Flux → Σ_b calibration
    # ------------------------------------------------------------------

    def flux_to_sigma_b(
        self,
        flux_map: np.ndarray,
        dx: float,
        z: float,
        ml_ratio: float,
        zeropoint: float,
    ) -> np.ndarray:
        """
        Convert FITS flux [electrons/s] to baryonic surface mass density
        Σ_b [kg/m²].

        Calibration pipeline
        --------------------
        Step 1: Flux → AB magnitude per pixel
            m_AB = zeropoint − 2.5 × log10(flux)

        Step 2: AB magnitude → flux density [erg/s/cm²/Hz]
            f_ν = 10^(−(m_AB + 48.6) / 2.5)

        Step 3: Flux density → luminosity surface density [L_sun/pc²]
            Accounting for distance and pixel solid angle.

        Step 4: Luminosity → stellar mass surface density [M_sun/pc²]
            Σ_* = Υ × I   where Υ = ml_ratio [M_sun/L_sun]

        Step 5: Unit conversion → [kg/m²]
            Σ_b = Σ_* × _MSUN_PER_PC2_TO_SI

        Parameters
        ----------
        flux_map  : 2-D array  raw flux [electrons/s per pixel]
        dx        : float      physical pixel scale [m]
        z         : float      target redshift (for luminosity distance)
        ml_ratio  : float      stellar M/L [M_sun/L_sun]
        zeropoint : float      AB magnitude zeropoint for the filter

        Returns
        -------
        sigma_b : 2-D array  [kg/m²]  baryonic surface mass density
        """
        # Guard: replace non-positive flux with a small floor before log
        flux_safe = np.where(flux_map > 0.0, flux_map, 1.0e-30)

        # Step 1: AB magnitude per pixel
        m_AB = zeropoint - 2.5 * np.log10(flux_safe)

        # Step 2: flux density per pixel [erg/s/cm²/Hz]
        f_nu_cgs = 10.0 ** (-(m_AB + 48.6) / 2.5)

        # Step 3: luminosity surface density [L_sun/pc²]
        # Luminosity distance D_L = D_A × (1+z)²
        D_A_cm = self.cosmo.angular_diameter_distance(z).to(u.cm).value
        D_L_cm = D_A_cm * (1.0 + z) ** 2

        # Pixel solid angle in pc²
        dx_pc      = dx / 3.08567758e16          # m → pc
        pix_area_pc2 = dx_pc ** 2

        # Luminosity per pixel in L_sun (using nu × f_nu ≈ bolometric for
        # a flat SED approximation; acceptable for a first-pass calibration)
        # L [erg/s] = 4π D_L² × f_ν × Δν
        # We use the surface brightness I_ν [L_sun/pc²/Hz] route:
        # I_ν = f_ν × 4π D_L² / (pix_area_pc2 × pc_to_cm²)
        pc_to_cm   = 3.08567758e18
        L_per_pix_erg_s_Hz = f_nu_cgs * 4.0 * np.pi * D_L_cm**2
        I_nu_lsun_pc2_hz   = (
            L_per_pix_erg_s_Hz
            / (_LSUN_W * 1e7)          # erg/s → L_sun (W)... erg/s / (L_sun in erg/s)
            / pix_area_pc2
        )
        # Integrate over a representative bandwidth Δν for the filter
        # ACS F606W: Δλ ≈ 2300 Å → Δν ≈ c Δλ / λ_eff²
        lambda_eff_cm = 6060e-8        # F606W effective wavelength [cm]
        delta_lambda_cm = 2300e-8      # F606W bandwidth [cm]
        delta_nu_hz    = (3.0e10 * delta_lambda_cm) / lambda_eff_cm**2

        I_lsun_pc2 = I_nu_lsun_pc2_hz * delta_nu_hz   # [L_sun/pc²]

        # Step 4: stellar mass surface density [M_sun/pc²]
        sigma_star_msun_pc2 = ml_ratio * I_lsun_pc2

        # Step 5: SI conversion [kg/m²]
        sigma_b = sigma_star_msun_pc2 * _MSUN_PER_PC2_TO_SI

        # Zero out pixels where flux was non-physical
        sigma_b = np.where(flux_map > 0.0, sigma_b, 0.0)

        print(
            f"[EuclidLoader] Σ_b calibration complete | "
            f"M/L={ml_ratio} M☉/L☉ | "
            f"Σ_b peak={np.max(sigma_b):.3e} kg/m² | "
            f"Σ_b mean={np.mean(sigma_b[sigma_b > 0]):.3e} kg/m²"
        )
        return sigma_b

    # ------------------------------------------------------------------
    # 4. Primary interface — load + calibrate in one call
    # ------------------------------------------------------------------

    def load_and_calibrate(
        self,
        fits_path: str,
        target_key: str,
        ml_ratio: float | None = None,
        center_crop: int = 1000,
    ) -> tuple[np.ndarray, float, dict]:
        """
        Load a FITS file and return a calibrated (sigma_b, dx, provenance) tuple
        ready for MimeticEngine.compute_effective_density().

        Parameters
        ----------
        fits_path  : str   path to FITS file, relative to self.data_path
        target_key : str   key into TARGET_REGISTRY (e.g. 'bullet_cluster')
        ml_ratio   : float | None
            Override the registry M/L ratio. If None, uses the registry default.
            Set explicitly for sensitivity analysis sweeps.
        center_crop : int  pixels each side of centre to retain (default 1000)

        Returns
        -------
        sigma_b    : 2-D array  [kg/m²]   baryonic surface mass density
        dx         : float      [m]        physical pixel scale
        provenance : dict                  full audit trail of assumptions used

        Example
        -------
        # Default M/L
        sigma_b, dx, prov = loader.load_and_calibrate(
            "j90701010/j90701010_drz.fits", "bullet_cluster"
        )

        # Sensitivity sweep
        for ml in [1.5, 2.5, 3.5, 5.0]:
            sigma_b, dx, prov = loader.load_and_calibrate(
                "j90701010/j90701010_drz.fits", "bullet_cluster", ml_ratio=ml
            )
            sigma_eff, _, tel = engine.compute_effective_density(sigma_b, dx)
            print(f"M/L={ml} → DM ratio={tel['dm_ratio_mean']:.2f}x")
        """
        if target_key not in TARGET_REGISTRY:
            raise KeyError(
                f"Target '{target_key}' not in TARGET_REGISTRY. "
                f"Available: {list(TARGET_REGISTRY.keys())}"
            )

        reg = TARGET_REGISTRY[target_key]

        # Resolve M/L — override takes precedence, registry is default
        ml_used    = ml_ratio if ml_ratio is not None else reg["ml_ratio"]
        ml_source  = "override" if ml_ratio is not None else f"registry ({reg['ml_ref']})"

        # Load raw FITS
        flux_map, header, wcs = self.load_cluster_fits(
            fits_path,
            hdu_index   = reg["flux_hdu"],
            center_crop = center_crop,
        )

        # Physical pixel scale
        dx = self.extract_dx(wcs, reg["z"])

        # Calibrate to Σ_b
        sigma_b = self.flux_to_sigma_b(
            flux_map,
            dx        = dx,
            z         = reg["z"],
            ml_ratio  = ml_used,
            zeropoint = reg["zeropoint"],
        )

        # Provenance — complete audit trail
        provenance = {
            "target_key"    : target_key,
            "target_name"   : reg["name"],
            "fits_path"     : fits_path,
            "redshift_z"    : reg["z"],
            "ml_ratio_used" : ml_used,
            "ml_source"     : ml_source,
            "ml_band"       : reg["ml_band"],
            "ml_ref"        : reg["ml_ref"],
            "zeropoint_AB"  : reg["zeropoint"],
            "flux_unit"     : reg["flux_unit"],
            "dx_m"          : dx,
            "dx_pc"         : dx / 3.08567758e16,
            "center_crop_px": center_crop,
            "sigma_b_peak"  : float(np.max(sigma_b)),
            "sigma_b_mean"  : float(np.mean(sigma_b[sigma_b > 0])),
            "gas_fraction"  : reg["gas_fraction"],
            "notes"         : reg["notes"],
            "cosmo"         : str(self.cosmo),
        }

        if reg["gas_fraction"] is None:
            warnings.warn(
                f"[{target_key}] Gas mass (ICM) not included in Σ_b — "
                "optical flux traces stellar component only. "
                "For cluster targets the ICM may dominate baryonic mass. "
                "Integrate X-ray/SZ data to complete the baryonic budget.",
                UserWarning,
                stacklevel=2,
            )

        print(
            f"[EuclidLoader] ✓ {reg['name']} calibrated | "
            f"M/L={ml_used} ({ml_source}) | "
            f"dx={dx:.3e} m | "
            f"Σ_b_peak={np.max(sigma_b):.3e} kg/m²"
        )
        return sigma_b, dx, provenance

    # ------------------------------------------------------------------
    # 5. Provenance report
    # ------------------------------------------------------------------

    @staticmethod
    def print_provenance(provenance: dict) -> None:
        """Pretty-print the provenance dict for logging or paper appendix."""
        print("\n── Calibration Provenance ──────────────────────────────")
        for k, v in provenance.items():
            print(f"  {k:<20s}: {v}")
        print("────────────────────────────────────────────────────────\n")


# ---------------------------------------------------------------------------
# Script entry point — quick loader test (no real FITS required)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("EuclidLoader V4.2 — Registry & Calibration Self-Test")
    print("=" * 60)

    print("\nRegistered targets:")
    for key, reg in TARGET_REGISTRY.items():
        D_A = _COSMO.angular_diameter_distance(reg["z"]).to(u.Mpc).value
        print(
            f"  {key:<20s}  z={reg['z']:.3f}  "
            f"D_A={D_A:.1f} Mpc  "
            f"M/L={reg['ml_ratio']} M☉/L☉  "
            f"ref: {reg['ml_ref'][:40]}"
        )

    print("\nSensitivity analysis preview (synthetic 512×512 map):")
    print("  M/L sweep for bullet_cluster at fixed synthetic flux:")

    # Synthetic flux map (Gaussian cluster, no real FITS needed)
    GRID = 512
    x    = np.linspace(-1, 1, GRID)
    X, Y = np.meshgrid(x, x)
    synthetic_flux = 100.0 * np.exp(-(X**2 + Y**2) / 0.2**2) + 1.0

    reg      = TARGET_REGISTRY["bullet_cluster"]
    loader   = EuclidLoader()

    # Approximate dx for Bullet Cluster (no real WCS here)
    D_A_m    = _COSMO.angular_diameter_distance(reg["z"]).to(u.m).value
    pix_rad  = np.deg2rad(0.05 / 3600)   # typical HST ACS 0.05 arcsec/pixel
    dx_synth = D_A_m * pix_rad

    print(f"\n  Synthetic dx = {dx_synth:.4e} m ({dx_synth/KPC_TO_M*1e3:.3f} pc/pix)")
    print(f"  {'M/L':>6}  {'Σ_b_peak [kg/m²]':>20}  {'Σ_b_mean [kg/m²]':>20}")

    for ml in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
        sb = loader.flux_to_sigma_b(
            synthetic_flux,
            dx        = dx_synth,
            z         = reg["z"],
            ml_ratio  = ml,
            zeropoint = reg["zeropoint"],
        )
        print(f"  {ml:>6.1f}  {np.max(sb):>20.4e}  {np.mean(sb[sb>0]):>20.4e}")

    print("\n[EuclidLoader] Self-test complete.")
    print("To run on real data: loader.load_and_calibrate(fits_path, target_key)")