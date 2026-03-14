import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import zoom, gaussian_filter
from core.mimetic_engine import MimeticEngine

# --- CONFIGURATION UPDATE ---
FITS_PATH = "data/raw_fits/mastDownload/HST/jbqz31010/jbqz31010_drz.fits"

def run_sentinel_on_fits():
    print(f"--- FINAL CALIBRATION: {FITS_PATH} ---")
    
    with fits.open(FITS_PATH) as hdul:
        data = np.nan_to_num(hdul[1].data, nan=0.0)
    
    data = np.clip(data, 0, np.percentile(data, 99))
    rho_raw = zoom(data, 512 / data.shape[0])
    # Smooth to represent the cluster-scale gas distribution
    rho_baryon = gaussian_filter(rho_raw, sigma=4.0)
    
    # 1. Physical Scale (kg/m^3)
    rho_baryon_phys = (rho_baryon / np.max(rho_baryon)) * 8e-23
    
    # 2. Field Scale (Normalized Q)
    # We target the Notch by setting the peak gradient response to Q=2.0
    gy, gx = np.gradient(rho_baryon / np.max(rho_baryon))
    grad_sq = (gx**2 + gy**2)
    Q_field = (grad_sq / np.max(grad_sq)) * 2.0
    
    # 3. Process
    engine = MimeticEngine(kappa=0.80)
    rho_total = engine.compute_effective_density(rho_baryon_phys, Q_field)
    
    ratio = np.max(rho_total) / np.max(rho_baryon_phys)
    print(f"\n[REPORT] DM Ratio: {ratio:.2f}x")
    print(f"[REPORT] Peak Gravity: {np.max(rho_total):.2e} kg/m^3")
    
    # 4. Save Result
    plt.figure(figsize=(12, 5))
    plt.subplot(1,2,1); plt.title("Baryons (1E 0657-558)"); plt.imshow(rho_baryon_phys, cmap='hot')
    plt.subplot(1,2,2); plt.title("Effective Lensing Map"); plt.imshow(rho_total, cmap='coolwarm')
    plt.savefig("Sentinel_FINAL_REPORT.png")
    print("ANALYSIS LOCKED. Plot saved to Sentinel_FINAL_REPORT.png")
    plt.show()

if __name__ == "__main__":
    run_sentinel_on_fits()
