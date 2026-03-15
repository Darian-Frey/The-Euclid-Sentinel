# debug_bullet.py
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from core.mimetic_engine import MimeticEngine

def run_diagnostic():
    # 1. Direct path to Bullet Cluster
    path = "data/raw_fits/mastDownload/HST/j90701010/j90701010_drz.fits"
    engine = MimeticEngine(kappa=0.80)
    
    print(f"🔭 DIAGNOSTIC: Opening {path}")
    
    with fits.open(path) as hdul:
        # Get actual SCI data (usually HDU 1 for this specific HST file)
        data = hdul[1].data.astype(np.float64)
        data = np.nan_to_num(data)
        
        # 2. CENTER CROP (Crucial for Hubble Mosaics)
        # Instead of the whole empty sky, let's look at the center 1000 pixels
        h, w = data.shape
        data = data[h//2-500:h//2+500, w//2-500:w//2+500]
        
    # 3. Physics Processing
    # Simple gradient for Q-field
    q_field = np.abs(np.gradient(data)[0]) + np.abs(np.gradient(data)[1])
    q_field = (q_field / (np.max(q_field) + 1e-10)) * 10.0
    rho_eff = engine.compute_effective_density(data, q_field)
    
    # 4. ROBUST VISUALIZATION
    fig, ax = plt.subplots(1, 2, figsize=(16, 8), facecolor='black')
    
    # Log scale helps see faint galaxies without blowing out stars
    # We add a tiny constant to avoid log(0)
    baryon_display = np.log10(data - np.min(data) + 1.0)
    mimetic_display = np.log10(rho_eff - np.min(rho_eff) + 1.0)

    ax[0].imshow(baryon_display, cmap='magma', origin='lower')
    ax[0].set_title("Hubble Baryonic Input (Center Crop)", color='white')
    
    ax[1].imshow(mimetic_display, cmap='viridis', origin='lower')
    ax[1].set_title("Mimetic Gravity Potential", color='white')

    for a in ax:
        a.axis('off')

    plt.tight_layout()
    plt.savefig("DIAGNOSTIC_BULLET.png", facecolor='black')
    print("✅ DIAGNOSTIC COMPLETE: Check 'DIAGNOSTIC_BULLET.png'")

if __name__ == "__main__":
    run_diagnostic()
