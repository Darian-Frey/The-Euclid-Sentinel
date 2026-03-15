# tools/run_full_survey.py
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import zoom, gaussian_filter

# Ensure 'core' is accessible
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.mimetic_engine import MimeticEngine

def run_survey():
    # Verified Path Registry
    survey_targets = {
        "Bullet_Cluster": "data/raw_fits/mastDownload/HST/j90701010/j90701010_drz.fits",
        "Abell_370": "data/raw_fits/mastDownload/HST/jabu01030/jabu01030_drz.fits",
        "El_Gordo": "data/raw_fits/mastDownload/HST/jbqz31010/jbqz31010_drz.fits",
        "HUDF_DeepField": "data/raw_fits/mastDownload/HST/j8wc7c010/j8wc7c010_drz.fits"
    }

    engine = MimeticEngine(kappa=0.80)
    print("🛰️ EUCLID SENTINEL: EXECUTING CALIBRATED PHYSICS BATCH\n" + "="*55)

    for name, path in survey_targets.items():
        if not os.path.exists(path): continue
        print(f"🔭 PROCESSING: {name}...")
        
        try:
            with fits.open(path) as hdul:
                # Target the SCI layer (HDU 1)
                raw_data = np.nan_to_num(hdul[1].data, nan=0.0)
            
            # 1. CROP & ZOOM (Mimics fits_to_mimetic logic)
            # Center crop first to avoid edge artifacts
            h, w = raw_data.shape
            crop = raw_data[h//2-1000:h//2+1000, w//2-1000:w//2+1000]
            
            # Clip and downsample to 512x512 for smooth field resolution
            crop = np.clip(crop, 0, np.percentile(crop, 99))
            rho_raw = zoom(crop, 512 / crop.shape[0])
            
            # 2. BARYONIC SMOOTHING
            rho_baryon = gaussian_filter(rho_raw, sigma=4.0)
            rho_baryon_phys = (rho_baryon / (np.max(rho_baryon) + 1e-10)) * 8e-23
            
            # 3. Q-FIELD GRADIENT (Corrected to Q=2.0 Peak)
            gy, gx = np.gradient(rho_baryon / (np.max(rho_baryon) + 1e-10))
            grad_sq = (gx**2 + gy**2)
            Q_field = (grad_sq / (np.max(grad_sq) + 1e-10)) * 2.0
            
            # 4. ENGINE EXECUTION
            rho_total = engine.compute_effective_density(rho_baryon_phys, Q_field)
            
            # 5. RENDERING (High-Contrast Cluster Style)
            fig, ax = plt.subplots(1, 2, figsize=(14, 6), facecolor='black')
            
            # Left: Baryons (Hot Map)
            ax[0].imshow(rho_baryon_phys, cmap='hot', origin='lower')
            ax[0].set_title(f"Baryons: {name}", color='white', fontsize=12)
            
            # Right: Total Effective Potential (Coolwarm Halo Style)
            ax[1].imshow(rho_total, cmap='coolwarm', origin='lower')
            ax[1].set_title(r"Mimetic Lensing Map ($\kappa=0.80$)", color='white', fontsize=12)

            for a in ax: a.axis('off')

            plt.tight_layout()
            output_name = f"SENTINEL_REPORT_{name}.png"
            plt.savefig(output_name, facecolor='black')
            plt.close()
            print(f"✅ ANALYSIS LOCKED: {output_name}")
            
        except Exception as e:
            print(f"[RECOVERY_NODE] Analysis failed for {name}: {e}")

    print("\n🏁 BATCH SURVEY COMPLETE. All targets aligned to fits_to_mimetic logic.")

if __name__ == "__main__":
    run_survey()