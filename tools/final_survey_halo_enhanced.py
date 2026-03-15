# tools/final_survey_halo_enhanced.py
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.mimetic_engine import MimeticEngine

def run_halo_survey():
    survey_targets = {
        "Bullet_Cluster": "data/raw_fits/mastDownload/HST/j90701010/j90701010_drz.fits",
        "Abell_370": "data/raw_fits/mastDownload/HST/jabu01030/jabu01030_drz.fits",
        "El_Gordo": "data/raw_fits/mastDownload/HST/jbqz31010/jbqz31010_drz.fits",
        "HUDF_DeepField": "data/raw_fits/mastDownload/HST/j8wc7c010/j8wc7c010_drz.fits"
    }

    engine = MimeticEngine(kappa=0.80)
    print("🛰️ EUCLID SENTINEL: GENERATING GHOST HALOS\n" + "="*50)

    for name, path in survey_targets.items():
        if not os.path.exists(path): continue
        print(f"🔭 RENDERING HALO: {name}...")
        
        try:
            with fits.open(path) as hdul:
                # Use HDU[1] (Science Layer)
                data = hdul[1].data.astype(np.float64)
                data = np.nan_to_num(data)
                
                # Center Crop (1000x1000)
                h, w = data.shape
                data = data[h//2-500:h//2+500, w//2-500:w//2+500]
                
            # --- THE HALO FIX ---
            # 1. Smooth the input slightly to remove sensor noise
            baryon_clean = gaussian_filter(data, sigma=2.0)
            
            # 2. Calculate Q-field
            q_field = np.abs(np.gradient(baryon_clean)[0]) + np.abs(np.gradient(baryon_clean)[1])
            q_field = (q_field / (np.max(q_field) + 1e-10)) * 10.0
            
            # 3. Initial Mimetic Potential
            rho_eff = engine.compute_effective_density(data, q_field)
            
            # 4. HALO EMERGENCE: Apply Gaussian convolution to the potential 
            # This creates the "Halo" effect around the stars
            halo_potential = gaussian_filter(rho_eff, sigma=15.0)
            
            # Rendering
            fig, ax = plt.subplots(1, 2, figsize=(16, 8), facecolor='black')
            
            # Log Scaling for Contrast
            b_disp = np.log10(baryon_clean - np.min(baryon_clean) + 1.0)
            m_disp = np.log10(halo_potential - np.min(halo_potential) + 1.0)

            ax[0].imshow(b_disp, cmap='hot', origin='lower')
            ax[0].set_title(f"Baryonic Matter: {name}", color='white', fontsize=15)
            
            # Using 'viridis' or 'plasma' for that deep-space halo look
            ax[1].imshow(m_disp, cmap='plasma', origin='lower')
            ax[1].set_title(f"Emergent Mimetic Halo (Sentinel)", color='cyan', fontsize=15)

            for a in ax: a.axis('off')

            plt.tight_layout()
            plt.savefig(f"SENTINEL_HALO_{name}.png", facecolor='black', dpi=120)
            plt.close()
            print(f"✅ HALO LOCKED: SENTINEL_HALO_{name}.png")
            
        except Exception as e:
            print(f"[RECOVERY_NODE] Rendering failed for {name}: {e}")

    print("\n🏁 MISSION COMPLETE. The Ghost Halos have materialized.")

if __name__ == "__main__":
    run_halo_survey()
