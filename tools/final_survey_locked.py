# tools/final_survey_locked.py
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.mimetic_engine import MimeticEngine

def run_locked_survey():
    # Final verified paths for your ThinkPad P15
    survey_targets = {
        "Bullet_Cluster": "data/raw_fits/mastDownload/HST/j90701010/j90701010_drz.fits",
        "Abell_370": "data/raw_fits/mastDownload/HST/jabu01030/jabu01030_drz.fits",
        "El_Gordo": "data/raw_fits/mastDownload/HST/jbqz31010/jbqz31010_drz.fits",
        "HUDF_DeepField": "data/raw_fits/mastDownload/HST/j8wc7c010/j8wc7c010_drz.fits"
    }

    engine = MimeticEngine(kappa=0.80)
    print("🛰️ EUCLID SENTINEL: EXECUTING FINAL LOCKED SURVEY\n" + "="*50)

    for name, path in survey_targets.items():
        if not os.path.exists(path): continue
        print(f"🔭 ANALYZING: {name}...")
        
        try:
            with fits.open(path) as hdul:
                # Use HDU[1] for SCI data
                data = hdul[1].data.astype(np.float64)
                data = np.nan_to_num(data)
                
                # Center Crop (1000x1000) to capture core lensing
                h, w = data.shape
                data = data[h//2-500:h//2+500, w//2-500:w//2+500]
                
            # Physics Calculation
            q_field = np.abs(np.gradient(data)[0]) + np.abs(np.gradient(data)[1])
            q_field = (q_field / (np.max(q_field) + 1e-10)) * 10.0
            rho_eff = engine.compute_effective_density(data, q_field)
            
            # Rendering
            fig, ax = plt.subplots(1, 2, figsize=(16, 8), facecolor='black')
            
            # Robust Log-Scale Rendering
            b_disp = np.log10(data - np.min(data) + 1.0)
            m_disp = np.log10(rho_eff - np.min(rho_eff) + 1.0)

            ax[0].imshow(b_disp, cmap='magma', origin='lower')
            ax[0].set_title(f"Hubble: {name}", color='white', fontsize=14)
            
            ax[1].imshow(m_disp, cmap='viridis', origin='lower')
            ax[1].set_title(f"Sentinel Mimetic Potential", color='white', fontsize=14)

            for a in ax: a.axis('off')

            plt.tight_layout()
            plt.savefig(f"SENTINEL_FINAL_{name}.png", facecolor='black')
            plt.close()
            print(f"✅ LOCKED: SENTINEL_FINAL_{name}.png")
            
        except Exception as e:
            print(f"[RECOVERY_NODE] Analysis failed for {name}: {e}")

    print("\n🏁 MISSION COMPLETE. All reports stabilized and locked.")

if __name__ == "__main__":
    run_locked_survey()
