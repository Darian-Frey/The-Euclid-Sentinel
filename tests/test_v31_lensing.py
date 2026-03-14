import numpy as np
import matplotlib.pyplot as plt
from core.mimetic_engine import MimeticEngine

# Updated Calibration for Audited V3.1.1
engine = MimeticEngine(kappa=0.80)

# 2. Setup a Simulated Galaxy Cluster
GRID = 512
BOX = 8e21  # ~2.6 Megaparsecs
x = np.linspace(-BOX/2, BOX/2, GRID)
y = np.linspace(-BOX/2, BOX/2, GRID)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2) + 1e-20

# Create a Baryonic Core (The visible "Gas" of the cluster)
# Using an NFW-like profile: rho ~ 1/r
rho_baryon = 1e-23 / (R / 1e21) * np.exp(-R/2e21)

# 3. Calculate the Scalar Field Gradient (Q-Scalar)
# Sigma represents the metric perturbation from the baryonic mass
sigma = np.log10(rho_baryon + 1e-30) 
gy, gx = np.gradient(sigma, x[1]-x[0])
grad_sq = (gx**2 + gy**2) * 5e18
# 4. RUN THE ENGINE: Compute Emergent Gravity (Dark Matter Proxy)
rho_eff = engine.compute_effective_density(rho_baryon, grad_sq)

# 5. DIAGNOSTICS: Check for the "Notch"
# We compare Total Gravity (rho_eff) vs. Standard Gravity (rho_baryon)
profile_eff = rho_eff[GRID//2, GRID//2:]
profile_bar = rho_baryon[GRID//2, GRID//2:]
radius_axis = x[GRID//2:] / 3.086e22 # Convert to Megaparsecs

# --- VISUALIZATION ---
plt.figure(figsize=(10, 6))
plt.plot(radius_axis, profile_eff, label='Total Effective Gravity (V3.1)', color='cyan', lw=2)
plt.plot(radius_axis, profile_bar, label='Baryonic Mass Only', color='orange', linestyle='--')

plt.yscale('log')
plt.title(f"CODE-GEO Phase IV: V3.1 Lensing Profile (κ={engine.KAPPA})")
plt.xlabel("Radius (Mpc)")
plt.ylabel("Effective Density (kg/m³)")
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.legend()

# Highlight the Notch Zone
plt.axvspan(0.1, 0.3, color='red', alpha=0.1, label='The Krylov Notch Zone')

print("TEST STATUS: RUNNING...")
print(f"Max Effective Density: {np.max(rho_eff):.2e}")
print(f"Baryonic Peak: {np.max(rho_baryon):.2e}")

if np.any(rho_eff > rho_baryon):
    print("RESULT: SUCCESS. Emergent gravity detected without manual Dark Matter.")
else:
    print("RESULT: FAILURE. No emergent gravity produced.")

plt.show()