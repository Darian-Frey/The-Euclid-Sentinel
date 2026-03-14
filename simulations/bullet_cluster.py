import numpy as np
import matplotlib.pyplot as plt
from core.mimetic_engine import MimeticEngine

# 1. Setup Simulation (Units in METERS)
GRID = 400
BOX = 1e22 
dx = BOX / GRID 
x = np.linspace(-BOX/2, BOX/2, GRID)
y = np.linspace(-BOX/2, BOX/2, GRID)
X, Y = np.meshgrid(x, y)

engine = MimeticEngine(kappa=0.80)

# 2. Define Baryons (Gas) - Units: kg/m^3
center_a = [0.15 * BOX, 0]
rho_a = 5e-23 * np.exp(-((X-center_a[0])**2 + Y**2) / (0.05 * BOX)**2)
center_b = [-0.1 * BOX, 0]
rho_b = 8e-23 * np.exp(-((X-center_b[0])**2 + Y**2) / (0.1 * BOX)**2)
rho_baryon = rho_a + rho_b

# The Sentinel Baseline: Calibrated to 8.5x DM Ratio
sigma_amplitude = 0.025
sigma_center_a = [0.22 * BOX, 0] 
sigma_center_b = [-0.15 * BOX, 0] 
sigma = sigma_amplitude * (np.exp(-((X-sigma_center_a[0])**2 + Y**2) / (0.08 * BOX)**2) + 
                           np.exp(-((X-sigma_center_b[0])**2 + Y**2) / (0.12 * BOX)**2))

# 4. Gradient with PHYSICAL SPACING (The fix for the 'Unit Bomb')
gy, gx = np.gradient(sigma, dx)
grad_sq = (gx**2 + gy**2)

# 5. Compute Emergent Gravity (Using the internal C_rho coupling)
rho_total = engine.compute_effective_density(rho_baryon, grad_sq)

# --- TERMINAL AUDIT ---
print("--- FINAL AUDITED BULLET CLUSTER DATA ---")
print(f"Baryon Max Density:  {np.max(rho_baryon):.2e} kg/m^3")
print(f"Gravity Max Density: {np.max(rho_total):.2e} kg/m^3")
print(f"Effective DM Ratio:  {np.max(rho_total)/np.max(rho_baryon):.2f}x")
print(f"Q-Scale (Peak):      {np.max((engine.ELL**2) * grad_sq):.2e}")
print("-----------------------------------------")

# 6. Plotting
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1); plt.title("Visible Baryons (Gas)"); plt.imshow(rho_baryon, cmap='hot')
plt.subplot(1, 2, 2); plt.title("Lensing Map (CODE-GEO)"); plt.imshow(rho_total, cmap='coolwarm')
plt.show()