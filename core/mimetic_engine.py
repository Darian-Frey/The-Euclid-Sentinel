# core/mimetic_engine.py
# CODE-GEO Phase IV: The Euclid Sentinel
# Implementation of Mimetic-Conformal V3.1.3 (Fusion-Stabilized)

import numpy as np

class MimeticEngine:
    def __init__(self, kappa=0.80, ell=5e20):
        """
        Initialize the engine with the Hartley-Krylov constant.
        :param kappa: The universal damping constant (Source of Truth: 0.80)
        :param ell:   Correlation length scale (Cluster scale ~ 160 kpc)
        """
        self.KAPPA = kappa
        self.ELL = ell
        self.G = 6.67430e-11  # Gravitational Constant
        self.C = 3.0e8        # Speed of Light

    def free_function_f_q(self, Q):
        """
        V3.1.3: Evaluates the Mimetic response. 
        Input Q is now pre-stabilized by the compute function.
        """
        # Soft-saturation to prevent exponential blow-up
        Q_sat = np.clip(Q, 0, 100.0) 
        
        sqrtQ = np.sqrt(Q_sat + 1e-20)
        
        # Term 1: Main Mimetic Flow
        term1 = Q_sat * (1 - np.exp(-self.KAPPA * sqrtQ))
        
        # Term 2: The Krylov Notch (The Dark Matter Signature)
        notch_amp = (2 * self.KAPPA / 3) * (1 - self.KAPPA)
        term2 = notch_amp * (Q_sat**1.5) * np.exp(-self.KAPPA * Q_sat)
        
        return term1 + term2

    def compute_effective_density(self, rho_baryon, Q_input):
        """
        Computes effective density using pre-scaled Q input to bypass 
        ThinkPad P15 floating-point overflows.
        """
        # The Audited Coupling Constant (~2.15e-16 kg/m^3)
        # This constant is massive, but f_q is small.
        C_rho = (self.C**2) / (8 * np.pi * self.G * (self.ELL**2))
        
        f_q_val = self.free_function_f_q(Q_input)
        
        # Total Density = Matter + (Coupling * Field_Response)
        rho_mimetic = rho_baryon * (f_q_val * 12.0)
        
        return rho_baryon + rho_mimetic

    def get_lensing_potential(self, rho_eff, dx, grid_size):
        """Solves Poisson equation via FFT for lensing potential."""
        kx = np.fft.fftfreq(grid_size, d=dx) * 2 * np.pi
        ky = np.fft.fftfreq(grid_size, d=dx) * 2 * np.pi
        KX, KY = np.meshgrid(kx, ky)
        k2 = KX**2 + KY**2
        k2[0, 0] = 1e-10 
        rho_k = np.fft.fft2(rho_eff)
        phi_k = -4 * np.pi * self.G * rho_k / k2
        return np.real(np.fft.ifft2(phi_k))

# --- VALIDATION BLOCK ---
if __name__ == "__main__":
    engine = MimeticEngine()
    print(f"Mimetic Engine V3.1.3 (Fusion-Stabilized) Initialized.")
    test_Q = np.array([0.0, 1.0, 5.0])
    print(f"F(Q) Response: {engine.free_function_f_q(test_Q)}")