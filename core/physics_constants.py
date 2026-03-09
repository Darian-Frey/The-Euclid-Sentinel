"""
core/physics_constants.py

Baseline physical constants for the The-Euclid-Sentinel simulation framework.

This module defines the standard physical constants and CODE-GEO V3.1–specific
parameters used throughout the Mimetic-Conformal gravity pipeline. In this
framework, dark matter phenomenology is replaced by an information-latency
scalar field governing gravitational dynamics.

Key Concept
-----------
The theory introduces a universal damping parameter, KAPPA, representing the
"Informational Refresh Efficiency" of the gravitational field. This constant
modulates the response of spacetime curvature to mass-energy distribution
within the Mimetic-Conformal V3.1 framework.

Constants defined here serve as the single authoritative reference for
all simulation modules.

Notes
-----
KAPPA = 0.80 is the universal informational refresh efficiency constant
in the CODE-GEO V3.1 Mimetic-Conformal gravity model.
"""

import math

# ---------------------------------------------------------------------
# Standard Physical Constants
# ---------------------------------------------------------------------

# Gravitational constant (m^3 kg^-1 s^-2)
G: float = 6.67430e-11

# Speed of light in vacuum (m/s)
C: float = 299_792_458.0


# ---------------------------------------------------------------------
# CODE-GEO V3.1 Specific Constants
# ---------------------------------------------------------------------

# MOND acceleration scale (m/s^2)
A0: float = 1.21e-10

# Hartley-Krylov damping constant
# Represents the universal informational refresh efficiency of gravity
#KAPPA: float = 0.80
LAMBDA = 0.05
BETA = 1.5


# MOND characteristic length scale
# ℓ = sqrt(c^2 / a0)
ELL: float = math.sqrt((C ** 2) / A0)


# ---------------------------------------------------------------------
# Astrophysical Mass Units
# ---------------------------------------------------------------------

# Solar mass (kg)
MSUN: float = 1.98847e30


# ---------------------------------------------------------------------
# Astrophysical Distance Units
# ---------------------------------------------------------------------

# Kiloparsec to meters
KPC_TO_M: float = 3.08567758e19


# ---------------------------------------------------------------------
# Velocity Conversion Factors
# ---------------------------------------------------------------------

# Kilometers per second to meters per second
KM_S_TO_M_S: float = 1000.0

# Meters per second to kilometers per second
M_S_TO_KM_S: float = 0.001


# ---------------------------------------------------------------------
# Natural Units (Scalar Field / High Energy Normalization)
# ---------------------------------------------------------------------

# Planck mass (kg) — useful for scalar field normalization
M_PL: float = 4.341e-9
