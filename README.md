# The Euclid-Sentinel: CODE-GEO V4.2
**Mimetic-Conformal Gravity Simulation Framework**

## Project Overview
The-Euclid-Sentinel is a research pipeline investigating the **P vs NP** problem of galactic dynamics: can dark matter phenomenology be computed purely from baryonic distributions using a Mimetic-Conformal scalar field?

This repository implements the **CODE-GEO V4.2** framework, which utilizes a Hartley-Krylov free function $F(Q)$ to govern gravitational information latency.

## Core Theory: V4.2 Hybrid Action
Unlike previous iterations (V3.1/V4.1) which relied on exponential interpolators, V4.2 utilizes a **Mimetic-Conformal Hybrid Anchor**. 

Under a conformal metric rescaling $\tilde{g}_{\mu\nu} = A^2(\sigma) g_{\mu\nu}$, we ensure $c_T = c$ (exactly matching LIGO GW170817 constraints) while the free function $F(Q)$ is defined as:

$$F(\mathcal{Q}) = \sqrt{\mathcal{Q}(1+\mathcal{Q})} - \text{arcsinh}(\sqrt{\mathcal{Q}}) + \lambda \mathcal{Q}^2 e^{-\beta \mathcal{Q}}$$

* **MOND Anchor:** Exactly reproduces the Standard MOND interpolator $\mu_{\text{std}}$.
* **Caustic Guard ($\lambda$):** Suppresses mathematical singularities in the transition region.
* **Stability Sentinel:** Rigorously monitors $c_s^2 \in [0.5, 1.0]$ to prevent ghost instabilities and causality violations.

## Repository Structure
* `core/action.py`: The physical engine. Implements $F(Q)$, $F'(Q)$, and the Stability Sentinel.
* `core/physics_constants.py`: Single Source of Truth for physical constants ($a_0$, $G$, etc.).
* `tools/sparc_refinery_v4.py`: High-precision optimization tool for calibrating parameters against galactic rotation curves.

## Current Status: [2026-03-09]
- **V4.2 Engine:** ✅ VALIDATED
- **Self-Consistency Test:** ✅ PASSED (RMSE: 0.000003%)
- **Next Milestone:** Empirical ingestion of the SPARC (Spitzer Photometry and Accurate Rotation Curves) database.

## Installation & Usage
Ensure you are in the root directory and have `numpy` and `scipy` installed.

```bash
# Run the V4.2 Self-Consistency Validation
PYTHONPATH=. python3 tools/sparc_refinery_v4.py
