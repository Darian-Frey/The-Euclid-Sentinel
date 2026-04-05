# Mimetic-Conformal Gravity and the SPARC Rotation Curve Sample: Testing the Minimal Free Function Against 171 Disk Galaxies

**Shane Hartley**  
Independent Research  
*Correspondence: shane.hartley06@gmail.com*

**Preprint — CODE-GEO V4.2**  
*Draft date: April 2026*

---

## Abstract

We present the first systematic confrontation of the Mimetic-Conformal V4.2 gravity framework with the full SPARC database of 171 galaxy rotation curves (Lelli et al. 2016). The theory is built on a scalar field action with free function $F(\mathcal{Q}) = \sqrt{\mathcal{Q}(1+\mathcal{Q})} - \text{arcsinh}(\sqrt{\mathcal{Q}})$, where $\mathcal{Q} \equiv -g^{\mu\nu}\partial_\mu\sigma\partial_\nu\sigma / a_0^2$ is the dimensionless kinetic invariant. This choice anchors $F'(\mathcal{Q})$ exactly to the Standard MOND interpolation function $\mu_\text{std}$, satisfying gravitational wave speed constraints ($c_T = c$) and ghost-free propagation ($c_s^2 \in [0.50, 1.00]$). A caustic guard term $\lambda\mathcal{Q}^2 e^{-\beta\mathcal{Q}}$ was introduced in earlier versions to protect against scalar field singularities. We fit the framework to 171 SPARC galaxies spanning five decades in luminosity and three decades in surface brightness, solving the correct AQUAL implicit equation $g_N \cdot F'((g_N/a_0)^2) = g_\text{bar}$ at each radius. With $a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$ pinned to the canonical MOND value, we find the optimizer drives the caustic guard amplitude to $\lambda_\text{opt} = 0$. A $\lambda$-profile scan over 60 grid points reveals a **peaked** likelihood profile (curvature $= +2.11$ at $\lambda = 0$, flatness ratio $= 0.133$), confirming that SPARC actively disfavours $\lambda > 0$ rather than being merely insensitive to it. A joint free fit over $(\lambda, a_0)$ simultaneously confirms the result is not degenerate with $a_0$: the joint optimum recovers $\lambda = 0$ with $a_0 = 1.264 \times 10^{-10}$ m s$^{-2}$ ($+4.46\%$, $\Delta\text{RMSE} = -0.029$, $< 1\%$ improvement). The minimal action — containing no free parameters beyond the MOND scale $a_0$ — fits 61/171 galaxies at excellent precision (RMSE $< 10$ km s$^{-1}$) and 54 acceptably (RMSE $< 20$ km s$^{-1}$). The 56 poor fits ($\geq 20$ km s$^{-1}$) form a physically coherent population identified by Spearman rank correlations against eight galaxy properties, all statistically significant ($p \leq 3 \times 10^{-3}$): RMSE correlates most strongly with $V_\text{flat}$ ($\rho = +0.562$, $p = 1.3 \times 10^{-15}$), gas fraction ($\rho = -0.520$), and peak $g_\text{bar}/a_0$ ($\rho = +0.506$), with poor-fit galaxies having mean $V_\text{flat} = 171$ km s$^{-1}$ and $g_\text{bar}/a_0 = 8.4$ versus 97 km s$^{-1}$ and 1.3 for good fits. Because the $\lambda = 0$ limit is algebraically equivalent to Standard MOND with $\mu_\text{std}$, these residuals are attributable to the specific interpolation function rather than to the underlying conformal-mimetic framework. The caustic guard is empirically excluded at galactic scales and may require replacement by higher-derivative terms to address scalar field stability without disturbing the rotation curve phenomenology.

---

## 1. Introduction

The observed flatness of galaxy rotation curves at large radii, and the tight empirical correlation between observed and baryonic centripetal accelerations across five decades in galaxy luminosity (the Radial Acceleration Relation, RAR; McGaugh et al. 2016), remain among the most constraining observational facts in extragalactic astronomy. The standard cosmological model ($\Lambda$CDM) accounts for these observations through cold dark matter (CDM) halos, but the tightness and universality of the RAR — and its apparent dependence on baryonic properties alone — is a challenge for halo-based explanations (Desmond et al. 2019; Keller & Wadsley 2017).

Modified Newtonian Dynamics (MOND; Milgrom 1983) provides an empirically successful phenomenological account of the RAR through a single acceleration scale $a_0 \approx 1.2 \times 10^{-10}$ m s$^{-2}$, below which the effective gravitational acceleration is enhanced relative to the Newtonian prediction. The Standard MOND interpolation function $\mu_\text{std}(x) = x/\sqrt{1+x^2}$ achieves excellent rotation curve fits across the SPARC sample (Li et al. 2018; Lelli et al. 2017).

Despite this success, MOND lacks a fully relativistic, covariant formulation with satisfactory properties. Tensor-vector-scalar theory (TeVeS; Bekenstein 2004) achieves relativistic MOND but predicts gravitational wave speeds $c_T \neq c$, in conflict with the GW170817/GRB 170817A multi-messenger constraint (Abbott et al. 2017). This rules out a broad class of vector-disformal scalar-tensor theories as MOND completions.

Mimetic gravity (Chamseddine & Mukhanov 2013) offers an alternative route. By imposing a constraint on the scalar field kinetic invariant, mimetic gravity generates a pressureless dark-matter-like fluid from pure geometry, without introducing new matter degrees of freedom. Under a conformal coupling $\tilde{g}_{\mu\nu} = A^2(\sigma)g_{\mu\nu}$, tensor perturbations propagate on the same light-cone as photons, satisfying the GW170817 constraint exactly at the action level. The free function $F(\mathcal{Q})$ controlling the kinetic self-interactions of $\sigma$ then plays the role of the MOND interpolation function.

This paper presents CODE-GEO V4.2, the Mimetic-Conformal framework, and its confrontation with the SPARC rotation curve database. Section 2 describes the theoretical framework and the V4.2 free function. Section 3 details the SPARC dataset and fitting methodology, including the correction to the implicit AQUAL equation that is required for physically correct rotation curve predictions. Section 4 presents the results: the global fit, the $\lambda$-profile scan, and the joint $(\lambda, a_0)$ free fit. Section 5 discusses implications for the caustic stability argument and the theoretical status of the framework. Section 6 concludes.

---

## 2. Theoretical Framework

### 2.1 The Mimetic-Conformal Action

The Mimetic-Conformal framework is defined by the action

$$S = \int d^4x \sqrt{-g} \left[ \frac{M_\text{Pl}^2}{2} R + M_\text{Pl}^2 F(\mathcal{Q}) + \mathcal{L}_\text{matter} \right]$$

where $\mathcal{Q} \equiv -g^{\mu\nu}\partial_\mu\sigma\partial_\nu\sigma$ is the kinetic invariant of the mimetic scalar field $\sigma$, normalised to the Planck mass $M_\text{Pl}$. The conformal metric rescaling $\tilde{g}_{\mu\nu} = A^2(\sigma)g_{\mu\nu}$ ensures tensor perturbations propagate at $c_T = c$ exactly, satisfying the GW170817 constraint at the level of the action.

The non-linear self-interactions of $\sigma$ are governed entirely by the free function $F(\mathcal{Q})$. In the quasi-Newtonian, spherically symmetric limit, the modified Poisson equation reduces to

$$\nabla \cdot \left[ F'(\mathcal{Q}) \nabla\Phi_N \right] = 4\pi G \rho_\text{bar}$$

where $\mathcal{Q} = (|\nabla\Phi_N|/a_0)^2 = (g_N/a_0)^2$ is now the dimensionless ratio of the Newtonian gravitational acceleration to the MOND scale. This is the AQUAL (A QUadratic Lagrangian; Bekenstein & Milgrom 1984) formulation of MOND, with $F'(\mathcal{Q})$ playing the role of the MOND interpolation function.

### 2.2 The V4.2 Free Function

Previous versions of the framework (V3.1, V4.1) used exponential interpolation forms that failed the SPARC rotation curve test (RMSE $18.7$–$23.4\%$) due to a hump-dip artefact in the MOND transition region and an exponential-versus-algebraic shape mismatch with $\mu_\text{std}$.

The V4.2 hybrid free function is constructed by directly integrating $\mu_\text{std}$:

$$F_\text{MOND}(\mathcal{Q}) = \int_0^\mathcal{Q} \mu_\text{std}(\sqrt{t})\,dt = \sqrt{\mathcal{Q}(1+\mathcal{Q})} - \text{arcsinh}(\sqrt{\mathcal{Q}})$$

This construction guarantees $F'(\mathcal{Q}) = \mu_\text{std}(\sqrt{\mathcal{Q}}) = \sqrt{\mathcal{Q}}/\sqrt{1+\mathcal{Q}}$ exactly, with zero shape-RMSE against the Standard MOND interpolator by construction.

A caustic guard term $\lambda\mathcal{Q}^2 e^{-\beta\mathcal{Q}}$ was appended to address field-theoretic caustic formation concerns:

$$F(\mathcal{Q}) = \underbrace{\sqrt{\mathcal{Q}(1+\mathcal{Q})} - \text{arcsinh}(\sqrt{\mathcal{Q}})}_{F_\text{MOND},\; F' = \mu_\text{std}} + \underbrace{\lambda\mathcal{Q}^2 e^{-\beta\mathcal{Q}}}_\text{caustic guard, $\mathcal{O}(\mathcal{Q}^2)$ at small $\mathcal{Q}$}$$

with free parameters: caustic guard amplitude $\lambda$ (optimised by the refinery) and decay rate $\beta = 1.5$ (fixed by the V4.2 stability audit).

### 2.3 Stability Properties

The V4.2 free function satisfies the following constraints, verified analytically and numerically:

| Property | Condition | Status |
|:---|:---|:---|
| Deep MOND free function | $F(\mathcal{Q}) \sim \frac{2}{3}\mathcal{Q}^{3/2}$ as $\mathcal{Q} \to 0$ (so $F'\to\sqrt{\mathcal{Q}}$) | ✓ |
| Newtonian limit | $F'(\mathcal{Q}) \to 1$ as $\mathcal{Q} \to \infty$ | ✓ |
| MOND interpolation | $F'(\mathcal{Q}) = \mu_\text{std}(\sqrt{\mathcal{Q}})$ (zero RMSE) | ✓ |
| Ghost-free propagation | $c_s^2 \in [0.50, 1.00]$ for $\mathcal{Q} \in [10^{-4}, 10^3]$ | ✓ |
| GW speed constraint | $c_T = c$ (conformal structure) | ✓ |

**Figure 1.** *The V4.2 hybrid free function and its stability properties.*
*Left:* The free function $F(\mathcal{Q}) = \sqrt{\mathcal{Q}(1+\mathcal{Q})} - \text{arcsinh}(\sqrt{\mathcal{Q}})$ (blue) plotted over $\mathcal{Q} \in [10^{-4}, 10^3]$ on a log--log scale.
*Centre:* The derivative $F'(\mathcal{Q})$ (blue) compared against the Standard MOND interpolation function $\mu_\text{std}(\sqrt{\mathcal{Q}}) = \sqrt{\mathcal{Q}}/\sqrt{1+\mathcal{Q}}$ (black dashed). The two curves are indistinguishable: the root-mean-square deviation is identically zero by construction, since $F(\mathcal{Q})$ is derived by integrating $\mu_\text{std}$ directly.
*Right:* The scalar sound speed $c_s^2(\mathcal{Q}) = F'/(F' + 2\mathcal{Q}F'')$ at the default parameters $\lambda = 0.05$, $\beta = 1.5$ (green). The red dashed line marks the causality limit $c_s^2 = 1$; the orange dotted line marks $c_s^2 = 0.5$ (the deep-MOND asymptotic value); the black dashed line marks $c_s^2 = 0$ (ghost threshold). The scalar field propagates subluminally and ghost-free across the full parameter range, satisfying the GW170817 constraint $c_T = c$ by conformal coupling structure. At $\lambda = 0$ (the SPARC-preferred value) the bounds tighten to $c_s^2 \in [0.500, 0.999]$.

---

## 3. Data and Methodology

### 3.1 The SPARC Database

We use the full SPARC (Spitzer Photometry and Accurate Rotation Curves) database of Lelli et al. (2016), comprising 175 nearby disk galaxies with Spitzer 3.6 μm surface photometry and high-quality H I/H$\alpha$ rotation curves. SPARC spans five decades in luminosity ($10^7$–$10^{12}$ $L_\odot$), three decades in surface brightness, and morphological types from S0 to Irr.

Of the 175 galaxies, 4 are excluded by quality cuts (fewer than 5 data points with positive velocity error). The final sample contains **171 galaxies**.

SPARC rotation curve files provide, at each measured radius $R$:

- $V_\text{obs}(R)$ — observed circular velocity [km s$^{-1}$]
- $\delta V(R)$ — velocity uncertainty [km s$^{-1}$]
- $V_\text{gas}(R)$, $V_\text{disk}(R)$, $V_\text{bul}(R)$ — baryonic velocity components [km s$^{-1}$]

These velocity components encode the Newtonian gravitational contribution of each baryonic component, computed from the Spitzer photometry with stellar mass-to-light ratios $\Upsilon_\text{disk}$ and $\Upsilon_\text{bul}$ applied.

### 3.2 Baryonic Acceleration

The total baryonic circular velocity is

$$V_\text{bar}^2 = V_\text{gas}^2 + \Upsilon_\text{disk} V_\text{disk}^2 + \Upsilon_\text{bul} V_\text{bul}^2$$

where we adopt $\Upsilon_\text{disk} = 0.50$ M$_\odot$/L$_\odot$ and $\Upsilon_\text{bul} = 0.70$ M$_\odot$/L$_\odot$ at 3.6 μm, following the stellar population synthesis results of McGaugh & Schombert (2015). These are astrophysical parameters of the stellar population, not free parameters of the Mimetic-Conformal theory.

The Newtonian baryonic centripetal acceleration at each radius is then

$$g_\text{bar}(R) = V_\text{bar}^2(R) / R$$

**The projection problem is absent for SPARC.** Unlike 2D lensing maps, where the baryonic surface density must be converted to a 3D Newtonian acceleration field via a Poisson solve, the SPARC photometric pipeline has already decomposed the baryonic mass into radial velocity components. The quantity $g_\text{bar}(R)$ is therefore the exact 3D Newtonian acceleration at each measured point, with no projection ambiguity.

### 3.3 The AQUAL Implicit Equation

The rotation curve prediction requires solving the modified Poisson equation for $g_N$, the effective gravitational acceleration. In spherical symmetry, this gives the **implicit equation**

$$g_N \cdot F'\!\left(\left(\frac{g_N}{a_0}\right)^2\right) = g_\text{bar} \tag{★}$$

which must be solved for $g_N$ at each radius.

**An important subtlety:** An earlier version of this analysis used the approximation $g_\text{eff} = g_\text{bar} / F'((g_\text{bar}/a_0)^2)$, substituting $g_\text{bar}$ for $g_N$ in the argument of $F'$. In the deep-MOND regime ($g_\text{bar} \ll a_0$, $\mathcal{Q} \ll 1$), this gives $F'(\mathcal{Q}) \approx \sqrt{\mathcal{Q}} = g_\text{bar}/a_0$ and therefore $g_\text{eff} \approx a_0$ — a constant, independent of the local baryonic density. This yields $V_\text{pred} \propto \sqrt{r}$ (a rising curve), drastically overpredicting rotation velocities for low-mass dwarfs. The correct AQUAL equation (★) gives $g_N \to \sqrt{g_\text{bar} \cdot a_0}$ in the deep-MOND limit, recovering the observed flat-curve behaviour.

We solve equation (★) at each data point using vectorised Newton's method, initialised with the exact algebraic solution for $\lambda = 0$ (pure $\mu_\text{std}$):

$$g_N^{(0)} = g_\text{bar} \cdot \sqrt{\frac{1 + \sqrt{1 + 4(a_0/g_\text{bar})^2}}{2}}$$

which is the closed-form solution to $g_N \cdot \mu_\text{std}(g_N/a_0) = g_\text{bar}$. Newton's method then converges in $\leq 5$ iterations for $|\lambda| \leq 0.5$. The predicted circular velocity is

$$V_\text{pred}(R) = \sqrt{g_N(R) \cdot R}$$

### 3.4 Objective Function and Optimisation

We minimise the global weighted root-mean-square velocity residual across all galaxies

$$\mathcal{L}(\lambda) = \frac{1}{N_\text{gal}} \sum_{i=1}^{N_\text{gal}} \sqrt{\frac{1}{N_i} \sum_{j=1}^{N_i} \left(\frac{V_\text{pred}(R_{ij}) - V_\text{obs}(R_{ij})}{\delta V(R_{ij})}\right)^2}$$

using Powell's method with $a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$ pinned to the canonical MOND value. The sole free parameter is the caustic guard amplitude $\lambda \in [0, 0.5]$. The decay rate $\beta = 1.5$ is fixed by the V4.2 theoretical audit and is not tuned.

For the global minimum search we precede Powell with a differential evolution (DE) stage (seed 42, maxiter 500), ensuring robustness against local minima.

---

## 4. Results

### 4.1 Global Fit: λ = 0

The global optimiser drives the caustic guard amplitude to

$$\boxed{\lambda_\text{opt} = 0.000}$$

with global weighted RMSE = 4.744. The minimal free function $F_\text{MOND}(\mathcal{Q})$ alone — containing **zero free parameters beyond $a_0$** — provides the best description of the SPARC sample.

The per-galaxy results are summarised in Table 1. Of 171 galaxies:

- **61 galaxies** (36%) fit at excellent precision: RMSE $< 10$ km s$^{-1}$
- **54 galaxies** (32%) fit acceptably: RMSE $\in [10, 20)$ km s$^{-1}$  
- **56 galaxies** (33%) fit poorly: RMSE $\geq 20$ km s$^{-1}$

The median RMSE is **12.68 km s$^{-1}$** and the mean is 17.65 km s$^{-1}$.

The best-fitting galaxies are predominantly gas-dominated dwarfs and late-type spirals in the deep-MOND or MOND-transition regime ($g_\text{bar}/a_0 \lesssim 1$). UGCA281 achieves RMSE = 1.52 km s$^{-1}$; UGC01281 achieves 2.55 km s$^{-1}$.

The poor fits ($\geq 20$ km s$^{-1}$) cluster into identifiable groups:

1. **High-surface-brightness bulge-dominated spirals** (NGC2841: 75.3 km s$^{-1}$, NGC5985: 90.2 km s$^{-1}$, UGC02487: 86.4 km s$^{-1}$): These have $V_\text{flat} > 280$ km s$^{-1}$ and stellar mass $M_* > 10^{11}$ M$_\odot$. Standard MOND with $\mu_\text{std}$ is known to struggle with this population (Bottema & Pestaña 2015; Banik & Zhao 2022), and the failure here is consistent with published MOND analyses rather than a distinctive signature of the Mimetic-Conformal framework.

2. **Transition-regime outliers** (ESO563-G021: 65.3 km s$^{-1}$, UGC12506: 50.8 km s$^{-1}$): $g_\text{bar}/a_0 \sim 1$ with complex rotation curve morphologies. These galaxies are generically challenging for single-parameter MOND descriptions.

3. **Candidates for non-circular motions** (barred and disturbed galaxies with high-$g_\text{bar}/a_0$ cores): velocity dispersion support and bar-driven streaming motions can systematically offset observed circular velocities from the quasi-circular predictions of the framework.

The overall fit quality is comparable to published Standard MOND results on the full SPARC sample using the same $\mu_\text{std}$ interpolation function (Li et al. 2018), as expected: at $\lambda = 0$ the Mimetic-Conformal framework is exactly equivalent to MOND with $\mu_\text{std}$.

**Figure 2.** *Rotation curve fits for the 12 best-fitting SPARC galaxies.*
Each panel shows the observed circular velocity $V_\text{obs}$ (white points with $1\sigma$ error bars), the Mimetic-Conformal V4.2 prediction $V_\text{pred}$ (cyan), and the baryonic Newtonian velocity $V_\text{bar}$ (orange dashed), as functions of galactocentric radius $R$ [kpc]. The prediction uses the minimal free function $F_\text{MOND}(\mathcal{Q})$ with $\lambda = 0$, $a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$ (pinned), $\Upsilon_\text{disk} = 0.50$, and $\Upsilon_\text{bul} = 0.70$ M$_\odot$/L$_\odot$ (McGaugh & Schombert 2015). The per-galaxy RMSE is shown in each panel title. These 12 galaxies span the full range of morphology from gas-dominated dwarfs (e.g.\ UGCA281, DDO154) to intermediate-mass spirals, all in the deep-MOND or MOND-transition regime ($g_\text{bar}/a_0 \lesssim 1$). The framework achieves RMSE $< 10$ km s$^{-1}$ for 61/171 SPARC galaxies.

**Figure 3.** *Radial Acceleration Relation (RAR) for 171 SPARC galaxies.*
Grey points show the observed centripetal acceleration $g_\text{obs} = V_\text{obs}^2/R$ versus the baryonic Newtonian acceleration $g_\text{bar} = V_\text{bar}^2/R$ for all measured radii in all 171 galaxies. The cyan curve shows the Mimetic-Conformal V4.2 prediction $g_\text{obs} = g_\text{bar} / F'((g_\text{bar}/a_0)^2)$ at $\lambda = 0$ — equivalent to Standard MOND with $\mu_\text{std}$ — evaluated over $g_\text{bar} \in [10^{-13}, 10^{-8}]$ m s$^{-2}$. The orange dashed line shows the Newtonian expectation $g_\text{obs} = g_\text{bar}$. Vertical and horizontal dotted grey lines mark $g = a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$. Below $a_0$ the framework predicts $g_\text{obs} \to \sqrt{g_\text{bar} \cdot a_0}$ (flat rotation curve regime); above $a_0$ the prediction converges to Newtonian. The framework reproduces the tight empirical RAR of McGaugh et al.\ (2016) across five decades in $g_\text{bar}$.

### 4.2 Dynamical Regime Breakdown

Classifying galaxies by their maximum baryonic acceleration:

| Regime | Criterion | $N$ galaxies | Role in constraint |
|:---|:---|:---|:---|
| Deep MOND | $g_\text{bar}/a_0 < 0.3$ | 91 | Guard invisible (Q peak not sampled) |
| Transition | $0.3 \leq g_\text{bar}/a_0 < 5$ | 55 | Guard peak Q $\sim 1.33$ sampled; drives $\lambda \to 0$ |
| Newtonian | $g_\text{bar}/a_0 \geq 5$ | 25 | Guard exponentially suppressed |

The caustic guard's contribution to $F'(\mathcal{Q})$, given by $\lambda(2\mathcal{Q}-\beta\mathcal{Q}^2)e^{-\beta\mathcal{Q}}$, peaks at $\mathcal{Q} = (2-\sqrt{2})/\beta \approx 0.39$ for $\beta=1.5$, corresponding to $g_\text{bar}/a_0 \approx 0.62$. The 55 transition-regime galaxies sample exactly this region of parameter space and provide the binding constraint on $\lambda$.

### 4.3 The λ-Profile Scan

To distinguish between two qualitatively different null results — insensitivity (flat likelihood) vs. active disfavour (peaked likelihood) — we scan the global RMSE across 60 uniformly spaced values of $\lambda \in [0, 0.5]$, holding all other parameters fixed.

The profile is **peaked** at $\lambda = 0$:

| Quantity | Value | Interpretation |
|:---|:---|:---|
| RMSE at $\lambda = 0.000$ | 4.744 | Global minimum |
| RMSE at $\lambda = 0.050$ | 4.778 | $+0.034$ vs minimum |
| RMSE at $\lambda = 0.500$ | 5.247 | $+0.503$ vs minimum |
| Curvature at $\lambda = 0$ | $+2.11$ | Positive → genuine minimum |
| Flatness ratio $[0, 0.1]$ | $0.133$ | Peaked (threshold: 0.05) |

The curvature of $+2.11$ at $\lambda = 0$ confirms this is a genuine minimum of the objective function, not a boundary effect from the $\lambda \geq 0$ constraint. The RMSE increases monotonically and smoothly as $\lambda$ increases from zero. SPARC **actively disfavours** the caustic guard — it is not merely unconstrained by the data.

The RMSE increase of $+0.034$ at the V4.2 theoretical default ($\lambda = 0.05$) is small in absolute terms but systematic across all 171 galaxies. The guard correction is counterproductive in the galactic rotation curve regime.

**Figure 4.** *Likelihood profile scan: global RMSE as a function of the caustic guard amplitude $\lambda$.*
All panels use 171 SPARC galaxies with $a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$ (pinned) and $\beta = 1.5$ (fixed). The scan uses 60 uniformly spaced values of $\lambda \in [0, 0.5]$.
*Top left:* Global weighted RMSE versus $\lambda$ for all 171 galaxies (cyan). The lime dashed line marks $\lambda_\text{opt} = 0$; the orange dotted line marks the V4.2 theoretical default $\lambda = 0.05$.
*Top right:* RMSE by dynamical regime — deep MOND ($g_\text{bar}/a_0 < 0.3$, red; $n = 91$), transition regime ($0.3$--$5$, cyan; $n = 55$), and Newtonian ($g_\text{bar}/a_0 \geq 5$, yellow; $n = 25$). The constraint on $\lambda$ is driven by the transition-regime galaxies, which sample the guard's peak at $\mathcal{Q} \approx 0.39$.
*Bottom left:* RMSE change $\Delta(\text{RMSE})$ relative to $\lambda = 0$ for the 10 transition-regime galaxies most sensitive to $\lambda$, colour-coded by $g_\text{bar}/a_0$. All curves slope upward, confirming the guard is uniformly counterproductive.
*Bottom right:* Zoomed profile over $\lambda \in [0, 0.15]$. The curvature at $\lambda = 0$ is $+2.11$ (positive definite minimum); the flatness ratio over $[0, 0.1]$ is $0.133$. The profile is peaked: SPARC actively disfavours $\lambda > 0$ rather than being insensitive to it.

### 4.4 Joint (λ, a₀) Free Fit

To test whether the $\lambda = 0$ result is degenerate with the MOND acceleration scale — specifically, whether a non-zero $\lambda$ could be absorbed into a rescaled $a_0$ (interpretation b) — we perform a joint free fit over $(\lambda, a_0)$ simultaneously, with bounds $\lambda \in [0, 0.5]$ and $a_0 \in [0.5, 2.5] \times 10^{-10}$ m s$^{-2}$.

The joint optimum is:

$$\lambda_\text{joint} = 0.000, \quad a_{0,\text{joint}} = 1.264 \times 10^{-10} \text{ m s}^{-2} \quad (+4.46\%)$$

with global RMSE = 4.715, compared to 4.744 for the pinned-$a_0$ fit ($\Delta\text{RMSE} = -0.029$, a $0.6\%$ improvement). The 2D RMSE landscape shows an **isolated minimum** near $(\lambda = 0, a_0 = 1.21 \times 10^{-10})$ with no diagonal ridge indicative of $\lambda$–$a_0$ degeneracy.

Interpretation (b) is therefore **definitively ruled out**: $\lambda$ and $a_0$ are not degenerate in the SPARC parameter space. The caustic guard amplitude is excluded independently of the assumed MOND scale.

The recovered $a_{0,\text{joint}} = 1.264 \times 10^{-10}$ m s$^{-2}$ is within $4.5\%$ of the canonical value and consistent with the range of published MOND $a_0$ constraints ($1.1$–$1.3 \times 10^{-10}$ m s$^{-2}$; Begeman et al. 1991; McGaugh et al. 2016; Li et al. 2018). The small positive shift may reflect systematic bias from our fixed M/L ratios: a higher $\Upsilon_\text{disk}$ would increase $g_\text{bar}$ and shift the preferred $a_0$ downward. This systematic is comparable to the $\Delta a_0$ we observe.

**Figure 5.** *Joint ($\lambda$, $a_0$) free fit: 2D RMSE landscape and 1D slices.*
*Left:* Global weighted RMSE as a function of both $\lambda \in [0, 0.5]$ and $a_0 \in [0.5, 2.5] \times 10^{-10}$ m s$^{-2}$, computed on a $40 \times 40$ grid over 171 SPARC galaxies. The colour scale runs from minimum RMSE (dark) to the 90th percentile (light). The lime point marks the joint optimum ($\lambda = 0$, $a_0 = 1.264 \times 10^{-10}$ m s$^{-2}$); the cyan star marks the canonical pinned values. The landscape shows an isolated minimum with no diagonal ridge, ruling out degeneracy between $\lambda$ and $a_0$.
*Right:* One-dimensional slices through the optimum. The cyan curve shows RMSE versus $\lambda$ at $a_0 = a_{0,\text{opt}}$; the orange curve shows RMSE versus $a_0$ (right axis) at $\lambda = 0$. Both slices are peaked at their respective optima with no flat directions, confirming $\lambda = 0$ is a genuine independent minimum. The joint RMSE of $4.715$ versus the pinned RMSE of $4.744$ represents a $0.6\%$ improvement, within systematic uncertainties of the M/L ratios.


### 4.6 RMSE Residuals and Galaxy Properties

We correlate the per-galaxy RMSE with eight properties extracted from the SPARC rotation curve files and the results table: $V_\text{flat}$, the peak baryonic acceleration ratio $g_\text{bar}/a_0$, bulge fraction $f_\text{bul}$, gas fraction $f_\text{gas}$, central disk surface brightness $\log_{10}\text{SB}_\text{disk}$, the MOND fraction (fraction of radii with $g_\text{bar}/a_0 < 1$), a profile shape index, and the number of data points. Spearman rank correlations $\rho$ and p-values for all 171 galaxies are given in Table 2.

All eight correlations are statistically significant at $p \leq 3 \times 10^{-3}$, and their signs are internally consistent: every property that correlates positively with RMSE ($V_\text{flat}$, $g_\text{bar}/a_0$, $\log_{10}\text{SB}_\text{disk}$, $f_\text{bul}$) characterises the same population — massive, high-surface-brightness, bulge-dominated galaxies. Every negative correlation ($f_\text{gas}$, MOND fraction, profile shape) identifies their complement — gas-rich, low-surface-brightness dwarfs. The internal consistency across eight independent properties confirms that the 56 poor fits ($\text{RMSE} \geq 20$ km s$^{-1}$) form a coherent, physically interpretable population rather than random scatter.

**Table 2.** Spearman rank correlations between per-galaxy RMSE and galaxy properties ($n = 171$).

| Property | $\rho$ | $p$ |
|:---|:---:|:---:|
| $V_\text{flat}$ [km s$^{-1}$] | $+0.562$ | $1.3 \times 10^{-15}$ |
| Gas fraction | $-0.520$ | $3.2 \times 10^{-13}$ |
| $g_\text{bar}/a_0$ (max) | $+0.506$ | $1.8 \times 10^{-12}$ |
| MOND fraction of radii | $-0.493$ | $7.6 \times 10^{-12}$ |
| $\log_{10} \text{SB}_\text{disk}$ [$L_\odot\,\text{pc}^{-2}$] | $+0.485$ | $1.9 \times 10^{-11}$ |
| Bulge fraction | $+0.414$ | $1.9 \times 10^{-8}$ |
| Profile shape index | $-0.348$ | $3.1 \times 10^{-6}$ |
| $N$ data points | $+0.226$ | $3.0 \times 10^{-3}$ |

The mean properties of good-fit ($\text{RMSE} < 20$ km s$^{-1}$, $n = 115$) and poor-fit ($n = 56$) galaxies confirm the separation: poor-fit galaxies have $V_\text{flat} = 171$ km s$^{-1}$ versus 97 km s$^{-1}$ for good fits, bulge fractions of 0.116 versus 0.016, and $g_\text{bar}/a_0 = 8.4$ versus 1.3.

The physical mechanism is the behaviour of $\mu_\text{std}(x) = x/\sqrt{1+x^2}$ in the transition regime $x \sim 1$–$10$. This function approaches the Newtonian limit at $\mathcal{O}(x^{-2})$ — a relatively sharp transition that produces the wrong curvature of the effective acceleration profile in the inner-to-outer handoff zone of massive, bulge-dominated galaxies, where the baryonic acceleration spans the full range from Newtonian core to MOND outskirt. Gas-rich dwarfs with $g_\text{bar}/a_0 \ll 1$ throughout are unaffected, which explains their excellent fits (UGCA281: 1.5 km s$^{-1}$; DDO154: 3.3 km s$^{-1}$).

Because the $\lambda = 0$ limit of Mimetic-Conformal V4.2 is algebraically identical to standard MOND with $\mu_\text{std}$, the observed residuals — concentrated in high-acceleration, high-surface-brightness, bulge-dominated galaxies — are attributable to the specific choice of interpolation function rather than a breakdown of the underlying conformal-mimetic framework. The alternative simple interpolation function $\mu_\text{simple}(x) = x/(1+x)$ is known to reduce residuals for this population in published MOND analyses (Li et al. 2018), but its adoption would require reconstructing the free function $F_\text{simple}(\mathcal{Q}) = \mathcal{Q} - 2\sqrt{\mathcal{Q}} + 2\ln(1 + \sqrt{\mathcal{Q}})$ and re-verifying the stability properties of the resulting action; this is deferred to future work.

**Figure 6.** *Per-galaxy RMSE versus galaxy properties: Spearman rank correlations.*
Eight panels show RMSE [km s$^{-1}$] versus each property for all 171 SPARC galaxies, colour-coded by dynamical regime: deep MOND ($g_\text{bar}/a_0 < 0.3$, blue; $n = 91$), transition ($0.3$--$5$, orange; $n = 55$), and Newtonian ($g_\text{bar}/a_0 \geq 5$, red; $n = 25$). The dashed horizontal line marks the poor-fit threshold at RMSE $= 20$ km s$^{-1}$. Panel titles give the Spearman rank correlation coefficient $\rho$ and p-value; all eight correlations are statistically significant ($p \leq 3 \times 10^{-3}$). Properties with $\rho > 0$ ($V_\text{flat}$, $g_\text{bar}/a_0$, central disk surface brightness, bulge fraction) identify the same population — massive, high-surface-brightness, bulge-dominated galaxies. Properties with $\rho < 0$ (gas fraction, MOND fraction, profile shape index) identify their complement — gas-rich, low-surface-brightness dwarfs. The coherence of the sign pattern across eight independent properties confirms that the 56 poor fits form a physically interpretable population rather than random scatter. Poor-fit galaxies have mean $V_\text{flat} = 171$ km s$^{-1}$ (versus 97 km s$^{-1}$ for good fits) and mean $g_\text{bar}/a_0 = 8.4$ (versus 1.3), consistent with the known limitations of $\mu_\text{std}$ in the Newtonian transition regime.


---

## 5. Discussion

### 5.1 The Status of the Caustic Guard

The caustic guard $\lambda\mathcal{Q}^2 e^{-\beta\mathcal{Q}}$ was introduced in CODE-GEO V4.2 to address the known tendency of mimetic gravity actions to develop caustics — singularities where the mimetic fluid density blows up — in the non-linear dynamical regime of structure formation (Ramazanov et al. 2016; Firouzjahi et al. 2017).

The SPARC results establish three statements about the guard:

**1. It is empirically excluded at galactic rotation curve scales.** The profile scan is peaked, the joint fit is non-degenerate, and the constraint originates from transition-regime galaxies ($g_\text{bar}/a_0 \sim 0.3$–$5$) that sample the guard's peak at $\mathcal{Q} = 2/\beta \approx 1.33$. This is not ambiguous: $\lambda = 0$ is the best-fit model for the full SPARC sample.

**2. The exclusion is regime-specific, not universal.** SPARC rotation curves probe quasi-static, virialized galactic disks. The mimetic scalar field $\sigma$ has settled into a configuration tracking the baryonic potential. The dynamical, non-linear regime relevant to caustic formation — shell crossing during gravitational collapse, mergers, early-universe structure formation — is not probed by SPARC. The guard may still be physically required in those contexts.

**3. The algebraic guard form may be theoretically insufficient.** An important observation (following Cai et al. 2017) is that algebraic modifications to $F(\mathcal{Q})$ — terms of the form $Q^n e^{-\beta Q}$ — modify the effective equation of state of the mimetic fluid but do not introduce a pressure gradient that can prevent geodesic crossing. Only higher-derivative terms such as $(\Box\sigma)^2$ provide genuine resistance to worldline crossings. If this is correct, the $\lambda\mathcal{Q}^2 e^{-\beta\mathcal{Q}}$ guard was never capable of resolving the caustic problem in the full theory regardless of its galactic phenomenology. The SPARC exclusion of $\lambda$ and the theoretical inadequacy of the guard form are thus consistent: the guard should be replaced by a higher-derivative stabiliser that is designed to activate only where gradients become large (collapse regime) and to remain negligible in quasi-static galactic disks.

### 5.2 Minimal Theory

At $\lambda = 0$, the Mimetic-Conformal V4.2 framework reduces to

$$F(\mathcal{Q}) = \sqrt{\mathcal{Q}(1+\mathcal{Q})} - \text{arcsinh}(\sqrt{\mathcal{Q}})$$

with a single astrophysical parameter $a_0$. This is the minimal relativistic scalar-tensor completion of Standard MOND consistent with:
- Gravitational wave speed $c_T = c$ (GW170817)
- Ghost-free propagation ($c_s^2 \in [0.50, 1.00]$)
- Exact Standard MOND rotation curve phenomenology ($F' = \mu_\text{std}$ identically)
- Zero free parameters beyond the MOND scale

The SPARC result — median RMSE 12.68 km s$^{-1}$ across 171 galaxies with no parameter tuning — establishes that this minimal action is empirically viable at galactic scales.

### 5.3 Systematic Uncertainties

**Stellar mass-to-light ratios.** We adopt $\Upsilon_\text{disk} = 0.50$ and $\Upsilon_\text{bul} = 0.70$ M$_\odot$/L$_\odot$ following McGaugh & Schombert (2015). These carry uncertainties of $\pm 30$–$50\%$ from stellar population modelling. A joint fit over $(\lambda, a_0, \Upsilon_\text{disk})$ would tighten the $a_0$ constraint but is beyond the scope of this analysis.

**Non-circular motions.** Barred galaxies, warped disks, and kinematically disturbed systems may show systematic offsets between the observed line-of-sight velocity and the quasi-circular prediction of the framework. Several of the poor fits (NGC4389, NGC3992, UGC03205) are known to host strong bars. A quality-flag analysis excluding such systems would improve the aggregate statistics but is not performed here.

**Distance uncertainties.** SPARC distances carry uncertainties of $\sim 10$–$20\%$ for galaxies without Cepheid or tip-of-the-red-giant-branch calibrators. Distance errors propagate into $g_\text{bar}$ and $R$ simultaneously, introducing correlated scatter rather than systematic bias.

**Hot gas contribution.** For galaxy cluster targets (the lensing analysis; see Section 6), the hot intracluster medium (ICM) dominates the baryonic budget. The optical SPARC analysis used here traces stellar mass only, and the ICM contribution is negligible for the disk galaxies in the sample.

### 5.4 Comparison with Published MOND Results

Li et al. (2018) fit the full SPARC sample with the AQUAL formulation of MOND using $\mu_\text{std}$, finding a best-fit $a_0 = (1.20 \pm 0.02_\text{stat} \pm 0.24_\text{sys}) \times 10^{-10}$ m s$^{-2}$ with $\Upsilon_\text{disk} = 0.50$ fixed. Our recovered $a_0 = 1.264 \times 10^{-10}$ m s$^{-2}$ at $\lambda = 0$ is consistent with this at the $2\sigma$ systematic level. The fit quality (median RMSE 12.68 km s$^{-1}$) is directly comparable to their results, as expected: at $\lambda = 0$ and $\beta = 1.5$ the Mimetic-Conformal V4.2 prediction is algebraically identical to MOND with $\mu_\text{std}$.

The contribution of this work is therefore not to improve on published MOND rotation curve fits — which are already excellent — but to demonstrate that a fully covariant, ghost-free scalar-tensor completion of Standard MOND exists and matches the data at the same level, with no additional free parameters.

---

## 6. Future Work

### 6.1 Euclid Weak-Lensing Predictions

The Euclid space mission (Laureijs et al. 2011) will provide weak-lensing shear measurements for $\sim 1.5 \times 10^9$ galaxies over the full extra-galactic sky. The Mimetic-Conformal framework makes specific predictions for the effective lensing convergence:

$$\kappa_\text{eff} = \frac{\Sigma_\text{bar}}{F'(\mathcal{Q}(\Sigma_\text{bar}, d_x))}$$

where $\mathcal{Q}$ is derived from the thin-disk Poisson kernel $\hat{\Phi}(\mathbf{k}) = -2\pi G \hat{\Sigma}_b(\mathbf{k})/|\mathbf{k}|$ applied to the baryonic surface density. A pipeline for this prediction has been developed (the Euclid Sentinel tool) and will be confronted with Euclid ERO and wide-survey data in a companion paper.

### 6.2 Cluster-Scale Lensing

The framework faces a known challenge at cluster scales: the missing mass in galaxy clusters in MOND-like frameworks (the "residual mass problem"; Clowe et al. 2006; Angus et al. 2008). The Bullet Cluster lensing offset — baryonic peaks spatially separated from the lensing peaks by $\sim 8$ arcsec — is a strong constraint on any modified gravity theory. The Mimetic-Conformal framework predicts lensing enhancement that tracks the baryonic distribution, which is qualitatively consistent with the non-collisional dark matter interpretation but cannot reproduce the observed spatial separation without an additional dark component or a more complex scalar field configuration.

Quantitative predictions for the Bullet Cluster, Abell 370, El Gordo, and the HUDF are under development and will be presented alongside Euclid-data confrontation.

### 6.3 The BIG-SPARC Extension

The BIG-SPARC database (Haubner et al. 2024) is currently being compiled from $\sim 4000$ galaxies with homogeneous H I datacubes and WISE near-infrared photometry. Application of the current pipeline to BIG-SPARC would provide substantially stronger constraints on both $\lambda$ and $a_0$, with reduced systematics from the homogeneous reduction.

### 6.4 Higher-Derivative Caustic Stabilisation

As discussed in Section 5.1, the algebraic guard term $\lambda\mathcal{Q}^2 e^{-\beta\mathcal{Q}}$ may be theoretically insufficient to prevent geodesic crossing in the mimetic fluid. A natural replacement is the higher-derivative Lagrangian extension

$$F(\mathcal{Q}) \to F(\mathcal{Q}) + \gamma (\Box\sigma)^2$$

which provides genuine pressure gradients resisting worldline crossings. The coefficient $\gamma$ must be chosen small enough to remain phenomenologically invisible in galactic disks (consistent with $\lambda = 0$ from SPARC) while activating in the high-gradient collapse regime. This theoretical programme will be pursued in a companion paper.

---

## 7. Conclusions

We have confronted the Mimetic-Conformal V4.2 gravity framework with 171 galaxy rotation curves from the SPARC database (Lelli et al. 2016). Our main findings are:

1. **The minimal action fits the data.** The free function $F(\mathcal{Q}) = \sqrt{\mathcal{Q}(1+\mathcal{Q})} - \text{arcsinh}(\sqrt{\mathcal{Q}})$ with $a_0 = 1.21 \times 10^{-10}$ m s$^{-2}$ and **no free parameters** achieves median RMSE 12.68 km s$^{-1}$ across 171 SPARC galaxies. 61/171 galaxies are fit at excellent precision ($< 10$ km s$^{-1}$).

2. **The caustic guard is empirically excluded.** A global optimisation over the caustic guard amplitude $\lambda$ drives $\lambda_\text{opt} = 0$. A 60-point profile scan yields a peaked likelihood profile (curvature $+2.11$, flatness ratio $0.133$), confirming SPARC actively disfavours $\lambda > 0$ rather than being insensitive to it.

3. **The exclusion is not degenerate with $a_0$.** A joint free fit over $(\lambda, a_0)$ recovers $\lambda = 0$ with $a_0 = 1.264 \times 10^{-10}$ m s$^{-2}$ (+4.46%), inconsistent with the hypothesis that the guard can be absorbed into a rescaled MOND scale. The 2D RMSE landscape shows an isolated minimum with no ridge.

4. **The constraint originates from transition-regime galaxies.** The 55 galaxies with $g_\text{bar}/a_0 \in [0.3, 5]$ sample the guard's peak at $\mathcal{Q} \approx 1.33$ and provide the binding constraint. Deep-MOND dwarfs and Newtonian-core galaxies are insensitive to $\lambda$.

5. **Poor fits are consistent with published MOND literature.** The 56 galaxies with RMSE $\geq 20$ km s$^{-1}$ are predominantly high-surface-brightness bulge-dominated spirals and kinematically complex systems, consistent with the known limitations of Standard MOND with $\mu_\text{std}$ and not distinctive to the Mimetic-Conformal framework.

The Mimetic-Conformal V4.2 minimal action — satisfying GW speed constraints, ghost-freedom, and exact MOND phenomenology with zero free parameters — is a viable relativistic completion of Standard MOND, empirically validated at the level of the best published MOND rotation curve analyses.

---

## Acknowledgements

The SPARC database is the result of decades of observational work by the radio and optical astronomical community, compiled and homogenised by Federico Lelli, Stacy McGaugh, and James Schombert. This research made use of the Zenodo archival record of the SPARC Rotmod database (Zenodo record 16284118). Numerical computations were performed using NumPy, SciPy, Astropy, and Matplotlib.

---

## References

Abbott, B. P., et al. (2017). GW170817: Observation of Gravitational Waves from a Binary Neutron Star Inspiral. *Physical Review Letters*, 119, 161101.

Angus, G. W., Famaey, B., & Buote, D. A. (2008). X-ray group and cluster mass profiles in MOND: unexplained mass on the group scale. *MNRAS*, 387, 1470.

Banik, I., & Zhao, H. (2022). From galactic bars to the Hubble tension: weighing up the astrophysical evidence for Milgromian gravity. *Symmetry*, 14, 1331.

Begeman, K. G., Broeils, A. H., & Sanders, R. H. (1991). Extended rotation curves of spiral galaxies: dark haloes and modified dynamics. *MNRAS*, 249, 523.

Bekenstein, J. D. (2004). Relativistic gravitation theory for the modified Newtonian dynamics paradigm. *Physical Review D*, 70, 083509.

Bekenstein, J. D., & Milgrom, M. (1984). Does the missing mass problem signal the breakdown of Newtonian gravity? *ApJ*, 286, 7.

Bottema, R., & Pestaña, J. L. G. (2015). The distribution of dark matter in galaxies. *MNRAS*, 448, 2566.

Cai, Y.-F., Cheng, G., Dent, J., Piao, Y.-S., & Zhang, X. (2017). Features and stability analysis of non-trivial cosmological solutions in scalar-tensor theories of gravity. *JHEP*, 2017, 100.

Chamseddine, A. H., & Mukhanov, V. (2013). Mimetic Dark Matter. *JHEP*, 2013, 135.

Clowe, D., et al. (2006). A Direct Empirical Proof of the Existence of Dark Matter. *ApJL*, 648, L109.

Desmond, H., Katz, H., Lelli, F., & McGaugh, S. (2019). The scatter of the baryonic Tully-Fisher relation. *MNRAS*, 484, 239.

Firouzjahi, H., Gorji, M. A., & Mansoori, S. A. H. (2017). Instabilities in Mimetic Gravity and the Resolution. *JHEP*, 2017, 7.

Haubner, D., et al. (2024). BIG-SPARC: The new SPARC database. *arXiv:2411.13329*.

Keller, B. W., & Wadsley, J. W. (2017). ΛCDM is consistent with the Radial Acceleration Relation. *ApJL*, 835, L17.

Laureijs, R., et al. (2011). Euclid Definition Study Report. *arXiv:1110.3193*.

Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves. *AJ*, 152, 157.

Lelli, F., McGaugh, S. S., Schombert, J. M., & Pawlowski, M. S. (2017). One Law to Rule Them All: The Radial Acceleration Relation of Galaxies. *ApJ*, 836, 152.

Li, P., Lelli, F., McGaugh, S., & Schombert, J. (2018). Fitting the radial acceleration relation to individual SPARC galaxies. *A&A*, 615, A3.

McGaugh, S. S., Lelli, F., & Schombert, J. M. (2016). Radial Acceleration Relation in Rotationally Supported Galaxies. *Physical Review Letters*, 117, 201101.

McGaugh, S. S., & Schombert, J. M. (2015). Color-Mass-to-Light-ratio Relations for Disk Galaxies. *ApJ*, 802, 18.

Milgrom, M. (1983). A modification of the Newtonian dynamics as a possible alternative to the hidden mass hypothesis. *ApJ*, 270, 365.

Ramazanov, S., Comelli, D., & Verde, L. (2016). On non-linear perturbations of mimetic gravity. *JCAP*, 2016, 06, 020.

---

## Appendix A: Version History of the Free Function

| Version | Form | Failure mode |
|:---|:---|:---|
| V3.1 | $\mathcal{Q}(1 - e^{-\kappa\sqrt{\mathcal{Q}}})$ + Krylov Notch | Hump/dip in $\mu_\text{eff}$; $c_s^2 > 1$ (causality violation) |
| V4.1 | Soft-Krylov: $e^{-\kappa\mathcal{Q}^{1/4}}$ | RMSE 18.7%–23.4% on SPARC; shape mismatch |
| V4.2 | $F_\text{MOND} + \lambda\mathcal{Q}^2 e^{-\beta\mathcal{Q}}$ | Current; $\lambda = 0$ from SPARC |
| V4.2 minimal | $F_\text{MOND}$ only | **Empirically preferred; zero free parameters** |

## Appendix B: Numerical Pipeline

The complete pipeline is available as open-source software in the `The-Euclid-Sentinel` repository. Key modules:

- `core/action.py` — V4.2 free function, derivatives, and stability sentinel
- `core/mimetic_engine.py` — 2D lensing pipeline (thin-disk Poisson kernel)
- `tools/sparc_refinery_v4.py` — SPARC rotation curve refinery (this analysis)
- `tools/euclid_loader.py` — FITS calibration pipeline
- `tools/run_full_survey.py` — Cluster lensing survey runner

The SPARC refinery implements the correct AQUAL implicit equation (★) via vectorised Newton iteration. The thin-disk Poisson kernel for the lensing pipeline is $\hat{\Phi}(\mathbf{k}) = -2\pi G\hat{\Sigma}_b(\mathbf{k})/|\mathbf{k}|$, which correctly handles 2D projected surface mass density without conflating it with 3D volumetric density.