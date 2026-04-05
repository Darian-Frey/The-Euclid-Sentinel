# The Euclid Sentinel
### Mimetic-Conformal Gravity — CODE-GEO V4.2

A computational pipeline for testing the Mimetic-Conformal scalar-tensor gravity framework against galaxy rotation curves and gravitational lensing observations. The framework provides a relativistic, gravitational-wave-compliant completion of Modified Newtonian Dynamics (MOND) with no cold dark matter.

**Status:** SPARC rotation curve analysis complete — preprint in preparation.

---

## The Physics

Standard general relativity requires dark matter to explain the flat rotation curves of galaxies and the observed gravitational lensing. This framework instead modifies gravity at low accelerations through a scalar field $\sigma$ governed by the free function:

$$F(\mathcal{Q}) = \sqrt{\mathcal{Q}(1+\mathcal{Q})} - \text{arcsinh}(\sqrt{\mathcal{Q}})$$

where $\mathcal{Q} = -g^{\mu\nu}\partial_\mu\sigma\partial_\nu\sigma / a_0^2$ is the dimensionless kinetic invariant. This choice anchors $F'(\mathcal{Q})$ exactly to the Standard MOND interpolation function $\mu_\text{std}$, satisfying:

- **Gravitational wave constraint** $c_T = c$ — from the GW170817/GRB 170817A multi-messenger event, satisfied by construction through the conformal coupling
- **Ghost-free propagation** — scalar sound speed $c_s^2 \in [0.50, 1.00]$ verified across the full parameter space
- **Correct MOND limits** — $F'(\mathcal{Q}) \to \sqrt{\mathcal{Q}}$ as $\mathcal{Q} \to 0$ (deep MOND), $F'(\mathcal{Q}) \to 1$ as $\mathcal{Q} \to \infty$ (Newtonian)

---

## Key Results

**SPARC rotation curve fit (171 galaxies, Lelli et al. 2016):**

| Metric | Value |
|:---|:---|
| Median RMSE | 12.68 km/s |
| Mean RMSE | 17.65 km/s |
| Excellent fits (RMSE < 10 km/s) | 61 / 171 |
| Acceptable fits (RMSE < 20 km/s) | 115 / 171 |
| Optimal caustic guard amplitude λ | **0.000** (data prefer the minimal action) |
| Profile curvature at λ=0 | +2.11 (peaked — not insensitive) |
| Joint (λ, a₀) fit — Δa₀ | +4.46% (non-degenerate with λ) |

The caustic guard term $\lambda\mathcal{Q}^2 e^{-\beta\mathcal{Q}}$ is empirically excluded by the SPARC data: the likelihood profile is peaked at $\lambda=0$ with curvature $+2.11$, and a joint free fit over $(\lambda, a_0)$ confirms the exclusion is not degenerate with the MOND scale. The minimal action — zero free parameters beyond $a_0$ — is the best-fit model.

Poor fits (RMSE ≥ 20 km/s, 56 galaxies) concentrate in high-mass, high-surface-brightness, bulge-dominated systems (Spearman $\rho = +0.56$ with $V_\text{flat}$, $p = 10^{-15}$) — the known limitation of the Standard MOND interpolation function $\mu_\text{std}$, not a distinctive failure of the conformal framework.

---

## Repository Structure

```
core/
  action.py              V4.2 free function, derivatives, stability sentinel
  mimetic_engine.py      2D lensing pipeline (thin-disk Poisson kernel)
  physics_constants.py   Single source of truth for all constants
  solvers.py             Stub — AQUAL solver migration target

tools/
  sparc_refinery_v4.py   SPARC rotation curve refinery (main analysis)
  euclid_loader.py       FITS calibration pipeline (flux → Σ_b [kg/m²])
  run_full_survey.py     Cluster lensing survey runner
  rmse_property_analysis.py  RMSE vs galaxy property correlation analysis
  generate_supplementary_table.py  Paper table generator

tests/
  test_gw_compliance.py  GW170817 constraint tests (12 tests)
  test_stability.py      Free function stability suite (13 tests)
  test_v31_lensing.py    V4.2 lensing pipeline tests (13 tests)

docs/
  PREPRINT.md            Paper draft
  ROADMAP.md             Development and publication roadmap

simulations/
  bullet_cluster.py      Synthetic Bullet Cluster validation

data/
  sparc/                 SPARC rotation curve files (gitignored — see below)
  raw_fits/              HST FITS data (gitignored — see below)

survey_outputs/          Generated results (gitignored — regenerate from pipeline)
```

---

## Installation

```bash
pip install -r requirements.txt
```

Requires Python 3.10+. Core dependencies: `numpy`, `scipy`, `matplotlib`, `astropy`.

---

## Data

**SPARC rotation curves** (Lelli et al. 2016):
```bash
# Automatic download (requires internet):
python3 tools/sparc_refinery_v4.py   # downloads on first run

# Manual download:
# https://zenodo.org/records/16284118 → Rotmod_LTG.zip → extract to data/sparc/
```

**HST FITS data** (for cluster lensing — optional):
```bash
python3 tools/fetch_eso_data.py bullet_cluster
python3 tools/fetch_eso_data.py abell370
python3 tools/fetch_eso_data.py elgordo
python3 tools/fetch_eso_data.py hudf
```

---

## Running the Analysis

**SPARC rotation curve fit:**
```bash
# Standard fit (171 galaxies):
python3 tools/sparc_refinery_v4.py

# Full analysis with profile scan and joint fit:
python3 tools/sparc_refinery_v4.py --global-search --profile-scan --joint-fit

# Optimizer self-test (no data required):
python3 tools/sparc_refinery_v4.py --self-test
```

**RMSE property correlation analysis:**
```bash
python3 tools/rmse_property_analysis.py
```

**Supplementary table generation:**
```bash
python3 tools/generate_supplementary_table.py
```

**Cluster lensing survey** (requires HST data):
```bash
python3 tools/run_full_survey.py
```

**Test suite:**
```bash
pytest tests/ -v
# or individually:
python3 tests/test_gw_compliance.py
python3 tests/test_stability.py
python3 tests/test_v31_lensing.py
```

---

## Test Suite

38 tests across three files, all passing:

| File | Tests | Coverage |
|:---|:---:|:---|
| `test_gw_compliance.py` | 12 | c_T=c structural proof, c_s² bounds, λ causality threshold |
| `test_stability.py` | 13 | F(Q) limits, MOND interpolation, positivity, guard localisation |
| `test_v31_lensing.py` | 13 | Thin-disk Poisson kernel, Q-field physics, engine interface |

Notable findings from the test suite:
- λ ≥ 0.44 violates $c_s^2 \leq 1$ (causality bound tightens `LAM_BOUNDS` to 0.35)
- $F(\mathcal{Q}) \sim \frac{2}{3}\mathcal{Q}^{3/2}$ as $\mathcal{Q} \to 0$ — the deep-MOND free function limit

---

## Physical Constants

All constants are defined in `core/physics_constants.py` — the single source of truth. Nothing is hard-coded elsewhere.

| Constant | Value | Role |
|:---|:---|:---|
| $a_0$ | $1.21 \times 10^{-10}$ m/s² | MOND acceleration scale |
| $\lambda$ | 0.000 (SPARC optimum) | Caustic guard amplitude |
| $\beta$ | 1.5 | Guard decay rate (fixed) |
| $\Upsilon_\text{disk}$ | 0.50 M☉/L☉ | Stellar M/L (McGaugh & Schombert 2015) |
| $\Upsilon_\text{bul}$ | 0.70 M☉/L☉ | Bulge M/L |

Note: `KAPPA` (the V3.1 "Hartley-Krylov damping constant", κ=0.80) has been removed. It was a parameter of a deprecated free function and has no role in V4.2. Any code importing `KAPPA` will raise `ImportError` intentionally.

---

## What Changed from Previous Versions

The original pipeline (Phase VIII) had several critical issues that have been corrected:

| Issue | Previous behaviour | V4.2 behaviour |
|:---|:---|:---|
| DM ratio | Hardcoded factor of 12.0 | Emerges from $F'(\mathcal{Q})$ |
| Q normalisation | Forced to peak of 2.0 (arbitrary) | $\mathcal{Q} = (g_N/a_0)^2$ from Poisson |
| Poisson kernel | 3D volumetric $k^2$ on 2D map | 2D thin-disk $\lvert\mathbf{k}\rvert$ |
| MOND formula | $g_\text{eff} = g_\text{bar}/F'(\mathcal{Q}(g_\text{bar}))$ (wrong deep-MOND limit) | Correct AQUAL implicit solve |
| Validation | Self-consistency test labelled "VERIFIED" | Separated: `--self-test` flag |

The Phase VIII "Numerical Universality" results (~12× DM ratio across all targets) were artefacts of the hardcoded factor, not physical predictions.

---

## Roadmap

- **Phase 1** (active): SPARC paper — preprint in preparation
- **Phase 2**: Cluster lensing — Bullet Cluster, Abell 370, El Gordo, HUDF
- **Phase 3**: Euclid weak-lensing confrontation (~2027 data release)
- **Phase 4**: Higher-derivative caustic stabilisation theory

See `docs/ROADMAP.md` for detailed task breakdown.

---

## References

- Lelli, McGaugh & Schombert (2016) — SPARC database, AJ 152 157
- McGaugh, Lelli & Schombert (2016) — Radial Acceleration Relation, PRL 117 201101
- McGaugh & Schombert (2015) — 3.6μm M/L ratios, ApJ 802 18
- Abbott et al. (2017) — GW170817, PRL 119 161101
- Chamseddine & Mukhanov (2013) — Mimetic gravity, JHEP 2013 135
- Bekenstein & Milgrom (1984) — AQUAL formulation, ApJ 286 7

---

**Lead developer:** Shane Frey  
**Hardware:** ThinkPad P15 Gen 2i  
**Preprint:** `docs/PREPRINT.md`