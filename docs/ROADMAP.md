# The Euclid Sentinel — Development & Publication Roadmap

**Framework:** Mimetic-Conformal V4.2 (CODE-GEO)  
**Last updated:** April 2026  
**Current phase:** Phase 1 — SPARC paper final preparation

---

## State of the Project

### What is complete and verified

| Component | Status | Notes |
|:---|:---|:---|
| V4.2 free function (`core/action.py`) | ✅ Complete | MOND anchor + stability sentinel verified |
| Mimetic engine V4.2 (`core/mimetic_engine.py`) | ✅ Complete | Thin-disk Poisson kernel; physical Σ_b units |
| Physics constants (`core/physics_constants.py`) | ✅ Complete | KAPPA removed; C_LIGHT alias; MSUN_PER_PC2_TO_SI added |
| SPARC refinery (`tools/sparc_refinery_v4.py`) | ✅ Complete | AQUAL implicit solver; 171-galaxy fit; profile scan; joint fit |
| FITS loader (`tools/euclid_loader.py`) | ✅ Complete | Physical calibration; per-target M/L registry; provenance dict |
| Survey runner (`tools/run_full_survey.py`) | ✅ Complete | Calibrated pipeline; M/L sensitivity sweep |
| RMSE property analysis (`tools/rmse_property_analysis.py`) | ✅ Complete | 8-property Spearman correlation; all significant |
| Supplementary table generator (`tools/generate_supplementary_table.py`) | ✅ Complete | Markdown + LaTeX + full CSV |
| Test suite — GW compliance | ✅ 12/12 | c_T=c structural; c_s² bounds; λ causality threshold |
| Test suite — Stability | ✅ 13/13 | F(Q) limits; MOND interpolation; guard localisation |
| Test suite — Lensing pipeline | ✅ 13/13 | Thin-disk Poisson; Q-field physics; engine interface |
| Preprint (`docs/PREPRINT.md`) | ✅ Draft | All SPARC sections complete including 4.6; lensing TBD |
| README.md | ✅ Updated | Actual V4.2 results; Phase VIII artefacts documented |
| MANIFEST.json | ✅ Updated | Accurate status; deprecation log |
| requirements.txt | ✅ Populated | All dependencies listed |
| .gitignore | ✅ Written | Data, outputs, logs, pycache excluded |
| Repo cleanup | ✅ Done | Stray files deleted; Phase VIII PNGs archived |

### Key scientific results

1. **λ_opt = 0.000** — SPARC drives caustic guard to zero across 171 galaxies
2. **Peaked profile** — curvature +2.11 at λ=0; data actively disfavour λ > 0
3. **Non-degenerate with a₀** — joint (λ, a₀) fit: λ stays zero, a₀ shifts +4.46%
4. **λ ≥ 0.44 violates causality** — c_s² > 1 discovered during testing; LAM_BOUNDS tightened to 0.35
5. **Median RMSE 12.68 km/s** — 61/171 excellent, 54 acceptable, 56 poor
6. **Poor fits are not random** — V_flat ρ=+0.562 (p=10⁻¹⁵); attributable to μ_std, not the conformal framework
7. **F(Q) ~ (2/3)Q^{3/2} as Q→0** — canonical kinetic limit claim in old docstring was incorrect; corrected

### What was wrong and is now fixed

| Issue | Old behaviour | Fixed behaviour |
|:---|:---|:---|
| Hardcoded 12.0 DM factor | DM ratio imposed | Emerges from F'(Q) |
| Wrong Poisson kernel | 3D volumetric k² on 2D map | 2D thin-disk \|k\| |
| Wrong MOND formula | g_eff = g_bar/F'(Q(g_bar)) — constant in deep MOND | Correct AQUAL implicit solve |
| V4.1 engine in survey tools | Deprecated, never wired to V4.2 | All tools use free_function_v42 |
| Self-consistency test mislabelled | "VALIDATED" from fitting own output | Separated into --self-test flag |
| Q normalised arbitrarily | Q_field = (grad_sq/max) × 2.0 | Q = (g_N/a₀)² from Poisson |
| KAPPA commented out | Silent fallback to κ=1.0 | Removed; ImportError intentional |
| Σ_b units wrong | rho [kg/m³] on projected map | sigma_b [kg/m²] thin-disk |
| Docstring: F(Q)→Q as Q→0 | False claim | F(Q)~(2/3)Q^{3/2} as Q→0 |

---

## Phase 1 — SPARC Paper

*Target: arXiv submission*  
*Status: final preparation — science complete, writing remains*

### 1.1 Remaining before submission

| Task | Priority | Status |
|:---|:---|:---|
| Run supplementary table script | High | ❌ `python3 tools/generate_supplementary_table.py` |
| Generate Fig 1 (F'(Q) vs μ_std) | High | ❌ `python3 core/action.py` → action_self_test_v42.png |
| Write figure captions (6 figures) | High | ❌ ~1–2 hours |
| Update abstract with 4.6 numbers | High | ❌ Written before Section 4.6 existed |
| Add author affiliation | High | ❌ Placeholder in PREPRINT.md |
| M/L sensitivity sweep (Υ_disk variation) | Medium | ❌ Quantifies +4.46% Δa₀ systematic |
| Internal review pass | Medium | ❌ Built section by section; coherence check needed |
| LaTeX conversion | Low | ❌ Mechanical; journal choice affects style |

### 1.2 Figures — status

| Figure | File | Status |
|:---|:---|:---|
| Fig 1: F'(Q) vs μ_std, c_s² | `action_self_test_v42.png` | ❌ Run `python3 core/action.py` |
| Fig 2: Rotation curve grid (best 12) | `SPARC_ROTATION_CURVES_V42.png` | ✅ Generated |
| Fig 3: Radial Acceleration Relation | `SPARC_RAR_V42.png` | ✅ Generated |
| Fig 4: λ likelihood profile (4 panels) | `SPARC_LAMBDA_PROFILE_V42.png` | ✅ Generated |
| Fig 5: Joint (λ, a₀) landscape | `SPARC_JOINT_FIT_V42.png` | ✅ Generated |
| Fig 6: RMSE vs galaxy properties | `RMSE_PROPERTY_ANALYSIS.png` | ✅ Generated |

### 1.3 New findings to incorporate into the paper

**Causality bound on λ:**  
Tests found λ ≥ 0.44 causes c_s² > 1 (superluminal scalar propagation). LAM_BOUNDS updated to 0.35 in the refinery. Mention briefly in Section 2.3 (stability properties).

**μ_simple alternative:**  
Three-LLM synthesis (Q-three session) concluded that switching the free function anchor to μ_simple = x/(1+x) would improve fits for the 56 HSB poor-fit galaxies. The corresponding free function is:
$$F_\text{simple}(\mathcal{Q}) = \mathcal{Q} - 2\sqrt{\mathcal{Q}} + 2\ln(1 + \sqrt{\mathcal{Q}})$$
Mentioned in Section 5.1 as future work. A follow-up SPARC run would produce one publishable figure.

**RMSE-property correlations (Section 4.6 — now written):**  
All 8 correlations significant (p ≤ 3×10⁻³). Top predictors: V_flat ρ=+0.562, gas fraction ρ=−0.520, g_bar/a₀ ρ=+0.506. Poor fits form a coherent, physically interpretable population — not random scatter.

### 1.4 Code cleanup — all complete

- ✅ `physics_constants.py` — KAPPA removed, C_LIGHT alias, new constants
- ✅ `tests/test_gw_compliance.py` — 12 tests, passing
- ✅ `tests/test_stability.py` — 13 tests, passing
- ✅ `tests/test_v31_lensing.py` — updated to V4.2, 13 tests passing
- ✅ `core/conformal_constraints.py` — stub with planned scope
- ✅ `core/solvers.py` — stub documenting AQUAL solver migration target
- ✅ `README.md` — fully rewritten
- ✅ `GALLERY.md` — Phase VIII outputs archived
- ✅ `MANIFEST.json` — updated

---

## Phase 2 — Cluster Lensing Analysis

*Target: companion paper or extended preprint*  
*Status: pipeline built; HST data downloaded; not yet run end-to-end*

### 2.1 Data status

| Target | HST FITS | Status |
|:---|:---|:---|
| Bullet Cluster | `j90701010_drz.fits` | ✅ Downloaded |
| Abell 370 | `jabu01030_drz.fits` | ✅ Downloaded |
| El Gordo | `jbqz31010_drz.fits` | ✅ Downloaded |
| HUDF | `j8wc7c010_drz.fits` | ✅ Downloaded |
| Bullet Cluster | Chandra X-ray | ❌ Needed for ICM systematic |
| Abell 370 | XMM-Newton X-ray | ❌ Optional |

### 2.2 Before the first pipeline run

**ICM gas is missing.** Optical pipeline traces stellar mass only. For all cluster targets the ICM outweighs the stellar component (~6:1 for Bullet Cluster). The lensing predictions will be systematically biased (g_bar underestimated → Q underestimated → DM ratio overestimated).

**Decision required:**
- **(a)** Run with stellar component only, document the systematic floor — valid for a first-pass paper
- **(b)** Obtain Chandra/XMM data and add X-ray Σ_gas map to `euclid_loader.py` — required for a complete result

Recommendation: (a) for the companion paper to establish the baseline.

**M/L validation first:** Before running, validate `euclid_loader.py` outputs against published baryonic mass maps. Run `python3 tools/fits_to_mimetic.py --target bullet_cluster` and check Σ_b peak is in the right range.

### 2.3 The Bullet Cluster tension

The framework predicts lensing enhancement tracking the baryonic distribution. The observed lensing peak is offset ~8 arcsec from the X-ray gas peak (Clowe et al. 2006) — the framework will not reproduce this without a non-trivial scalar field configuration. This must be confronted honestly in the lensing paper, not avoided.

### 2.4 Q-derivation — resolved

Hub-and-spoke session (April 2026) concluded:
- ✅ Thin-disk Poisson kernel Φ(k) = -2πG Σ_b(k)/|k| — correct procedure
- ✅ Input: Σ_b [kg/m²] — not ∇flux, not ∇κ
- ✅ Q computed post-deprojection
- ⚠️ Thin-disk underestimates line-of-sight mass in triaxial merging clusters — documented systematic

---

## Phase 3 — Euclid Confrontation

*Target: publication timed with Euclid wide-survey data release (~2027)*  
*Status: forward planning only*

Euclid will provide weak-lensing shear for ~10⁹ galaxies. The framework predicts a tight RAR with scatter at the M/L level — falsifiable at statistical power far beyond SPARC.

Tools not yet built:
- `tools/euclid_shear_loader.py` — parse VIS + NISP shear catalogue
- `tools/rar_euclid.py` — stack lensing signal by baryonic mass bin
- `tools/redshift_rar.py` — RAR as function of lens redshift (tests a₀(z) = const)

---

## Phase 4 — Theoretical Programme

*Running in parallel with Phases 1–3*

### 4.1 Higher-derivative caustic stabilisation

The algebraic guard `λQ²e^{-βQ}` is empirically excluded by SPARC and theoretically insufficient — algebraic F(Q) terms modify the equation of state but cannot prevent geodesic crossing. The natural replacement:

$$S_\text{stab} = \gamma \int d^4x \sqrt{-g}\, (\Box\sigma)^2$$

γ must be negligible in quasi-static galactic disks (consistent with λ=0) and active in the collapse regime.

**Hub-and-spoke question — not yet sent:**
> *In mimetic gravity, algebraic F(Q) terms cannot prevent geodesic crossing — only (□σ)² provides genuine caustic resistance. (a) Is this correct at the equations-of-motion level? (b) What coefficient γ makes (□σ)² invisible in galactic disks but active in the collapse regime? (c) Does (□σ)² break the c_T=c constraint, and what coupling structure avoids this?*

### 4.2 μ_simple free function

Implementation path:
1. Add `free_function_simple()` and derivatives to `core/action.py`
2. Verify c_s² ∈ [0, 1] for the new functions
3. Add `--anchor simple` flag to `sparc_refinery_v4.py`
4. Run 171-galaxy fit; report RMSE improvement for HSB galaxies
5. Check whether the guard remains excluded under μ_simple

### 4.3 Gravitational lensing at the action level

The lensing pipeline uses the quasi-Newtonian AQUAL limit. A complete weak-lensing analysis requires the full post-Newtonian deflection from the conformal scalar-tensor action — whether Ψ differs from Φ in this theory. Theoretical calculation for the lensing paper.

### 4.4 Cluster residual mass problem

MOND-like frameworks predict insufficient cluster mass without additional dark matter. Candidate resolutions: massive neutrinos (Σmν ~ 2 eV), non-trivial scalar field configurations in mergers, different cluster-regime physics. Confront honestly in the lensing paper.

---

## Timeline — Updated

| Milestone | Target | Status |
|:---|:---|:---|
| Test suite 38/38 passing | April 2026 | ✅ Done |
| Section 4.6 (RMSE-property) written | April 2026 | ✅ Done |
| Repo cleanup and gitignore | April 2026 | ✅ Done |
| README, MANIFEST, ROADMAP updated | April 2026 | ✅ Done |
| action.py docstring fix | April 2026 | ✅ Done |
| Run supplementary table script | April 2026 | ❌ 5 minutes |
| Generate Fig 1 | April 2026 | ❌ 2 minutes |
| Figure captions written | April 2026 | ❌ 1–2 hours |
| Abstract updated | April 2026 | ❌ 30 minutes |
| SPARC paper submitted to arXiv | May 2026 | ❌ LaTeX conversion |
| Cluster pipeline first run | May 2026 | ❌ HST data ready |
| μ_simple SPARC fit | May 2026 | ❌ One session |
| Higher-derivative hub-and-spoke | May 2026 | ❌ Send the question |
| Lensing paper draft | Aug 2026 | ❌ |
| Euclid confrontation | 2027 | ❌ |

---

## Epistemic Honesty

- **SPARC fit**: robust — 171 real galaxies, correct physics, 38 passing tests, reproducible
- **Cluster lensing**: pipeline built and calibrated, HST data downloaded, not yet run end-to-end
- **Euclid predictions**: forward projections only — no data confrontation yet
- **Caustic stability**: empirically excluded at galactic scales; theoretical status (higher-derivative question) not yet resolved

The Phase VIII "Numerical Universality" results were artefacts of a hardcoded amplification factor, not physical predictions. Archived in `archive/phase_08/`.