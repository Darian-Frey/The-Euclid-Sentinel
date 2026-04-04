# The Euclid Sentinel — Development & Publication Roadmap

**Framework:** Mimetic-Conformal V4.2 (CODE-GEO)  
**Last updated:** April 2026  
**Current phase:** Post-SPARC analysis — preprint in progress

---

## State of the Project

### What is real and working

| Component | Status | Notes |
|:---|:---|:---|
| V4.2 free function (`core/action.py`) | ✅ Complete | MOND anchor + stability sentinel verified |
| Mimetic engine V4.2 (`core/mimetic_engine.py`) | ✅ Complete | Thin-disk Poisson kernel; physical Σ_b units |
| SPARC refinery (`tools/sparc_refinery_v4.py`) | ✅ Complete | Correct AQUAL implicit solver; 171-galaxy fit |
| FITS loader (`tools/euclid_loader.py`) | ✅ Complete | Physical calibration; per-target M/L registry |
| Survey runner (`tools/run_full_survey.py`) | ✅ Complete | Calibrated pipeline; M/L sensitivity sweep |
| Bullet Cluster simulation (`simulations/bullet_cluster.py`) | ✅ Complete | Synthetic validation; correct Σ_b units |
| Survey rendering tools | ✅ Complete | Three-panel + four-panel + five-panel outputs |
| Preprint draft (`docs/PREPRINT.md`) | ✅ Draft | SPARC results complete; lensing TBD |

### What was wrong and is now fixed

| Issue | Old behaviour | Fixed behaviour |
|:---|:---|:---|
| Hardcoded 12.0 DM factor | DM ratio *imposed*, not predicted | Emerges from F'(Q) |
| Wrong Poisson kernel | 3D volumetric `k²` on 2D map | 2D thin-disk `\|k\|` |
| Wrong MOND formula | `g_eff = g_bar/F'(Q(g_bar))` — gives g_eff→a₀ in deep MOND | Correct AQUAL implicit solve |
| V4.1 engine in survey tools | Deprecated Soft-Krylov, never wired to V4.2 | All tools now use `free_function_v42` |
| Self-consistency test labelled as physics | "VALIDATED" from fitting own output | Separated into `--self-test` flag |
| Q normalised to arbitrary peak 2.0 | `Q_field = (grad_sq/max) * 2.0` | `Q = (g_N/a₀)²` from Poisson |
| `KAPPA` commented out in constants | Silent fallback to κ=1.0 | Removed; V4.2 uses λ, β, not κ |

### Key scientific results to date

1. **λ_opt = 0.000** — SPARC drives caustic guard amplitude to zero (171 galaxies)
2. **Peaked profile** — curvature +2.11 at λ=0; SPARC actively disfavours λ > 0
3. **Non-degenerate with a₀** — joint (λ, a₀) fit confirms independence; Δa₀ = +4.46%
4. **Constraint from transition regime** — 55 galaxies with g_bar/a₀ ∈ [0.3, 5] do the work
5. **Median RMSE 12.68 km/s** — comparable to published MOND results at λ = 0

---

## Phase 1 — SPARC Paper Completion

*Target: arXiv submission*  
*Estimated effort: 2–4 weeks*

### 1.1 Remaining analysis tasks

**Joint a₀ sensitivity with M/L** *(medium priority)*  
Run the joint fit with Υ_disk swept over [0.3, 0.5, 0.7, 1.0] and report how a₀_opt shifts. This quantifies whether the +4.46% Δa₀ is a stellar population systematic or a genuine tension. Expected: a₀_opt shifts downward toward canonical as Υ_disk increases. Likely a one-session code addition to `sparc_refinery_v4.py`.

**RMSE correlation with galaxy properties** *(medium priority)*  
The 56 poor fits need systematic diagnosis. Correlate per-galaxy RMSE against: inclination, bar presence (morphological flag), effective surface brightness, bulge fraction ($V_\text{bul}/V_\text{obs}$), and distance method. Identifies whether failures are data quality or genuine framework tension. Expected: poor fits cluster in high-inclination, barred, HSB galaxies — consistent with known MOND limitations. Output: one panel for the paper.

**BIG-SPARC preview** *(low priority, stretch goal)*  
If BIG-SPARC data becomes publicly available before submission, a preliminary fit would substantially strengthen the paper. Currently in preparation (Haubner et al. 2024) — monitor arXiv.

### 1.2 Figures needed for submission

| Figure | Source | Status |
|:---|:---|:---|
| Fig 1: F'(Q) vs μ_std, c_s² panel | `core/action.py` self-test | Generate from existing code |
| Fig 2: Rotation curve grid (best 12 fits) | `SPARC_ROTATION_CURVES_V42.png` | ✅ Generated |
| Fig 3: Radial Acceleration Relation | `SPARC_RAR_V42.png` | ✅ Generated |
| Fig 4: λ likelihood profile (4 panels) | `SPARC_LAMBDA_PROFILE_V42.png` | ✅ Generated |
| Fig 5: Joint (λ, a₀) 2D landscape | `SPARC_JOINT_FIT_V42.png` | ✅ Generated |
| Fig 6: RMSE vs galaxy property (diagnostic) | New — from property correlation analysis | ❌ Not yet |

### 1.3 Paper tasks

- [ ] Add author affiliation and contact information to `PREPRINT.md`
- [ ] Write supplementary Table 1 from `SPARC_RESULTS_V42.csv` (per-galaxy results, all 171)
- [ ] Write figure captions for all 6 figures
- [ ] Add Section 4.5: M/L sensitivity (from task 1.1 above)
- [ ] Add Section 4.6: RMSE–property correlation (from task 1.2 above)
- [ ] Review Section 5.1 against Reply 1 (higher-derivative caustic argument) — cite Cai et al. 2017 and Firouzjahi et al. 2017 in body text, not just references
- [ ] Convert to LaTeX for journal submission (target: *Physical Review D* or *A&A*)
- [ ] Internal review pass — flag any claims that exceed what the data support

### 1.4 Code cleanup before paper release

- [ ] `physics_constants.py` — remove commented-out `KAPPA` line; add missing `KAPPA` note
- [ ] `tests/test_gw_compliance.py` — write GW170817 constraint test (c_T = c from action structure)
- [ ] `tests/test_stability.py` — write ghost/causality check using `check_stability_v42()`
- [ ] `tests/test_v31_lensing.py` — update or remove (currently tests deprecated V3.1)
- [ ] `tools/lensing_audit.py` — currently empty; either implement or remove
- [ ] `tools/lensing_sim_v31.py` — currently empty; rename to `lensing_sim_v42.py` and implement
- [ ] `core/conformal_constraints.py` — currently empty; implement conformal constraint checker or remove
- [ ] `core/solvers.py` — currently empty; consider moving implicit AQUAL solver here from refinery
- [ ] README.md — update Phase VIII claim to reflect actual validated state vs. artefacted state
- [ ] GALLERY.md — update or archive (current images were produced with hardcoded 12.0 factor)
- [ ] MANIFEST.json — update status from `VALIDATED_AND_LOCKED` to accurate description

---

## Phase 2 — Cluster Lensing Analysis

*Target: companion paper or extended preprint*  
*Prerequisite: real HST FITS data downloaded*  
*Estimated effort: 4–8 weeks*

### 2.1 The open Q-derivation question

The hub-and-spoke session (April 2026) reached consensus on the thin-disk procedure:

```
Φ(k) = -2πG Σ_b(k) / |k|   →   g_N = |∇Φ|   →   Q = (g_N/a₀)²
```

But a systematic remains flagged in the codebase: the thin-disk approximation is exact for face-on disk galaxies and is a motivated approximation for clusters, but underestimates line-of-sight mass in triaxial merging systems (Bullet Cluster, El Gordo).

**Before the lensing paper can be written, this systematic must be quantified.** Three approaches under evaluation (hub-and-spoke session):

- **(a) X-ray hydrostatic mass** — use Chandra/XMM temperature profiles to derive g_N independently of the optical Σ_b map. Resolves the projection problem for relaxed clusters; fails for mergers.
- **(b) SZ decrement** — Planck/ACT/SPT y-maps provide an independent pressure proxy. Requires ICM modelling.
- **(c) Weak lensing shear** — total (baryons + mimetic enhancement) g_N upper bound from the shear catalogue itself. Requires separating the baryonic and mimetic contributions, which is what we're trying to predict — circular unless done carefully.

**Recommended approach for the paper:** Use approach (a) for relaxed clusters (Abell 370) and document the thin-disk bias for merging systems as a systematic with a bracketed uncertainty.

### 2.2 Data requirements

| Target | Required data | Source | Status |
|:---|:---|:---|:---|
| Bullet Cluster | HST DRZ FITS | MAST (`j90701010_drz.fits`) | ❌ Download needed |
| Abell 370 | HST DRZ FITS | MAST (`jabu01030_drz.fits`) | ❌ Download needed |
| El Gordo | HST DRZ FITS | MAST (`jbqz31010_drz.fits`) | ❌ Download needed |
| HUDF | HST DRZ FITS | MAST (`j8wc7c010_drz.fits`) | ❌ Download needed |
| Bullet Cluster | Chandra X-ray | CXO archive | ❌ Optional; needed for systematic |
| Abell 370 | XMM-Newton X-ray | ESA archive | ❌ Optional |

Download all four HST targets via `tools/fetch_eso_data.py` or directly from MAST (https://mast.stsci.edu).

### 2.3 M/L calibration validation

The `euclid_loader.py` flux-to-Σ_b pipeline makes assumptions (AB magnitude, filter bandwidth, M/L) that need validation against published baryonic mass maps for the Bullet Cluster (Clowe et al. 2006) and Abell 370 (Hoekstra et al. 2011). Specifically:

- Compare our recovered Σ_b peak against published values. Expected: Bullet Cluster Σ_b peak ≈ 0.3–0.5 g cm⁻² in the stellar component (ICM adds ~6× more).
- The ICM gas is the dominant baryonic component for all cluster targets and is completely absent from the optical pipeline. **A companion X-ray Σ_gas map is required for physically meaningful cluster lensing predictions.** This is currently flagged as a `TODO [CLUSTER-PROJECTION]` in `mimetic_engine.py`.

### 2.4 Lensing prediction validation strategy

The key falsifiable prediction of the framework for clusters: the effective lensing convergence

$$\kappa_\text{eff}(R) = \frac{\Sigma_b(R)}{\Sigma_\text{crit} \cdot F'(\mathcal{Q}(R))}$$

should match the observed convergence map from weak lensing (Σ_obs/Σ_crit). For the Bullet Cluster, the framework predicts gravitational enhancement that tracks the baryonic distribution — which will *not* reproduce the observed spatial offset between the lensing peak and the X-ray gas peak. This offset is the core evidence for dark matter in Clowe et al. (2006).

The lensing paper must honestly confront this prediction: either the framework predicts the offset through a non-trivial scalar field configuration, or it predicts no offset and is in tension with the data. Neither outcome is known until the pipeline runs on real data.

---

## Phase 3 — Euclid Confrontation

*Target: publication timed with Euclid wide-survey data release*  
*Prerequisite: Phases 1 and 2 complete*  
*Estimated effort: 3–6 months*

### 3.1 The Euclid opportunity

Euclid's weak-lensing shear catalogue will provide g_obs(R) for ~10^9 galaxies across the extragalactic sky. The RAR can be tested statistically across the full mass and redshift range. For the Mimetic-Conformal framework, the key prediction is:

$$g_\text{obs} = \nu(g_\text{bar}/a_0) \cdot g_\text{bar}$$

where ν is the QUMOND nu function corresponding to F'(Q) = μ_std. This should produce a tight RAR with scatter at the level of the M/L uncertainty — a specific, falsifiable prediction that Euclid weak lensing can test at statistical power far beyond SPARC.

### 3.2 Redshift evolution

The MOND scale a₀ is a constant in the current framework. If Euclid reveals a redshift-dependent RAR — i.e., the relation between g_obs and g_bar shifts with z — this would either:
- Require a modified theory with evolving a₀(z), or
- Be attributable to selection effects and baryonic evolution

The `euclid_loader.py` infrastructure (with per-target redshift registry) is designed to support this analysis.

### 3.3 Required tools not yet built

- `tools/euclid_shear_loader.py` — parse Euclid shear catalogue format (VIS + NISP)
- `tools/rar_euclid.py` — stack weak lensing signal by baryonic mass bin
- `tools/redshift_rar.py` — compute RAR as function of lens redshift

---

## Phase 4 — Theoretical Programme

*Running in parallel with Phases 1–3*

### 4.1 Higher-derivative caustic stabilisation

**Priority: high** — needed to address Reply 1 critique and complete the theoretical picture.

The SPARC result shows the algebraic caustic guard `λQ²e^{-βQ}` is empirically excluded at galactic scales and may be theoretically insufficient (algebraic terms modify EoS but cannot prevent geodesic crossing). The natural replacement is:

$$S_\text{stab} = \gamma \int d^4x \sqrt{-g}\, (\Box\sigma)^2$$

This introduces genuine pressure gradients that resist worldline crossings in the mimetic fluid. The key design requirement: γ must be chosen such that the modification is:
- Negligible in quasi-static galactic disks (consistent with SPARC λ = 0 result)
- Active in high-gradient collapse regime (prevents caustics in structure formation)

The action-level analysis of this extension is the next theoretical paper.

### 4.2 Open question: the higher-derivative question for the hub-and-spoke

*Recommended next hub-and-spoke question for the three-LLM ensemble:*

> *In mimetic gravity, algebraic modifications to F(Q) — including terms of the form λQ²e^{-βQ} — modify the effective equation of state of the mimetic fluid but do not introduce pressure gradients that prevent geodesic crossing. Only higher-derivative terms like (□σ)² provide genuine caustic resistance. (a) Is this claim correct at the level of the equations of motion? (b) What is the correct coefficient γ such that (□σ)² is phenomenologically invisible in galactic disks (g_N ~ 10⁻¹¹–10⁻⁸ m/s²) but activates in the collapse regime (δ ~ 1, virialization)? (c) Does the addition of (□σ)² break the c_T = c constraint from GW170817, and if so what coupling structure avoids this?*

### 4.3 Gravitational lensing at the action level

The current lensing pipeline uses the quasi-Newtonian AQUAL limit. A proper weak lensing analysis requires the full post-Newtonian deflection angle from the conformal scalar-tensor action, including:
- Gravitomagnetic corrections at first post-Newtonian order
- The relationship between the mimetic scalar field configuration and the lensing potential Ψ vs. the Newtonian potential Φ (which can differ in scalar-tensor theories)

This is a theoretical calculation for the companion paper.

### 4.4 Cluster mass problem

MOND-like frameworks generically predict insufficient mass in galaxy clusters (the "residual mass problem"). The Mimetic-Conformal framework inherits this problem. Candidate resolutions:
- Massive neutrinos (Σmν ~ 2 eV; Angus et al. 2008) as the cluster dark matter
- Non-trivial scalar field configurations in merging systems
- A genuinely different cluster regime requiring additional physics

This should be confronted honestly in the lensing paper, not avoided.

---

## Immediate Next Steps (This Week)

In priority order:

1. **`physics_constants.py` cleanup** — remove commented `KAPPA`, fix missing constants, add version note documenting V4.2 parameter set
2. **Write the three test files** — `test_gw_compliance.py`, `test_stability.py`, `test_v31_lensing.py` (update to V4.2)
3. **Download HST FITS data** — run `tools/fetch_eso_data.py` for all four cluster targets
4. **Hub-and-spoke session** — send the higher-derivative question (Section 4.2 above) to three LLMs
5. **M/L sensitivity run** — `python3 tools/sparc_refinery_v4.py --joint-fit` extended to sweep Υ_disk

---

## Timeline

| Milestone | Target | Blocker |
|:---|:---|:---|
| SPARC paper complete | May 2026 | Figures, table, LaTeX conversion |
| HST data downloaded | April 2026 | Network / MAST access |
| Cluster lensing pipeline tested | June 2026 | HST data + X-ray data |
| Higher-derivative theory drafted | June 2026 | Hub-and-spoke session |
| Lensing paper draft | Aug 2026 | Cluster pipeline + theory |
| Euclid confrontation | 2027 | Euclid wide-survey data |

---

## Notes on Epistemic Honesty

The original Phase VIII "VALIDATED" status was incorrect: all survey results were artefacts of a hardcoded 12.0 amplification factor and an incorrect Q normalisation. The pipeline has been rebuilt from first principles. Current results are physically grounded but should be treated as:

- **SPARC fit**: robust — 171 real galaxies, correct physics, reproducible
- **Cluster lensing**: not yet tested on real data — the pipeline is built and calibrated but FITS data has not been downloaded and run end-to-end
- **Euclid predictions**: forward projections only — no confrontation with Euclid data yet

Any future claim of "validation" must be grounded in real observational data processed through the V4.2 pipeline, not in the pipeline's ability to produce plausible-looking images from hardcoded parameters.