# The Euclid Sentinel: Phase IV

"The universe does not hide its secrets; it simply waits for a better lens."

## THE MISSION
The Euclid Sentinel is a specialized diagnostic pipeline designed to ingest lensing maps from the Euclid Space Telescope and compare them against the predictions of the **Mimetic-Conformal V3.1** framework. 

We are searching for the **Krylov Notch**: a specific deficit in shear excess predicted to occur at approximately 8.5 arcminutes in stacked cluster profiles. This is a "smoking gun"—Standard Dark Matter (ΛCDM) cannot produce this feature.

## CORE SPECIFICATIONS (ROM: SCHEMA_V5)
- **Primary Constant:** κ = 0.80 (Hartley-Krylov Damping)
- **Sector:** Conformal-Mimetic (Ghost-free, c_s² → 1)
- **Target:** N ≥ 1,200 Galaxy Clusters (Euclid DR1/DR2)

## REPO STRUCTURE
- `/core`: The F(Q) Lagrangian engine where the "Source of Truth" resides.
- `/tools`: Data loaders for ESA’s Euclid archive.
- `/simulations`: High-resolution Bullet Cluster re-runs.

## STATUS
- [x] Mimetic-Conformal V3.1 Stability Audit
- [x] Caustic Regularization Verification
- [ ] Real-world Data Ingestion (Awaiting Euclid DR1)
