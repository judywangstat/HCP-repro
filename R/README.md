# R helper functions

Helpers for the scripts in `simulation/` and `real_data/`.
Source them from a result script rather than running them directly.

---

## 🧩 Data preparation helpers

| File | Purpose |
|------|---------|
| `cd4_data_helpers.R` | Load and preprocess the CD4 real-data example. |
| `gallstones_data_helpers.R` | Load and preprocess the gallstones real-data example. |

---

## ⚙️ Data-generating process and tuning

| File | Purpose |
|------|---------|
| `simulation_dgp.R` | Define data-generating mechanisms used in the simulation studies. |
| `imbalanced_covariates.R` | Five common X1 values repeated once or ten times for Table S.1. |
| `tuning_helpers.R` | Provide helper functions for selecting and tuning model parameters. |

---

## 🧠 Proposed HCP methods

| File | Purpose |
|------|---------|
| `hcp_region.R` | Implement the main hierarchical conformal prediction (HCP) procedure. |
| `hcp_finite_region.R` | Event-based finite-domain set representation for Table 1, Tables S.1/S.3, and Figure S.10; keeps all accepted components. |
| `hcp_localized_region.R` | Implement the localized HCP procedure. |
| `build_local_groups.R` | Construct local groups used in the localized HCP procedure. |
| `hcp_score_compare.R` | Implement HCP score-comparison variants used in supplementary simulations. |

---

## 📊 Benchmark methods

| File | Purpose |
|------|---------|
| `dwr_region.R` | Implement the DWR benchmark prediction region used in simulations. |
| `lc_region.R` | Implement the LC benchmark prediction region used in simulations. |
| `lmem_region.R` | LMEM fitting and Monte Carlo new-subject prediction for simulations, CD4, and gallstones. |
| `oracle_region.R` | True conditional HPD Oracle: distribution, density threshold, all component intervals, and total measure. |
| `dwr_region_realdata.R` | DWR prediction for the CD4 and gallstones analyses. |
| `lc_region_realdata.R` | LC prediction for the CD4 and gallstones analyses, including optional time-matched subsampling. |

---

## 🔧 Model fitting and evaluation

| File | Purpose |
|------|---------|
| `fit_propensity_model.R` | Fit the propensity/missingness model. |
| `fit_cond_density_qp.R` | Fit the conditional density model. |
| `evaluation_helpers.R` | Compute evaluation metrics such as coverage and prediction-region length. |

---

## 📊 Simulation drivers

| File | Purpose |
|------|---------|
| `run_one_replicate.R` | Run one simulation replicate for the benchmark comparisons (Table 1 and Tables S.1/S.3). |
| `run_one_replicate_conditional.R` | Run one replicate for conditional-coverage simulation studies. |
| `run_one_replicate_local.R` | Run one replicate for localized prediction-region simulations and local-coverage evaluation. |
| `run_one_replicate_pointwise.R` | Run one replicate for pointwise evaluation (Figure S.10). |
| `run_one_replicate_score_compare.R` | Run one replicate for score-comparison simulation studies. |
| `run_one_replicate_sensitivity.R` | Run one replicate for sensitivity analyses. |
| `run_one_replicate_simultaneous.R` | Run one replicate for simultaneous prediction-region simulation studies. |

---

## 🧪 Real-data drivers

| File | Purpose |
|------|---------|
| `run_one_leaveout_cd4.R` | Perform leave-one-subject-out analysis for the CD4 dataset. |
| `run_one_leaveout_gallstones.R` | Perform leave-one-subject-out analysis for the gallstones dataset. |


## ⚙️ How these helpers are used

Result scripts in `simulation/` and `real_data/` source the helpers they need.
Shared simulation setup and parallel execution are handled by
`load_simulation_benchmarks.R` and `run_benchmark_settings.R`.
`simulation_method_helpers.R` collects method-level diagnostics during runs.

Benchmark replicates return coverage and length for HCP, DWR, LC, LMEM, and
Oracle. Table summaries average within the held-out subject first, then across
replicates. HCP and Oracle lengths sum all components. Empty DWR/LC regions
have coverage 0 and missing length; length means omit these missing values.

`lmem_region()` handles simulation and real-data LMEM fits. Simulation models
use a fixed intercept, four fixed slopes, and four independent random slopes,
without a random intercept. The `cd4` and `gallstones` settings use an independent
random intercept and time slope with dataset-specific fixed effects. New-subject
prediction includes fixed-effect uncertainty, random effects, and residual error.

`oracle_region()` constructs the true conditional HPD set after integrating over
subject-specific coefficients. Homo and Heter have analytic Gaussian HPD
intervals; the other settings use the known conditional density. Disconnected
components are kept.

`hcp_region.R` implements the HCP score and calibration procedure.
`hcp_finite_region.R` converts the fitted rule into the finite-domain sets used
in Table 1, Tables S.1/S.3, and Figure S.10. The DGP, density/propensity
estimators, and DWR/LC functions are shared across simulation scripts.
