# R helper functions

This folder contains core helper functions used by the reproducibility scripts in the `simulation/` and `real_data/` folders.

These scripts are not intended to be run directly. Instead, they are sourced by higher-level scripts to perform simulation studies and real-data analyses.

The files are organized into the following categories:

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
| `hcp_finite_region.R` | Event-based finite-domain set representation for Table 1, Tables S.1/S.3, and Figure S.10; retains all accepted components. |
| `hcp_localized_region.R` | Implement the localized HCP procedure. |
| `build_local_groups.R` | Construct local groups used in the localized HCP procedure. |
| `hcp_score_compare.R` | Implement HCP score-comparison variants used in supplementary simulations. |

---

## 📊 Benchmark methods

| File | Purpose |
|------|---------|
| `dwr_region.R` | Implement the DWR benchmark prediction region used in simulations. |
| `lc_region.R` | Implement the LC benchmark prediction region used in simulations. |
| `lmem_region.R` | Canonical LMEM fitting and explicit Monte Carlo new-subject prediction for simulations, CD4, and gallstones. |
| `oracle_region.R` | Canonical true conditional HPD Oracle: distribution, density threshold, all component intervals, and total measure. |
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


## ⚙️ Calling the helpers

Source helpers from the corresponding result script; they do not start an
experiment by themselves. `load_simulation_benchmarks.R` loads the same files
on the master process and PSOCK workers. `run_benchmark_settings.R` performs
parallel replicates and aggregation, while `simulation_method_helpers.R` records
method failures, warnings, singularity and convergence without dropping successful
peer methods. Each result still has its own script in `simulation/` or `real_data/`.

Simulation replicate outputs have ten columns: coverage and length for HCP, DWR,
LC, LMEM, and Oracle. The imbalanced design is selected explicitly by the Table
S.1 driver. All primary metrics are subject-weighted. Oracle and finite-domain
HCP lengths sum all components; empty DWR/LC regions use coverage zero and missing
length, excluded from the length mean.

`lmem_region(setting="simulation")` fits a fixed intercept, four fixed slopes,
and four independent random slopes, without a random intercept. The `cd4` and
`gallstones` settings fit independent random intercept and time slope with their
dataset-specific fixed effects. The same function constructs new-subject
prediction draws from beta, random effects, and residual errors, using 500 draws
for simulation and 10,000 for real data. The caller's RNG state is restored.

`oracle_region()` uses the known conditional distribution after integrating over
subject-specific coefficients. Homo/Heter use analytic Gaussian HPD intervals.
Asym uses deterministic Gaussian–Gamma integration, and Bimo/Bimo-Fix use the
true Gaussian mixtures. Numerical density thresholds retain disjoint components
and enforce mass/boundary accuracy checks.

`hcp_region.R` fits the statistical score rule. Its optional passive observer
allows `hcp_finite_region.R` to retain fitted score quantities without refitting or
consuming random numbers. The latter evaluates strict p > 0.1 within the supplied
finite domain, with the specified tie/rank/CCT rules. Other result scripts use
their documented evaluators. The main DGP, density/propensity estimators, and
DWR/LC method definitions are shared across their callers.
