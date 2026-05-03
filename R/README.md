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
| `tuning_helpers.R` | Provide helper functions for selecting and tuning model parameters. |

---

## 🧠 Proposed HCP methods

| File | Purpose |
|------|---------|
| `hcp_region.R` | Implement the main hierarchical conformal prediction (HCP) procedure. |
| `hcp_localized_region.R` | Implement the localized HCP procedure. |
| `build_local_groups.R` | Construct local groups used in the localized HCP procedure. |
| `hcp_score_compare.R` | Implement HCP score-comparison variants used in supplementary simulations. |

---

## 📊 Benchmark methods

| File | Purpose |
|------|---------|
| `dwr_region.R` | Implement the DWR benchmark prediction region used in simulations. |
| `lc_region.R` | Implement the LC benchmark prediction region used in simulations. |
| `lmem_region.R` | Implement the linear mixed-effects model benchmark prediction region. |
| `oracle_region.R` | Implement the oracle benchmark prediction region used in simulations. |
| `dwr_region_legacy.R` | Legacy DWR implementation for real-data reproduction. |
| `lc_region_legacy.R` | Legacy LC implementation for real-data reproduction. |

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
| `run_one_replicate.R` | Run one simulation replicate for the main benchmark comparison (Table 1). |
| `run_one_replicate_conditional.R` | Run one replicate for conditional-coverage simulation studies. |
| `run_one_replicate_local.R` | Run one replicate for localized prediction-region simulations and local-coverage evaluation. |
| `run_one_replicate_pointwise.R` | Run one replicate for pointwise evaluation (Figure S10). |
| `run_one_replicate_score_compare.R` | Run one replicate for score-comparison simulation studies. |
| `run_one_replicate_sensitivity.R` | Run one replicate for sensitivity analyses. |
| `run_one_replicate_simultaneous.R` | Run one replicate for simultaneous prediction-region simulation studies. |

---

## 🧪 Real-data drivers

| File | Purpose |
|------|---------|
| `run_one_leaveout_cd4.R` | Perform leave-one-subject-out analysis for the CD4 dataset. |
| `run_one_leaveout_gallstones.R` | Perform leave-one-subject-out analysis for the gallstones dataset. |
