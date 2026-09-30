# Simulation scripts

These scripts reproduce the simulation results for **Hierarchical Conformal
Prediction for Clustered Data with Missing Responses**. Run from the repository
root, or set `HCP_REPO_DIR`; each script exposes a user-adjustable `repo_dir`.
All result files below are in `results/`.

## 📊 Tables

| Result | Script | Output |
|--------|--------|--------|
| Table 1: approximately 20% missing | `reproduce_table1.R` | `table1_final_results_missing20.csv` |
| Table S.1: cluster sizes 5/50 | `reproduce_tableS1_imbalanced.R` | `tableS1_imbalanced_cluster_sizes.csv` |
| Table S.2: simultaneous HCP | `reproduce_tableS2.R` | `tableS2_final_results.csv` |
| Table S.3: 50% missing | `reproduce_tableS3.R` | `tableS3_final_results_missing50.csv` |
| Tables S.5/S.6: density/residual scores | `reproduce_tableS5S6_score_compare.R` | `tableS5_score_compare_missing20.csv`, `tableS6_score_compare_missing50.csv` |

Table S.4 is reproduced by `real_data/reproduce_tableS4_gallstones.R`.
Run a table with `Rscript simulation/<script-name>.R --run`. Sourcing a table
script defines functions without starting its full experiment. Each table has a
named runner with replicate-count and worker-count arguments for small checks.

## 📦 Figure data and plotting

| Result | Data script | Plot script | Data → PDF in `results/` |
|--------|-------------|-------------|--------------------------|
| Figures S.1–S.3 | `reproduce_figureS1S2S3_data.R` | `plot_figureS1S2S3.R` | `figureS1_conditional_data.csv` … `figureS3_conditional_data.csv` → corresponding `figureS1_conditional_plot.pdf` … `figureS3_conditional_plot.pdf` |
| Figures S.4–S.6 | `reproduce_figureS4S5S6_data.R` | `plot_figureS4S5S6.R` | `figureS4_local_data.csv` … `figureS6_local_data.csv` → corresponding `figureS4_local_plot.pdf` … `figureS6_local_plot.pdf` |
| Figure S.7 | `reproduce_figureS7_data.R` | `plot_figureS7.R` | `figureS7_<scenario>_results.csv` → `figureS7.pdf` |
| Figures S.8–S.9 | `reproduce_figureS8S9_data.R` | `plot_figureS8S9.R` | `figureS8_sensitivity_B_data.csv`, `figureS9_sensitivity_S_data.csv` → corresponding PDFs |
| Figure S.10 | `reproduce_figureS10_data.R` | `plot_figureS10.R` | `figureS10_<method>_final_mat_n500.csv` → `figureS10_n500.pdf` |

Run a data script first, then its plot script. Figure S.10 data generation requires
`--run`; it computes HCP, DWR, LC, LMEM, and Oracle together on each of the same
1,000 replicate datasets per DGP and writes all five method CSVs. Alternatively,
source the script and call `run_figureS10_data()`. The plot reads these five
canonical method files and uses 40 groups of 125 pointwise rows, grouping seed
123, the manuscript method order, and a 6.5 × 9.5 inch PDF. Coverage labels run
from 0.50 to 1.00; length panels scale to their data. Missing group means are not
imputed.

## ⚙️ Designs and reproducibility

- Table 1, Tables S.1/S.3, and Figure S.10 compare HCP, DWR, LC, LMEM, and Oracle.
  They use the same five DGPs, with 1,000 replicate seeds per setting.
- Table S.1 samples each cluster size independently from 5/50 with probability
  1/2. Its dedicated `R/imbalanced_covariates.R` helper repeats each of the five
  X1 values equally often and permutes their order. Outcome and missingness
  generation use `R/simulation_dgp.R`; balanced simulations keep five rows per subject.
- These four comparisons use `R/hcp_finite_region.R`, with the same finite
  candidate domain, density fits, weights, splits, B/S, and calibration settings
  specified by their replicate driver. Table S.2, Tables S.5/S.6, and other figures
  use their documented HCP evaluators and tuning settings.
- Table 1 and Tables S.1/S.3 use L'Ecuyer-CMRG for data generation. Figure S.10
  uses Mersenne-Twister. Seeds are 1–1000. LMEM predictive draws use a separate
  L'Ecuyer stream, restored after each call; Oracle is deterministic.
- Table 1/Table S.1/Figure S.10 use weight cap 30. Table S.3 uses its documented
  sample-size-dependent cap. No parameter is chosen from the resulting coverage.
- The named result scripts expose the design explicitly. Shared loading,
  method diagnostics, and parallel iteration reside in `R/`; there is no separate
  configuration framework. Optional audit outputs contain per-method failures,
  warnings, valid counts, fit diagnostics, and pointwise results.

Minor last-digit differences may occur when reproducing the numerical summaries.
