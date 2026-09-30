# Simulation scripts

Simulation scripts for **Hierarchical Conformal Prediction for Clustered Data
with Missing Responses**. Run from the repository root or set `HCP_REPO_DIR`.
Each script has a `repo_dir` setting. Results are saved in `results/`.

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
script defines functions without starting the simulation. For small checks,
set the replicate and worker counts in its runner function.

## 📦 Figure data and plotting

| Result | Data script | Plot script | Data → PDF in `results/` |
|--------|-------------|-------------|--------------------------|
| Figures S.1–S.3 | `reproduce_figureS1S2S3_data.R` | `plot_figureS1S2S3.R` | `figureS1_conditional_data.csv` … `figureS3_conditional_data.csv` → corresponding `figureS1_conditional_plot.pdf` … `figureS3_conditional_plot.pdf` |
| Figures S.4–S.6 | `reproduce_figureS4S5S6_data.R` | `plot_figureS4S5S6.R` | `figureS4_local_data.csv` … `figureS6_local_data.csv` → corresponding `figureS4_local_plot.pdf` … `figureS6_local_plot.pdf` |
| Figure S.7 | `reproduce_figureS7_data.R` | `plot_figureS7.R` | `figureS7_<scenario>_results.csv` → `figureS7.pdf` |
| Figures S.8–S.9 | `reproduce_figureS8S9_data.R` | `plot_figureS8S9.R` | `figureS8_sensitivity_B_data.csv`, `figureS9_sensitivity_S_data.csv` → corresponding PDFs |
| Figure S.10 | `reproduce_figureS10_data.R` | `plot_figureS10.R` | `figureS10_<method>_final_mat_n500.csv` → `figureS10_n500.pdf` |

Run the data script first, then the plotting script. For Figure S.10, the data
script runs HCP, DWR, LC, LMEM, and Oracle on the same 1,000 replicate datasets
and held-out subjects for each DGP, with n=500. Use `--run`, or source the script
and call `run_figureS10_data()`. Then run `plot_figureS10.R`, which reads the five
method CSV files.

## ⚙️ Designs and reproducibility

- Table 1, Tables S.1/S.3, and Figure S.10 compare HCP, DWR, LC, LMEM, and Oracle.
  All use the paper's five DGPs and 1,000 replicates per setting.
- Table S.1 samples cluster sizes independently: 5 or 50, each with probability
  1/2. `R/imbalanced_covariates.R` uses the same five X1 values for both sizes.
  Each appears once for m=5 and ten times for m=50, in random order.
  `R/simulation_dgp.R` generates outcomes and missingness. Balanced simulations
  have five rows per subject.
- Seeds are 1–1000. Table 1 and Tables S.1/S.3 use L'Ecuyer-CMRG for data
  generation; Figure S.10 uses Mersenne-Twister.
- Table 1, Table S.1, and Figure S.10 use weight cap 30. Table S.3 uses the
  sample-size-dependent cap in its script.

The four comparisons use `R/hcp_finite_region.R`. It keeps all components where
p > 0.1, using the fitted density, weights, splits, B/S, and tie/rank/CCT rules.
The finite domain runs from the minimum to the maximum response among non-test
subjects, including simulated responses marked missing. Table S.2, Tables
S.5/S.6, and other figures use the grid evaluators in their scripts.

Simulation LMEM prediction uses 500 draws from a separate L'Ecuyer stream and
restores the caller's RNG state. Oracle is deterministic.

Minor last-digit differences may occur when reproducing the numerical summaries.
