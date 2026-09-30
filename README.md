# Reproducibility Materials

This repository contains code, input data, and numerical results for:

**Hierarchical Conformal Prediction for Clustered Data with Missing Responses**

Tables and figures without an S prefix belong to the main paper. Tables S.1–S.6
and Figures S.1–S.11 belong to the Supplementary Materials.

---

## 📁 Repository structure

| Folder | Contents |
|--------|----------|
| [R/](R/README.md) | Method implementations, data preparation, and replicate-level helpers. |
| [simulation/](simulation/README.md) | Clearly named scripts for each simulation table or figure. |
| [real_data/](real_data/README.md) | CD4 and gallstones reproduction and plotting scripts. |
| [results/](results/) | Manuscript numerical tables, figure data, and PDF figures. |
| [data/](data/) | CD4 and gallstones input datasets. |

## ⚙️ Software requirements

- R >= 4.2.
- R packages: **doParallel, doRNG, dplyr, foreach, gbm, ggplot2, ggh4x, grf,
  gridExtra, lme4, MASS, quantreg, quantregForest, randomForest, tidyr, viridis,
  xgboost**. Base/recommended packages include parallel and stats.
- The simulation calculations used R 4.4.3, lme4 1.1-36, Matrix 1.7-1,
  MASS 7.3-65, and quantreg 6.1. The real-data LMEM calculations used R 4.4.2
  with the same lme4 and Matrix versions. Figure S.10 was rendered using
  ggplot2 4.0.2 and ggh4x 0.3.1.9000.

Minor last-digit numerical variation may occur across reruns or software
environments. Seeds and RNG conventions are specified in the scripts.

## 🚀 How to reproduce a result

The supplied results can be viewed directly in `results/`. To run the code:

1. Set the repository location. Each script has a user-adjustable `repo_dir`:

   ```r
   repo_dir <- normalizePath(Sys.getenv("HCP_REPO_DIR", "."), mustWork = TRUE)
   ```

   Run from the repository root, or set `HCP_REPO_DIR` to its absolute path.

2. Find the target table/figure in the tables below and run its script. For example:

   ```sh
   Rscript simulation/reproduce_table1.R --run
   ```

   This writes `results/table1_final_results_missing20.csv`. The table scripts
   and Figure S.10 data script define functions when sourced; `--run` starts
   their full calculation. Full simulations use 1,000 replicates per setting.
   A small check can be run without overwriting supplied results:

   ```r
   source("simulation/reproduce_table1.R")
   run_table1_parallel(num_simulations_run = 2, num_cores = 2,
                       sample_sizes = 100, scenarios = "Homo", output_file = NULL)
   ```

3. For a figure, run the data script and then the plotting script. For example,
   to reproduce Figure S.10 from start to finish:

   ```sh
   Rscript simulation/reproduce_figureS10_data.R --run
   Rscript simulation/plot_figureS10.R
   ```

   To redraw a figure from its supplied CSV files, run only its plot script.

The Figure S.10 data script computes HCP, DWR, LC, LMEM, and Oracle together on
the same simulated dataset and held-out subject in each replicate: five DGPs,
n=500, and seeds 1–1000 per DGP. It writes all five method CSVs read by the plot
script.

## 📊 Main-paper results

| Result | Reproduction script | Output in `results/` |
|--------|---------------------|----------------------|
| Table 1: simulation, approximately 20% missing | [simulation/reproduce_table1.R](simulation/reproduce_table1.R) | `table1_final_results_missing20.csv` |
| Table 2: CD4 | [real_data/reproduce_table2_cd4.R](real_data/reproduce_table2_cd4.R) | `table2_cd4_results.csv` |
| Figure 2: CD4 prediction bands | [real_data/reproduce_figure2_cd4_data.R](real_data/reproduce_figure2_cd4_data.R), then [real_data/plot_figure2_cd4.R](real_data/plot_figure2_cd4.R) | `figure2_cd4_band_data.csv`, `figure2_cd4_band.pdf` |
| Figure 3: CD4 conditional density | [real_data/reproduce_figure3_cd4_density_data.R](real_data/reproduce_figure3_cd4_density_data.R), then [real_data/plot_figure3_cd4_density.R](real_data/plot_figure3_cd4_density.R) | `figure3_cd4_density_first100_data.csv`, `figure3_cd4_density.pdf` |

Figure 1 is the manuscript's method schematic and has no numerical reproduction script.

## 📊 Supplementary tables

| Result | Reproduction script | Output in `results/` |
|--------|---------------------|----------------------|
| Table S.1: independent cluster sizes 5 or 50 | [simulation/reproduce_tableS1_imbalanced.R](simulation/reproduce_tableS1_imbalanced.R) | `tableS1_imbalanced_cluster_sizes.csv` |
| Table S.2: simultaneous HCP prediction | [simulation/reproduce_tableS2.R](simulation/reproduce_tableS2.R) | `tableS2_final_results.csv` |
| Table S.3: simulation, 50% missing | [simulation/reproduce_tableS3.R](simulation/reproduce_tableS3.R) | `tableS3_final_results_missing50.csv` |
| Table S.4: gallstones | [real_data/reproduce_tableS4_gallstones.R](real_data/reproduce_tableS4_gallstones.R) | `tableS4_gallstones_results.csv` |
| Table S.5: density/residual scores, 20% missing | [simulation/reproduce_tableS5S6_score_compare.R](simulation/reproduce_tableS5S6_score_compare.R) | `tableS5_score_compare_missing20.csv` |
| Table S.6: density/residual scores, 50% missing | [simulation/reproduce_tableS5S6_score_compare.R](simulation/reproduce_tableS5S6_score_compare.R) | `tableS6_score_compare_missing50.csv` |

## 📈 Supplementary figures

| Result | Data script | Plot script | Output PDF(s) in `results/` |
|--------|-------------|-------------|-----------------------------|
| Figures S.1–S.3: conditional coverage | [simulation/reproduce_figureS1S2S3_data.R](simulation/reproduce_figureS1S2S3_data.R) | [simulation/plot_figureS1S2S3.R](simulation/plot_figureS1S2S3.R) | `figureS1_conditional_plot.pdf`, `figureS2_conditional_plot.pdf`, `figureS3_conditional_plot.pdf` |
| Figures S.4–S.6: local coverage | [simulation/reproduce_figureS4S5S6_data.R](simulation/reproduce_figureS4S5S6_data.R) | [simulation/plot_figureS4S5S6.R](simulation/plot_figureS4S5S6.R) | `figureS4_local_plot.pdf`, `figureS5_local_plot.pdf`, `figureS6_local_plot.pdf` |
| Figure S.7: increasing sample size | [simulation/reproduce_figureS7_data.R](simulation/reproduce_figureS7_data.R) | [simulation/plot_figureS7.R](simulation/plot_figureS7.R) | `figureS7.pdf` |
| Figures S.8–S.9: B/S sensitivity | [simulation/reproduce_figureS8S9_data.R](simulation/reproduce_figureS8S9_data.R) | [simulation/plot_figureS8S9.R](simulation/plot_figureS8S9.R) | `figureS8_sensitivity_B.pdf`, `figureS9_sensitivity_S.pdf` |
| Figure S.10: method comparison | [simulation/reproduce_figureS10_data.R](simulation/reproduce_figureS10_data.R) | [simulation/plot_figureS10.R](simulation/plot_figureS10.R) | `figureS10_n500.pdf` |
| Figure S.11: gallstones bands | [real_data/reproduce_figureS11_gallstones_data.R](real_data/reproduce_figureS11_gallstones_data.R) | [real_data/plot_figureS11_gallstones.R](real_data/plot_figureS11_gallstones.R) | `figureS11_gallstones_band.pdf` |

## Method implementations and numerical conventions

- **HCP:** [R/hcp_region.R](R/hcp_region.R) fits the score rule.
  [R/hcp_finite_region.R](R/hcp_finite_region.R) represents that rule as all
  accepted components on the specified finite domain for Table 1, Tables S.1/S.3,
  and Figure S.10. The simulation domain is the minimum and maximum response
  among non-test subjects, including simulated responses whose observation
  indicator is zero. Real-data and other supplementary HCP outputs use their
  documented grid evaluators.
- **DWR / LC:** simulations use [R/dwr_region.R](R/dwr_region.R) and
  [R/lc_region.R](R/lc_region.R); real-data analyses use
  [R/dwr_region_realdata.R](R/dwr_region_realdata.R) and
  [R/lc_region_realdata.R](R/lc_region_realdata.R).
- **LMEM:** [R/lmem_region.R](R/lmem_region.R) is the single implementation for
  simulations and real data. Simulation fits have a fixed intercept, four fixed
  slopes, and four independent random slopes. Real-data fits use an independent
  random intercept and time slope. New-subject intervals include fixed-effect
  estimation uncertainty, random-effect variation, and residual error.
- **Oracle:** [R/oracle_region.R](R/oracle_region.R) is the true conditional HPD
  benchmark. It integrates over the new subject's random coefficients and retains
  every component of the density level set.
- HCP and Oracle lengths are **total Lebesgue measure across represented
  components**, rather than convex-hull width. Table coverage and length are
  averaged within each test subject and then across subjects/replicates. Empty
  DWR/LC regions have coverage zero and missing length; length means omit these
  missing values, as specified by the evaluation helpers.

Details of inputs, models, RNG, and helper files are in the folder README files.
