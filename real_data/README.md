# Real-data analysis scripts

These scripts reproduce the CD4 and gallstones analyses. Input datasets are in
`data/`, helper functions in `R/`, and supplied numerical outputs in `results/`.
Run from the repository root or set `HCP_REPO_DIR`; `repo_dir` remains adjustable
at the top of every script.

## 📊 Tables

| Result | Script | Input | Output in `results/` |
|--------|--------|-------|----------------------|
| Table 2: CD4 | `reproduce_table2_cd4.R` | `data/CD4_data.txt` | `table2_cd4_results.csv` |
| Table S.4: gallstones | `reproduce_tableS4_gallstones.R` | `data/gallstones.txt` | `tableS4_gallstones_results.csv` |

Run `Rscript real_data/reproduce_table2_cd4.R --run` or
`Rscript real_data/reproduce_tableS4_gallstones.R --run`. Sourcing either file
only defines functions. The `subject_indices` argument selects a few original
leave-one-subject-out folds for a small check; it does not subset the training
cohort. For example:

```r
source("real_data/reproduce_table2_cd4.R")
run_cd4_table(methods="LMEM", missing_rate="20", prediction_type="pointwise",
             num_cores=2, subject_indices=1:2)
```

## 📈 Figures

| Result | Data script | Plot script | Output in `results/` |
|--------|-------------|-------------|----------------------|
| Figure 2: CD4 bands | `reproduce_figure2_cd4_data.R` | `plot_figure2_cd4.R` | `figure2_cd4_band_data.csv`, `figure2_cd4_band.pdf` |
| Figure 3: CD4 density | `reproduce_figure3_cd4_density_data.R` | `plot_figure3_cd4_density.R` | `figure3_cd4_density_first100_data.csv`, `figure3_cd4_density.pdf` |
| Figure S.11: gallstones bands | `reproduce_figureS11_gallstones_data.R` | `plot_figureS11_gallstones.R` | `figureS11_gallstones_band_data.csv`, `figureS11_gallstones_band.pdf` |

Run each data script followed by its plot script. Band figures display HCP,
DWR, and LC; their supplied data and PDFs are separate from the LMEM tables.

## Method and evaluation details

The single official LMEM implementation is `R/lmem_region.R`:

- CD4: `Y ~ time + age + smoke + drug + partners + cesd + (1 + time || id)`.
- Gallstones: `t = (time - 6)/18`, then
  `Y ~ T1 + T2 + T3 + Treat + (1 + t || id)`.

Both use independent random intercept and time slope, REML on observed training
responses, and 10,000 predictive draws for each new subject. Each draw combines
fixed-effect estimation uncertainty, one random-effect vector shared by all
rows of that subject, and independent residual errors. Variance components are
plug-in estimates. Pointwise intervals use 5%/95% empirical quantiles;
simultaneous intervals use the subject-level Bonferroni allocation.

Data preparation uses Mersenne-Twister: CD4 missingness seed 123 and gallstones
imputation seed 12345. LMEM prediction uses a separate L'Ecuyer stream with seed
123 plus the held-out subject's sorted index. Gallstones evaluation uses the
completed-data target, including the saved procedure's imputed missing outcomes.
Coverage and length are averaged within each subject and then across subjects;
simultaneous coverage requires all outcomes of a held-out subject to be covered.

HCP uses `R/hcp_region.R` and its real-data grid evaluator. Real-data DWR and LC
use `R/dwr_region_realdata.R` and `R/lc_region_realdata.R`, respectively. No Oracle
benchmark is used for real data.
