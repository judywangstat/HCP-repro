# Real-data analysis scripts

CD4 and gallstones analysis scripts. Data are in `data/`, helpers in `R/`, and
results in `results/`. Run from the repository root or set `HCP_REPO_DIR`.
Each script has a `repo_dir` setting.

## 📊 Tables

| Result | Script | Input | Output in `results/` |
|--------|--------|-------|----------------------|
| Table 2: CD4 | `reproduce_table2_cd4.R` | `data/CD4_data.txt` | `table2_cd4_results.csv` |
| Table S.4: gallstones | `reproduce_tableS4_gallstones.R` | `data/gallstones.txt` | `tableS4_gallstones_results.csv` |

Run `Rscript real_data/reproduce_table2_cd4.R --run` or
`Rscript real_data/reproduce_tableS4_gallstones.R --run`. Sourcing either file
only defines functions. Use `subject_indices` for a small check of selected
leave-one-subject-out folds. Each fold still uses the other subjects for training:

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

Run the data script first, then the plotting script. Band figures show HCP,
DWR, and LC. Their data and PDFs are separate from the LMEM tables.

## Analysis details

LMEM is implemented in `R/lmem_region.R`:

- CD4: `Y ~ time + age + smoke + drug + partners + cesd + (1 + time || id)`.
- Gallstones: `t = (time - 6)/18`, then
  `Y ~ T1 + T2 + T3 + Treat + (1 + t || id)`.

Both models use independent random intercept and time slope and are fit by REML
using observed training responses. New-subject intervals use 10,000 draws with
fixed-effect uncertainty, subject-level random effects, and residual error.
Pointwise intervals use the 5th and 95th percentiles. Simultaneous intervals use
the subject-level Bonferroni adjustment.

CD4 missingness uses seed 123; gallstones imputation uses seed 12345. Both use
Mersenne-Twister. LMEM prediction uses a separate L'Ecuyer stream and restores
the caller's RNG state. Gallstones evaluation includes imputed outcomes.
Coverage and length are averaged within each subject, then across subjects.
Simultaneous coverage requires all outcomes of a held-out subject to be covered.

HCP uses `R/hcp_region.R` with its grid evaluator. DWR and LC use
`R/dwr_region_realdata.R` and `R/lc_region_realdata.R`. Oracle is not used for
real-data analyses.
