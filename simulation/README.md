# Simulation scripts

This folder contains scripts used to reproduce the simulation results in the paper:

**Hierarchical Conformal Prediction for Clustered Data with Missing Responses**

Each script defines a local repository path near the top of the file. Please update `repo_dir` to the local path where this repository is stored before running the script.

---

## Table scripts

- `reproduce_table1.R`  
  Reproduce Table 1.

- `reproduce_tableS1.R`  
  Reproduce Table S1.

- `reproduce_tableS2.R`  
  Reproduce Table S2.

- `reproduce_tableS4S5_score_compare.R`  
  Reproduce Table S4 and Table S5.

---

## Figure data-generation scripts

- `reproduce_figureS1S2S3_data.R`  
  Generate data for Figures S1--S3.

- `reproduce_figureS4S5S6_data.R`  
  Generate data for Figures S4--S6.

- `reproduce_figureS7_data.R`  
  Generate data for Figure S7.

- `reproduce_figureS8S9_data.R`  
  Generate data for Figures S8--S9.

- `reproduce_figureS10_data.R`  
  Generate data for Figure S10.

---

## Figure plotting scripts

- `plot_figureS1S2S3.R`  
  Plot Figures S1--S3.

- `plot_figureS4S5S6.R`  
  Plot Figures S4--S6.

- `plot_figureS7.R`  
  Plot Figure S7.

- `plot_figureS8S9.R`  
  Plot Figures S8--S9.

- `plot_figureS10.R`  
  Plot Figure S10.

---

## Notes

- Scripts with names beginning with `reproduce_` generate numerical results or intermediate data files.
- Scripts with names beginning with `plot_` generate PDF figures from saved results.
- All generated outputs are saved in the `results/` folder.
- Tables and figures with an `S` prefix refer to supplementary materials.