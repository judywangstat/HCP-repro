# Simulation scripts

This folder contains scripts used to reproduce the simulation results in the paper:

**Hierarchical Conformal Prediction for Clustered Data with Missing Responses**

Each script defines a local repository path near the top of the file. Please update `repo_dir` to the local path where this repository is stored before running the script.

---

## 📊 Table scripts

| Script | Output |
|--------|--------|
| `reproduce_table1.R` | Table 1 |
| `reproduce_tableS1.R` | Table S1 |
| `reproduce_tableS2.R` | Table S2 |
| `reproduce_tableS4S5_score_compare.R` | Table S4 and Table S5 |

---

## 📦 Figure data-generation scripts

| Script | Output |
|--------|--------|
| `reproduce_figureS1S2S3_data.R` | Data for Figures S1--S3 |
| `reproduce_figureS4S5S6_data.R` | Data for Figures S4--S6 |
| `reproduce_figureS7_data.R` | Data for Figure S7 |
| `reproduce_figureS8S9_data.R` | Data for Figures S8--S9 |
| `reproduce_figureS10_data.R` | Data for Figure S10 |

---

## 📈 Figure plotting scripts

| Script | Output |
|--------|--------|
| `plot_figureS1S2S3.R` | Figures S1--S3 |
| `plot_figureS4S5S6.R` | Figures S4--S6 |
| `plot_figureS7.R` | Figure S7 |
| `plot_figureS8S9.R` | Figures S8--S9 |
| `plot_figureS10.R` | Figure S10 |

---

## 📝 Notes

- Scripts beginning with `reproduce_` generate numerical results or intermediate data files.
- Scripts beginning with `plot_` generate PDF figures from saved results.
- All generated outputs are saved in the `results/` folder.
- Tables and figures with an `S` prefix refer to supplementary materials.
