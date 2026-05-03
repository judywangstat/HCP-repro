# Real-data analysis scripts

This folder contains scripts used to reproduce the real-data analyses in the paper:

**Hierarchical Conformal Prediction for Clustered Data with Missing Responses**

Each script defines a local repository path near the top of the file. Please update `repo_dir` to the local path where this repository is stored before running the script.

---

## 📊 Table scripts

| Script | Output |
|--------|--------|
| `reproduce_table2_cd4.R` | Table 2 (CD4 dataset) |
| `reproduce_tableS3_gallstones.R` | Table S3 (gallstones dataset) |

---

## 📦 Figure data-generation scripts

| Script | Output |
|--------|--------|
| `reproduce_figure2_cd4_data.R` | Data for Figure 2 (CD4 prediction bands) |
| `reproduce_figure3_cd4_density_data.R` | Data for Figure 3 (CD4 conditional density) |
| `reproduce_figureS11_gallstones_data.R` | Data for Figure S11 (gallstones prediction bands) |

---

## 📈 Figure plotting scripts

| Script | Output |
|--------|--------|
| `plot_figure2_cd4.R` | Figure 2 (CD4 prediction bands) |
| `plot_figure3_cd4_density.R` | Figure 3 (CD4 conditional density) |
| `plot_figureS11_gallstones.R` | Figure S11 (gallstones prediction bands) |

---

## 📝 Notes

- Scripts beginning with `reproduce_` generate intermediate data files.
- Scripts beginning with `plot_` generate final figures (PDF format).
- All outputs are saved in the `results/` folder.
- Tables and figures with an `S` prefix refer to supplementary materials.
