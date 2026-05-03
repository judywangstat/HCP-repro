# ============================================================
# Plot Figure 2
# Selected CD4 real-data prediction bands
# HCPclust-repro
#
# This script reads the processed CD4 prediction-band data and
# plots selected subjects for Figure 2 in the paper. The figure
# compares HCP, DWR, and LC prediction bands against the observed
# CD4 trajectories.
#
# Input:
#   results/figure2_cd4_band_data.csv
#
# Output:
#   results/figure2_cd4_band.pdf
#
# Note:
#   Update repo_dir below to the local path of this repository
#   before running the script.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# User-adjustable paths
# ------------------------------------------------------------
repo_dir <- normalizePath("~/Desktop/HCPclust-repro", mustWork = TRUE)
results_dir <- file.path(repo_dir, "results")

input_file <- file.path(results_dir, "figure2_cd4_band_data.csv")
output_file <- file.path(results_dir, "figure2_cd4_band.pdf")

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
fig_dat <- read.csv(input_file, stringsAsFactors = FALSE)

required_cols <- c(
  "id", "time", "true", "lower", "upper",
  "method", "missing_rate", "prediction_type"
)

missing_cols <- setdiff(required_cols, names(fig_dat))
if (length(missing_cols) > 0L) {
  stop("Input file is missing columns: ", paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------
# Selected IDs
# ------------------------------------------------------------
selected_ids <- c(41566, 20344)

plot_dat <- fig_dat %>%
  filter(id %in% selected_ids, method %in% c("HCP", "DWR", "LC")) %>%
  mutate(
    id = factor(id, levels = selected_ids),
    id_label = paste0("ID = ", id),
    method = factor(method, levels = c("HCP", "DWR", "LC")),
    time = as.numeric(time),
    true = as.numeric(true),
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  ) %>%
  arrange(id, method, time)

if (nrow(plot_dat) == 0L) {
  stop("No data found for selected_ids.")
}

missing_ids <- setdiff(selected_ids, unique(as.numeric(as.character(plot_dat$id))))
if (length(missing_ids) > 0L) {
  warning("No data found for selected IDs: ", paste(missing_ids, collapse = ", "))
}

true_dat <- plot_dat %>%
  distinct(id, id_label, time, true)

# ------------------------------------------------------------
# Plot all selected IDs
# ------------------------------------------------------------
p <- ggplot() +
  geom_ribbon(
    data = filter(plot_dat, method == "HCP"),
    aes(x = time, ymin = lower, ymax = upper, fill = method),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(
    data = filter(plot_dat, method == "HCP"),
    aes(x = time, y = lower),
    color = "steelblue4",
    linewidth = 0.6
  ) +
  geom_line(
    data = filter(plot_dat, method == "HCP"),
    aes(x = time, y = upper),
    color = "steelblue4",
    linewidth = 0.6
  ) +
  geom_line(
    data = filter(plot_dat, method == "DWR"),
    aes(x = time, y = lower, color = method, linetype = method),
    linewidth = 0.8
  ) +
  geom_line(
    data = filter(plot_dat, method == "DWR"),
    aes(x = time, y = upper, color = method, linetype = method),
    linewidth = 0.8
  ) +
  geom_line(
    data = filter(plot_dat, method == "LC"),
    aes(x = time, y = lower, color = method, linetype = method),
    linewidth = 0.8
  ) +
  geom_line(
    data = filter(plot_dat, method == "LC"),
    aes(x = time, y = upper, color = method, linetype = method),
    linewidth = 0.8
  ) +
  geom_point(
    data = true_dat,
    aes(x = time, y = true),
    color = "gray20",
    size = 1.5
  ) +
  facet_wrap(~ id_label, ncol = 2, scales = "free_x") +
  labs(
    x = "Years since seroconversion",
    y = "CD4+ cell counts"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(
    values = c("HCP" = "steelblue3"),
    labels = c("HCP")
  ) +
  scale_color_manual(
    values = c(
      "DWR" = "firebrick4",
      "LC"  = "darkolivegreen4"
    ),
    labels = c("DWR", "LC")
  ) +
  scale_linetype_manual(
    values = c(
      "DWR" = "dashed",
      "LC"  = "dotdash"
    )
  ) +
  coord_cartesian(ylim = c(-300, 2600)) +
  theme(
    legend.position = "right",
    legend.key.width = unit(0.9, "cm"),
    legend.title = element_blank(),
    legend.spacing.y = unit(0.1, "lines"),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(4, "lines"),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.12, "cm"),
    plot.margin = margin(6, 8, 6, 6),
    strip.text = element_text(face = "bold")
  ) +
  guides(
    fill = guide_legend(
      title = NULL,
      order = 1,
      override.aes = list(
        fill = "steelblue3",
        alpha = 0.2,
        color = "steelblue4",
        linewidth = 0.6
      )
    ),
    color = guide_legend(title = NULL, order = 2),
    linetype = guide_legend(title = NULL, order = 2)
  )

print(p)

# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------
ggsave(
  filename = output_file,
  plot = p,
  width = 9,
  height = 3.5,
  device = "pdf"
)

cat("Saved selected-ID figure to:\n", output_file, "\n")