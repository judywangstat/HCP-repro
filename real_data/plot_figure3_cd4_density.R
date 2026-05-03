# ============================================================
# Plot Figure 3
# CD4 conditional density example
# HCPclust-repro
#
# This script plots the selected CD4 conditional density example:
#   Subject index = 17
#   id = 10191
#   time index = 7
#
# Input:
#   results/figure3_cd4_density_first100_data.csv
#
# Output:
#   results/figure3_cd4_density.pdf
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

input_file <- file.path(results_dir, "figure3_cd4_density_first100_data.csv")
output_file <- file.path(results_dir, "figure3_cd4_density.pdf")

# ------------------------------------------------------------
# Selected example
# ------------------------------------------------------------
selected_subject_index <- 17
selected_id <- 10191
selected_time_index <- 7

# ------------------------------------------------------------
# Read density data
# ------------------------------------------------------------
dens_dat <- read.csv(input_file, stringsAsFactors = FALSE)

plot_dat <- dens_dat %>%
  filter(
    subject_index == selected_subject_index,
    id == selected_id,
    time_index == selected_time_index
  ) %>%
  arrange(y_grid)

if (nrow(plot_dat) == 0L) {
  stop(
    "No density data found for subject_index = ", selected_subject_index,
    ", id = ", selected_id,
    ", time_index = ", selected_time_index
  )
}

# ------------------------------------------------------------
# Smooth and normalize density
# ------------------------------------------------------------
smoothed <- smooth.spline(
  x = plot_dat$y_grid,
  y = plot_dat$density,
  spar = 0.6
)

smoothed_y <- pmax(smoothed$y, 1e-5)
smoothed_y <- smoothed_y / sum(smoothed_y)

plot_smooth <- data.frame(
  y_grid = smoothed$x,
  density = smoothed_y
)

# ------------------------------------------------------------
# Construct density cutoff, shaded prediction region, and endpoints
# ------------------------------------------------------------
h_pool <- smoothed_y[smoothed_y > 4.2e-03]

if (length(h_pool) == 0L) {
  h_value <- quantile(smoothed_y, 0.75, na.rm = TRUE)
} else {
  h_value <- quantile(h_pool, 0.1, na.rm = TRUE)
}

above_idx <- which(smoothed_y > h_value)

if (length(above_idx) == 0L) {
  stop("No density values exceed the horizontal cutoff.")
}

solution1 <- above_idx[1]
solution2 <- above_idx[length(above_idx)]

ribbon_dat <- data.frame(
  y_grid = smoothed$x[solution1:solution2],
  ymin = 0,
  ymax = smoothed_y[solution1:solution2]
)

boundary_dat <- data.frame(
  x = c(smoothed$x[solution1], smoothed$x[solution2]),
  y0 = c(0, 0),
  y1 = c(h_value, h_value)
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(plot_smooth, aes(x = y_grid, y = density)) +
  geom_ribbon(
    data = ribbon_dat,
    aes(x = y_grid, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey75",
    alpha = 0.65
  ) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_hline(
    yintercept = h_value,
    linetype = "dashed",
    linewidth = 0.55,
    color = "black"
  ) +
  geom_segment(
    data = boundary_dat,
    aes(x = x, xend = x, y = y0, yend = y1),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.55,
    color = "black"
  ) +
  geom_point(
    data = boundary_dat,
    aes(x = x, y = 0),
    inherit.aes = FALSE,
    size = 1.7,
    color = "black"
  ) +
  labs(
    x = "CD4+ cell counts",
    y = "Conditional density of CD4+ cell counts"
  ) +
  coord_cartesian(
    xlim = c(-50, 1900),
    ylim = c(0, 0.03)
  ) +
  scale_x_continuous(
    breaks = c(0, 500, 1000, 1500)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.01, 0.02, 0.03)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey88"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

print(p)

# ------------------------------------------------------------
# Save output
# ------------------------------------------------------------
ggsave(
  filename = output_file,
  plot = p,
  width = 6.5,
  height = 4,
  device = grDevices::pdf
)

cat("Saved Figure 3 to:\n", output_file, "\n")
cat("Subject index:", selected_subject_index, "\n")
cat("Subject id:", selected_id, "\n")
cat("Time index:", selected_time_index, "\n")
cat("h value:", h_value, "\n")