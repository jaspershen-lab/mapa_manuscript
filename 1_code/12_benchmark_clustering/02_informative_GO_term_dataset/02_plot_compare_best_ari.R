## Revised benchmarking: compare default and optimized ARI values for the
## informative GO-term dataset.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")
script_dir <- dirname(script_file)
root <- normalizePath(file.path(script_dir, "../../.."), winslash = "/")
output_dir <- file.path(
  root,
  "3_data_analysis/12_benchmark_clustering/02_informative_GO_term_dataset"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

summary_candidates <- file.path(
  output_dir,
  c(
    "go_term_set_default_optimized_ari_settings.csv",
    "default_optimized_ari_settings.csv"
  )
)
summary_path <- summary_candidates[file.exists(summary_candidates)][[1]]
benchmark_summary <- read.csv(summary_path, check.names = FALSE)

tool_order <- c("MAPA", "PAVER", "aPEAR", "enrichplot")
tool_colors <- c(
  MAPA = "#DD5129",
  PAVER = "#43B284",
  aPEAR = "#0F7BA2",
  enrichplot = "#FAB255"
)

plot_data <- benchmark_summary |>
  filter(method %in% tool_order) |>
  transmute(
    method = factor(method, levels = rev(tool_order)),
    default_ari,
    optimized_ari
  )

stopifnot(
  nrow(plot_data) == length(tool_order),
  !anyNA(plot_data$default_ari),
  !anyNA(plot_data$optimized_ari)
)

point_data <- plot_data |>
  pivot_longer(
    cols = c(default_ari, optimized_ari),
    names_to = "type",
    values_to = "ARI"
  ) |>
  mutate(type = recode(
    type,
    default_ari = "Default",
    optimized_ari = "Optimized"
  ))

direction <- sign(plot_data$optimized_ari - plot_data$default_ari)
plot_data <- plot_data |>
  mutate(
    arrow_start = default_ari + 0.02 * direction,
    arrow_end = optimized_ari - 0.018 * direction
  )

p <- ggplot(plot_data, aes(y = method)) +
  geom_segment(
    aes(
      yend = method,
      x = arrow_start,
      xend = arrow_end,
      color = method
    ),
    arrow = arrow(length = grid::unit(0.08, "inches"), type = "closed"),
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_point(
    data = point_data,
    aes(x = ARI, color = method, fill = method, shape = type),
    size = 3,
    stroke = 0.5
  ) +
  geom_text(
    aes(x = default_ari, label = sprintf("%.3f", default_ari)),
    hjust = 1.35,
    size = 2.5
  ) +
  geom_text(
    aes(x = optimized_ari, label = sprintf("%.3f", optimized_ari)),
    hjust = -0.35,
    size = 2.5
  ) +
  scale_shape_manual(values = c(Default = 21, Optimized = 16)) +
  scale_color_manual(values = tool_colors, guide = "none") +
  scale_fill_manual(values = scales::alpha(tool_colors, 0), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.13, 0.13))) +
  guides(shape = guide_legend(
    title = NULL,
    override.aes = list(color = "black", fill = NA)
  )) +
  labs(x = "Adjusted Rand Index (ARI)", y = "Methods") +
  theme_bw() +
  theme(legend.position = "top")

ggsave(
  plot = p,
  filename = file.path(output_dir, "compare_best_ari.pdf"),
  width = 5,
  height = 6
)

