## Benchmarking on the curated pathway dataset.
## Step 2: Plot the optimized clustering consistency heatmap.

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(cowplot)
library(patchwork)

output_dir <-
  "3_data_analysis/12_benchmark_clustering/01_curated_pathway_dataset"

hm_data <-
  readr::read_csv(
    file.path(output_dir, "optimized_jaccard_heatmap_data.csv"),
    show_col_types = FALSE
  )
benchmark_summary <-
  readr::read_csv(
    file.path(output_dir, "default_optimized_ari_settings.csv"),
    show_col_types = FALSE
  )

expected_tools <- c("MAPA", "PAVER", "aPEAR", "enrichplot")
ari_by_tool <- stats::setNames(
  benchmark_summary$optimized_ari,
  benchmark_summary$method
)
stopifnot(
  all(expected_tools %in% hm_data$tool),
  all(expected_tools %in% names(ari_by_tool)),
  all(is.finite(ari_by_tool[expected_tools]))
)

tool_order <- names(sort(ari_by_tool[expected_tools], decreasing = TRUE))
tool_colors <-
  c(
    MAPA = methods[["mapa_cluster_label"]],
    PAVER = methods[["paver_cluster_label"]],
    aPEAR = methods[["apear_cluster_label"]],
    enrichplot = methods[["enrich_plot_cluster_label"]]
  )
gt_levels <- sort(unique(hm_data$ground_truth_label))
plot_list <- list()
row_counts <- integer()

for (current_tool in tool_order) {
  df <- hm_data |>
    filter(.data$tool == .env$current_tool) |>
    mutate(
      gt_factor = factor(ground_truth_label, levels = gt_levels),
      row_factor = factor(
        display_row,
        levels = rev(sort(unique(display_row)))
      )
    )

  stopifnot(
    nrow(df) == sum(hm_data$tool == current_tool),
    n_distinct(df$tool) == 1L
  )

  n_rows <- n_distinct(df$display_row)
  row_counts <- c(row_counts, n_rows)
  show_y <- n_rows <= 80
  y_labels <- setNames(
    paste0(
      "M",
      df$best_ground_truth[
        match(sort(unique(df$display_row)), df$display_row)
      ]
    ),
    as.character(sort(unique(df$display_row)))
  )
  panel_title <- sprintf(
    "%s (ARI = %.3f)",
    current_tool,
    ari_by_tool[[current_tool]]
  )

  plot_list[[current_tool]] <-
    ggplot(df, aes(gt_factor, row_factor, fill = jaccard)) +
    geom_tile(
      color = "white",
      linewidth = if (n_rows <= 80) 0.25 else 0
    ) +
    scale_fill_gradient(
      low = "white",
      high = tool_colors[[current_tool]],
      limits = c(0, 1),
      guide = "none"
    ) +
    scale_y_discrete(labels = if (show_y) y_labels else NULL) +
    labs(title = panel_title, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(
        color = tool_colors[[current_tool]],
        face = "bold",
        hjust = 0.5,
        size = 11
      ),
      axis.text.y = if (show_y) element_text(size = 5) else element_blank(),
      axis.ticks.y = if (show_y) {
        element_line(linewidth = 0.2)
      } else {
        element_blank()
      },
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 4, 2, 4)
    )
}

bottom_tool <- tool_order[[length(tool_order)]]
plot_list[[bottom_tool]] <- plot_list[[bottom_tool]] +
  labs(x = "Ground-truth module") +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = if (length(gt_levels) > 30) 5 else 7
    ),
    axis.ticks.x = element_line(linewidth = 0.2)
  )

stacked <- patchwork::wrap_plots(
  plot_list[tool_order],
  ncol = 1,
  heights = pmax(row_counts, 4)
) +
  patchwork::plot_annotation(
    title = "Curated pathway dataset: optimized clustering consistency"
  )

legend_plot <-
  ggplot(data.frame(x = 1, y = 1, z = 0.5), aes(x, y, fill = z)) +
  geom_point(shape = 22, size = 4) +
  scale_fill_gradient(
    low = "white",
    high = "grey20",
    limits = c(0, 1),
    name = "Jaccard\nindex",
    breaks = seq(0, 1, 0.25)
  ) +
  theme_void() +
  theme(legend.position = "right")

legend <- cowplot::get_legend(legend_plot)
combined <- cowplot::plot_grid(
  stacked,
  legend,
  nrow = 1,
  rel_widths = c(0.94, 0.06)
)

total_rows <- sum(row_counts)
plot_height <- min(30, max(9, 4 + total_rows * 0.055))

combined

ggsave(
  filename = file.path(output_dir, "optimized_heatmap_consistency.pdf"),
  plot = combined,
  width = 6,
  height = plot_height-1,
  limitsize = FALSE
)

message("Method order: ", paste(tool_order, collapse = " > "))
