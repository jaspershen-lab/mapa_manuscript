library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/02_control_data/05_benchmarking/comparison_result/all_clustering_res.rda")

jaccard_index <- function(ids_a, ids_b) {
  n_inter <- length(intersect(ids_a, ids_b))
  n_union <- length(union(ids_a, ids_b))
  if (n_union == 0) return(0)
  n_inter / n_union
}

# For one tool: compute Jaccard matrix, label each tool module by the GT module
# data ready for ggplot. Each GT module gets one column (its best-matched tool module).
# Rows = all GT modules; column label = GT module number the tool module is assigned to.
build_heatmap_data <- function(data, gt_col, tool_col, tool_name) {
  gt_mods   <- sort(unique(data[[gt_col]]))
  tool_mods <- sort(unique(data[[tool_col]]))

  jac_mat <- matrix(
    0,
    nrow = length(gt_mods),
    ncol = length(tool_mods),
    dimnames = list(as.character(gt_mods), as.character(tool_mods))
  )

  for (g in gt_mods) {
    ids_g <- data$ID[data[[gt_col]] == g]
    for (t in tool_mods) {
      ids_t <- data$ID[data[[tool_col]] == t]
      jac_mat[as.character(g), as.character(t)] <- jaccard_index(ids_g, ids_t)
    }
  }

  # All tool modules are kept; label each by the GT module it best matches.
  best_gt <- apply(jac_mat, 2, function(col) {
    as.integer(rownames(jac_mat)[which.max(col)])
  })

  result <- lapply(as.character(tool_mods), function(t_chr) {
    data.frame(
      gt_row     = as.integer(rownames(jac_mat)),
      tool_mod   = as.integer(t_chr),
      col_gt_mod = best_gt[t_chr],
      tool       = tool_name,
      jaccard    = as.numeric(jac_mat[, t_chr]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, result)
}

# Build data for all tools
hm_data <- bind_rows(
  build_heatmap_data(all_clustering_res, "ground_truth_label", "mapa_label",         "MAPA"),
  build_heatmap_data(all_clustering_res, "ground_truth_label", "paver_label",         "PAVER"),
  build_heatmap_data(all_clustering_res, "ground_truth_label", "apear_cluster",       "aPEAR"),
  build_heatmap_data(all_clustering_res, "ground_truth_label", "enrichplot_cluster",  "enrichplot")
)

# Factor levels for ordered display.
# col_id is a unique key per column (col_gt_mod + tool_mod); columns are ordered
# by GT-aligned label first, then by original tool module number as tiebreaker.
# The displayed tick label shows the GT-aligned module number.
gt_levels <- sort(unique(hm_data$gt_row))

col_order <- hm_data %>%
  distinct(tool, col_gt_mod, tool_mod) %>%
  arrange(tool, col_gt_mod, tool_mod) %>%
  mutate(col_id = paste0(tool, "_", col_gt_mod, "_", tool_mod))

col_id_levels <- col_order$col_id

col_id_labels <- setNames(
  paste0("M", col_order$col_gt_mod),
  col_order$col_id
)

hm_data <- hm_data %>%
  mutate(
    gt_row_label = factor(
      paste0("Module ", gt_row),
      levels = paste0("Module ", gt_levels)   # left-to-right on x-axis
    ),
    col_id = factor(
      paste0(tool, "_", col_gt_mod, "_", tool_mod),
      levels = col_id_levels
    ),
    tool = factor(tool, levels = c("MAPA", "PAVER", "aPEAR", "enrichplot"))
  )

# Heatmap — separate panel per tool so each carries its own fill color.
# A single shared gray-gradient legend represents the light-to-dark scale.
library(patchwork)

tool_colors <- c(
  MAPA        = "#DD5129FF",
  PAVER       = "#43B284FF",
  aPEAR       = "#0F7BA2FF",
  enrichplot  = "#FAB255FF"
)

make_tool_plot <- function(tool_name, show_x_axis) {
  df     <- dplyr::filter(hm_data, tool == tool_name)
  hi_col <- tool_colors[tool_name]

  ggplot(df, aes(x = gt_row_label, y = col_id, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.4) +
    # geom_text(
    #   data = dplyr::filter(df, jaccard > 0),
    #   aes(
    #     label = round(jaccard, 2),
    #     color = ifelse(jaccard > 0.5, "white", "black")
    #   ),
    #   size = 2.5, show.legend = FALSE
    # ) +
    scale_color_identity() +
    scale_fill_gradient(
      low = "#FFFFFF", high = hi_col,
      limits = c(0, 1), guide = "none"
    ) +
    scale_y_discrete(
      limits = rev,
      labels = col_id_labels
    ) +
    labs(
      title = tool_name,
      y     = NULL,
      x     = if (show_x_axis) "Ground truth module" else NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x      = if (show_x_axis) {
        element_text(angle = 45, hjust = 1, size = 8)
      } else {
        element_blank()
      },
      axis.ticks.x     = if (show_x_axis) element_line() else element_blank(),
      axis.text.y      = element_text(size = 8),
      strip.text       = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "#F5F5F5", color = "grey70"),
      legend.position  = "none",
      plot.title       = element_text(
        size = 11, face = "bold", hjust = 0.5, color = hi_col
      )
    )
}

# Extract a shared gray-gradient legend (white = 0, dark gray = 1)
p_leg_src <- ggplot(
  data.frame(x = 1, y = 1, z = 0.5),
  aes(x, y, fill = z)
) +
  geom_point(shape = 22, size = 5) +
  scale_fill_gradient(
    low    = "#FFFFFF",
    high   = "gray20",
    limits = c(0, 1),
    name   = "Jaccard\nIndex",
    breaks = seq(0, 1, 0.25),
    guide  = guide_colorbar(
      barwidth  = unit(0.4, "cm"),
      barheight = unit(3.5, "cm"),
      ticks.colour = "grey40"
    )
  ) +
  theme_void() +
  theme(legend.position = "right")

leg_gtable <- ggplotGrob(p_leg_src)
leg_idx    <- which(grepl("guide-box", leg_gtable$layout$name))[1]
legend_grob <- leg_gtable$grobs[[leg_idx]]

# Row counts per tool determine relative panel heights (tools are stacked)
tool_nrows <- hm_data %>%
  dplyr::distinct(tool, col_id) %>%
  dplyr::count(tool) %>%
  dplyr::arrange(match(tool, names(tool_colors)))

# Only the bottom panel shows the x-axis (GT module) labels
p_mapa       <- make_tool_plot("MAPA",       show_x_axis = FALSE)
p_paver      <- make_tool_plot("PAVER",      show_x_axis = FALSE)
p_apear      <- make_tool_plot("aPEAR",      show_x_axis = FALSE)
p_enrichplot <- make_tool_plot("enrichplot", show_x_axis = TRUE)

# Use area() so the legend grob spans all tool-panel rows on the right
n_gt <- length(unique(hm_data$gt_row))
nr   <- tool_nrows$n
cr   <- cumsum(c(0, nr))
tot  <- sum(nr)

design <- c(
  area(cr[1] + 1, 1, cr[2], n_gt),
  area(cr[2] + 1, 1, cr[3], n_gt),
  area(cr[3] + 1, 1, cr[4], n_gt),
  area(cr[4] + 1, 1, cr[5], n_gt),
  area(1, n_gt + 1, tot,    n_gt + 2)
)

p <- p_mapa + p_paver + p_apear + p_enrichplot +
  wrap_elements(legend_grob) +
  plot_layout(design = design)
  # plot_annotation(
  #   title = paste0(
  #     "Clustering consistency: Jaccard index",
  #     " between ground truth and tool modules"
  #   ),
  #   theme = theme(
  #     plot.title = element_text(size = 12, face = "bold", hjust = 0)
  #   )
  # )

p

ggsave(
  plot     = p,
  filename = "3_data_analysis/02_control_data/05_benchmarking/comparison_result/jaccard_heatmap.pdf",
  width    = 10,
  height   = 18
)
