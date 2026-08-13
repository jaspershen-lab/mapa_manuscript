## Benchmarking on the informative GO term dataset (k = 20)
## Step 11: Visualize the ontology-derived ground-truth modules.
##
## Each reference module is displayed as a star-shaped network. The central
## diamond is the shared direct-parent GO term (module anchor), which is used
## only to define the reference label and is not part of the benchmark input.
## The surrounding circles are the 313 informative GO terms evaluated by the
## clustering methods.

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

library(tidygraph)
library(ggraph)
library(ggforce)

input_dir <- "2_data/informative_go_set/k_20"
output_dir <- "3_data_analysis/11_informative_go_benchmarking"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ground_truth_dt <-
  readr::read_csv(file.path(input_dir, "ground_truth.csv"),
                  show_col_types = FALSE)

module_dt <-
  readr::read_csv(file.path(input_dir, "modules.csv"),
                  show_col_types = FALSE) |>
  dplyr::select(module_id, anchor_go_id, anchor_name, module_size)

## Dataset-level checks recorded in the k = 20 manifest.
stopifnot(
  nrow(ground_truth_dt) == 313,
  dplyr::n_distinct(ground_truth_dt$module_id) == 57,
  !anyDuplicated(ground_truth_dt$go_id),
  nrow(module_dt) == 57,
  sum(module_dt$module_size) == 313,
  setequal(ground_truth_dt$module_id, module_dt$module_id)
)

## Arrange the 57 disconnected modules on a fixed grid. This avoids the
## run-to-run variation and component overlap that can occur with force-directed
## layouts, while retaining the raw-module-network appearance used in Fig. 3c.
module_layout <-
  module_dt |>
  dplyr::arrange(dplyr::desc(module_size), anchor_go_id) |>
  dplyr::mutate(
    module_order = dplyr::row_number(),
    grid_column = (module_order - 1L) %% 8L,
    grid_row = (module_order - 1L) %/% 8L,
    centre_x = grid_column * 5.2,
    centre_y = -grid_row * 4.6,
    module_label = paste0(
      stringr::str_wrap(anchor_name, width = 23),
      "\n", anchor_go_id, " (n = ", module_size, ")"
    )
  )

member_layout <-
  ground_truth_dt |>
  dplyr::left_join(
    module_layout |>
      dplyr::select(module_id, module_size, centre_x, centre_y),
    by = "module_id"
  ) |>
  dplyr::group_by(module_id) |>
  dplyr::arrange(go_id, .by_group = TRUE) |>
  dplyr::mutate(
    member_order = dplyr::row_number(),
    angle = 2 * pi * (member_order - 1) / dplyr::n() + pi / 2,
    radius = 0.82 + 0.018 * dplyr::n(),
    x = centre_x + radius * cos(angle),
    y = centre_y + radius * sin(angle),
    node_id = go_id,
    node_label = term_name,
    node_type = "Informative GO term"
  ) |>
  dplyr::ungroup() |>
  dplyr::select(node_id, node_label, node_type, module_id, module_size, x, y)

anchor_layout <-
  module_layout |>
  dplyr::transmute(
    node_id = anchor_go_id,
    node_label = anchor_name,
    node_type = "Reference module anchor",
    module_id,
    module_size,
    x = centre_x,
    y = centre_y
  )

node_layout <-
  dplyr::bind_rows(anchor_layout, member_layout) |>
  dplyr::mutate(
    node_type = factor(
      node_type,
      levels = c("Reference module anchor", "Informative GO term")
    )
  )

edge_dt <-
  ground_truth_dt |>
  dplyr::transmute(
    from = anchor_go_id,
    to = go_id,
    module_id
  )

ground_truth_graph <-
  tidygraph::tbl_graph(
    nodes = node_layout,
    edges = edge_dt,
    directed = FALSE,
    node_key = "node_id"
  )

## Use a continuous qualitative palette to distinguish adjacent modules. The
## module colors are intentionally not shown as a 57-item legend; each module is
## identified directly by its anchor label.
module_ids <- module_layout$module_id
module_colors <- grDevices::hcl.colors(
  n = length(module_ids), palette = "Dynamic", alpha = 1
)
names(module_colors) <- module_ids

ground_truth_plot <-
  ggraph::ggraph(
    ground_truth_graph,
    layout = "manual",
    x = x,
    y = y
  ) +
  ggforce::geom_mark_ellipse(
    data = node_layout,
    aes(x = x, y = y, group = module_id,
        color = module_id, fill = module_id),
    expand = unit(2.0, "mm"),
    linewidth = 0.35,
    alpha = 0.055,
    show.legend = FALSE
  ) +
  ggraph::geom_edge_link(
    aes(edge_color = module_id),
    edge_width = 0.35,
    edge_alpha = 0.55,
    show.legend = FALSE
  ) +
  ggraph::geom_node_point(
    aes(shape = node_type, size = node_type, fill = module_id),
    color = "white",
    stroke = 0.35,
    show.legend = TRUE
  ) +
  geom_text(
    data = module_layout,
    aes(x = centre_x, y = centre_y - 1.48, label = module_label),
    size = 2.15,
    lineheight = 0.88,
    color = "black",
    vjust = 1
  ) +
  scale_color_manual(values = module_colors, guide = "none") +
  scale_fill_manual(values = module_colors, guide = "none") +
  ggraph::scale_edge_color_manual(values = module_colors, guide = "none") +
  scale_shape_manual(
    name = NULL,
    values = c(
      "Reference module anchor" = 23,
      "Informative GO term" = 21
    )
  ) +
  scale_size_manual(
    name = NULL,
    values = c(
      "Reference module anchor" = 4.0,
      "Informative GO term" = 2.15
    )
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        size = c(4.0, 2.5),
        fill = "grey55",
        color = "white"
      )
    ),
    size = "none"
  ) +
  coord_equal(clip = "off") +
  ggraph::theme_graph(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 8, 8, 8)
  )

ground_truth_plot

## Save coordinates so the exact layout can be reused in later panels.
readr::write_csv(
  node_layout |>
    dplyr::mutate(node_type = as.character(node_type)),
  file.path(output_dir, "informative_go_ground_truth_network_layout.csv")
)

ggsave(
  plot = ground_truth_plot,
  filename = file.path(output_dir, "informative_go_ground_truth_network.pdf"),
  width = 20,
  height = 18,
  units = "in",
  device = cairo_pdf,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  plot = ground_truth_plot,
  filename = file.path(output_dir, "informative_go_ground_truth_network.png"),
  width = 20,
  height = 18,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
