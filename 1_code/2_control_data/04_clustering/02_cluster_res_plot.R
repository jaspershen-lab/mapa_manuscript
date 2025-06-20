library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)
module_info <- readxl::read_excel("2_data/control_data_edge_data.xlsx", sheet = 1)
load("3_data_analysis/02_control_data/04_clustering/gn_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/bc_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/hc_graph_data.rda")

setwd("3_data_analysis/02_control_data/04_clustering/")

module_info <-
  module_info %>%
  dplyr::left_join(control_dt[, c("id", "database")], by = c("ID1" = "id")) %>%
  dplyr::left_join(control_dt[, c("id", "database")], by = c("ID2" = "id"))

###network for the raw modules
node_data <-
  control_dt %>%
  dplyr::select(id, everything()) %>%
  dplyr::mutate(expected_module =
                  stringr::str_replace(expected_module, "Functional_module_", "Module ")) %>%
  dplyr::mutate(expected_module = factor(expected_module, levels = stringr::str_sort(unique(expected_module), numeric = TRUE)))

edge_data <-
  module_info %>%
  dplyr::mutate(
    relationship = case_when(
      database.x == database.y ~ "within_database",
      database.x != database.y ~ "cross_database"
    )
  )

library(tidygraph)
library(ggraph)
library(igraph)

raw_graph_data <-
  tidygraph::tbl_graph(
    nodes = node_data,
    edges = edge_data,
    directed = FALSE,
    node_key = "id"
  ) %>%
  tidygraph::mutate(degree = centrality_degree())

# layout_df <- create_layout(raw_graph_data, layout = "fr")
# save(layout_df, file = "raw_graph_layout.rda")
load("raw_graph_layout.rda")

raw_graph_data <- raw_graph_data %>%
  activate(nodes) %>%
  mutate(x = layout_df$x, y = layout_df$y)

colors <- colorRampPalette(ggsci::pal_bmj()(n = 9))(12)

plot1 <-
  ggraph(raw_graph_data,
         layout = 'manual',
         x = x,
         y = y) +
  ggraph::geom_edge_link(aes(color = relationship),
                         alpha = 0.5,
                         show.legend = TRUE) +
  ggraph::geom_node_point(
    aes(
      size = degree,
      color = expected_module,
      shape = database
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  ggforce::geom_mark_ellipse(
    aes(
      x = x,
      y = y,
      group = expected_module,
      color = expected_module
    ),
    alpha = 1,
    expand = unit(5, "mm"),
    linewidth = 1,
    label.fontsize = 9,
    con.cap = 0,
    fill = NA,
    con.type = "straight",
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(
    aes(x = x, y = y, label = name),
    check_overlap = TRUE,
    size = 3,
    repel = TRUE
  ) +
  scale_color_manual(values = colors) +
  scale_size_continuous(range = c(3, 7)) +
  ggraph::scale_edge_color_manual(values = c(
    "within_database" = "black",
    "cross_database" = "red"
  )) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "left",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
library(extrafont)
loadfonts()
plot1
# ggsave(
#   plot = plot1,
#   filename = "raw_modules_network.pdf",
#   width = 10,
#   height = 8
# )


# Create network for Girvan-Newman clustering result
gn_graph_data <-
  gn_graph_data %>%
  activate(nodes) %>%
  dplyr::left_join(layout_df[, c("id", "x", "y")], by = c("node" = "id")) %>%
  dplyr::mutate(gn_result = paste("Module", gn_result, sep = " ")) %>%
  dplyr::mutate(gn_result = factor(gn_result, levels = stringr::str_sort(unique(gn_result), numeric = TRUE)))

plot2 <-
  ggraph(gn_graph_data,
         layout = 'manual',
         x = x,
         y = y) +
  ggraph::geom_edge_link(aes(size = sim),
                         show.legend = TRUE,
                         color = "black") +
  ggraph::geom_node_point(
    aes(
      size = degree,
      color = gn_result,
      shape = database
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  ggforce::geom_mark_ellipse(
    aes(
      x = x,
      y = y,
      group = gn_result,
      color = gn_result
    ),
    alpha = 1,
    expand = unit(5, "mm"),
    linewidth = 1,
    label.fontsize = 9,
    con.cap = 0,
    fill = NA,
    con.type = "straight",
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(
    aes(x = x, y = y, label = name),
    check_overlap = TRUE,
    size = 3,
    repel = TRUE
  ) +
  # scale_color_manual(values = colors) +
  scale_size_continuous(range = c(3, 7)) +
  ggraph::scale_edge_size_continuous(range = c(1, 3)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "left",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot2

library(extrafont)
loadfonts()

ggsave(
  plot = plot2,
  filename = "gn_clustering_network.pdf",
  width = 10,
  height = 8
)




# Create network for binary cut clustering result
bc_graph_data <-
  bc_graph_data %>%
  activate(nodes) %>%
  dplyr::left_join(layout_df[, c("id", "x", "y")], by = c("node" = "id")) %>%
  dplyr::mutate(binary_cut_result = paste("Module", binary_cut_result, sep = " ")) %>%
  dplyr::mutate(binary_cut_result = factor(binary_cut_result, levels = stringr::str_sort(unique(binary_cut_result), numeric = TRUE)))

plot2 <-
  ggraph(bc_graph_data,
         layout = 'manual',
         x = x,
         y = y) +
  ggraph::geom_edge_link(aes(size = sim),
                         show.legend = TRUE,
                         color = "black") +
  ggraph::geom_node_point(
    aes(
      size = degree,
      color = binary_cut_result,
      shape = database
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  ggforce::geom_mark_ellipse(
    aes(
      x = x,
      y = y,
      group = binary_cut_result,
      color = binary_cut_result
    ),
    alpha = 1,
    expand = unit(5, "mm"),
    linewidth = 1,
    label.fontsize = 9,
    con.cap = 0,
    fill = NA,
    con.type = "straight",
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(
    aes(x = x, y = y, label = name),
    check_overlap = TRUE,
    size = 3,
    repel = TRUE
  ) +
  # scale_color_manual(values = colors) +
  scale_size_continuous(range = c(3, 7)) +
  ggraph::scale_edge_size_continuous(range = c(1, 3)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "left",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot2

library(extrafont)
loadfonts()

ggsave(
  plot = plot2,
  filename = "bc_clustering_network.pdf",
  width = 10,
  height = 8
)






# Create network for binary cut clustering result
hc_graph_data <-
  hc_graph_data %>%
  activate(nodes) %>%
  dplyr::left_join(layout_df[, c("id", "x", "y")], by = c("node" = "id")) %>%
  dplyr::mutate(hc_result = paste("Module", hc_result, sep = " ")) %>%
  dplyr::mutate(hc_result = factor(hc_result, levels = stringr::str_sort(unique(hc_result), numeric = TRUE)))

plot3 <-
  ggraph(hc_graph_data,
         layout = 'manual',
         x = x,
         y = y) +
  ggraph::geom_edge_link(aes(size = sim),
                         show.legend = TRUE,
                         color = "black") +
  ggraph::geom_node_point(
    aes(
      size = degree,
      color = hc_result,
      shape = database
    ),
    alpha = 1,
    show.legend = TRUE
  ) +
  ggforce::geom_mark_ellipse(
    aes(
      x = x,
      y = y,
      group = hc_result,
      color = hc_result
    ),
    alpha = 1,
    expand = unit(5, "mm"),
    linewidth = 1,
    label.fontsize = 9,
    con.cap = 0,
    fill = NA,
    con.type = "straight",
    show.legend = TRUE
  ) +
  ggraph::geom_node_text(
    aes(x = x, y = y, label = name),
    check_overlap = TRUE,
    size = 3,
    repel = TRUE
  ) +
  # scale_color_manual(values = colors) +
  scale_size_continuous(range = c(3, 7)) +
  ggraph::scale_edge_size_continuous(range = c(1, 3)) +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "left",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot3

library(extrafont)
loadfonts()

ggsave(
  plot = plot3,
  filename = "hc_clustering_network.pdf",
  width = 10,
  height = 8
)
