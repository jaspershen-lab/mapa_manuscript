library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Create network for Girvan-Newman clustering result
load("3_data_analysis/02_control_data/04_clustering/gn_graph_data.rda")

colors <- colorRampPalette(c('#0ca9ce', '#78cfe5', '#c6ecf1', '#ff6f81', '#ff9c8f', '#ffc2c0','#d386bf',
                             '#cdb1d2', '#fae6f0', '#eb6fa6', '#ff88b5', '#00b1a5',"#ffa68f","#ffca75","#97bc83","#acd295",
                             "#00ada1","#009f93","#ace2da","#448c99","#00b3bc","#b8d8c9","#db888e","#e397a4","#ead0c7",
                             "#8f9898","#bfcfcb"))(12)

gn_plot <-
  gn_graph_data %>%
  ggraph::ggraph(layout = 'fr',
                 circular = FALSE) +
  ggforce::geom_mark_ellipse(
    aes(x = x,
        y = y,
        group = expected_module,
        color = expected_module),
    alpha = 1,
    expand = unit(5, "mm"),
    linewidth = 1,
    label.fontsize = 9,
    con.cap = 0,
    fill = NA,
    con.type = "straight",
    show.legend = TRUE
  ) +
  ggraph::geom_edge_link(
    aes(width = sim),
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = database),
    size = 6,
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = database_color) +
  guides(
    color = guide_legend(order = 1, ncol = 1, title = "Expected_functional_module", override.aes = list(linewidth = 0.5)),
    fill = guide_legend(order = 2, ncol = 1, title = "Database", override.aes = list(linewidth = 0)),
    edge_width = guide_legend(order = 3, title = "Similarity", override.aes = list(linewidth = 0))) +
  ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  labs(fill = "Database") +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "left",
    legend.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggraph::geom_node_text(
    aes(x = x,
        y = y,
        label = stringr::str_wrap(paste0("Functional_module_", as.character(gn_result), ":", name), width = 30)),
    check_overlap = TRUE,
    size = 3,
    repel = TRUE)

library(Cairo)
CairoPDF("3_data_analysis/02_control_data/04_clustering/gn_clustering_network.pdf", width=18, height=18)
gn_plot
dev.off()
