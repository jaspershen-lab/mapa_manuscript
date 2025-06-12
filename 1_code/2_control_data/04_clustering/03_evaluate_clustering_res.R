library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(mclust)
library(tidygraph)
library(igraph)

load("3_data_analysis/02_control_data/04_clustering/bc_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/hc_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/gn_graph_data.rda")

hc_cluster_re <- data.frame(node = names(hc_cluster), hc_result = unname(hc_cluster))
spectrum_cluster_re <- data.frame(node = rownames(sim_matrix), spectrum_result = spectrum_clusters$assignments)
binary_cut_cluster_re <- data.frame(node = rownames(sim_matrix), binary_cut_result = binary_cut_clusters)
hc_cluster_dynamic <- data.frame(node = hc$labels, hc_dynamic_result = clustering)

gn_node_data <- gn_graph_data |> activate(what = "nodes") |> as_tibble()

hc_node_date <- hc_graph_data |> activate(what = "nodes") |> as_tibble()
hc_cluster_re <- hc_node_date |> select(node, hc_result)

bc_node_date <- bc_graph_data |> activate(what = "nodes") |> as_tibble()
bc_cluster_re <- bc_node_date |> select(node, binary_cut_result)

node_data <-
  gn_node_data |>
  dplyr::mutate(true_label = as.numeric(sub(pattern = "Functional_module_", replacement = "", x = expected_module))) |>
  dplyr::left_join(hc_cluster_re, by = "node") |>
  dplyr::left_join(bc_cluster_re, by = "node")

# Calculate Adjusted Rand Index ====
ari_gn <- mclust::adjustedRandIndex(node_data$true_label, node_data$gn_result)
ari_hc <- mclust::adjustedRandIndex(node_data$true_label, node_data$hc_result)
ari_bc <- mclust::adjustedRandIndex(node_data$true_label, node_data$binary_cut_result)

# aris <- c(ari_gn, ari_hc, ari_markov, ari_spectrum, ari_bc)
aris <- c("Girvan Newman" = ari_gn, "Hierarchical" = ari_hc, "Binary cut" = ari_bc)
aris_sorted <- sort(aris, decreasing = TRUE)

ari_result <- data.frame("Clustering methods" = names(aris_sorted),
                         "Adjusted Rand Index" = unname(aris_sorted))
p <-
  ari_result |>
  ggplot(aes(x = Clustering.methods, y = Adjusted.Rand.Index)) +
  geom_col(width = 0.6) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12)   # ← set x‐axis tick label size
  )  +
  ggplot2::annotate(
    "text",
    x = ari_result$Clustering.methods,
    y = ari_result$Adjusted.Rand.Index + 0.02,
    label = sprintf(
      "%.2f",
      ari_result$Adjusted.Rand.Index
    ),
    hjust = 0.5,
    vjust = 0,
    fontface = "bold",
    size = 4
  ) +
  labs(
    x = "Clustering Methods",
    y = "Adjusted Rand Index"
  )

ggsave(plot = p, filename = "3_data_analysis/02_control_data/04_clustering/clustering_evaluation_plot_ARI.pdf",
       width = 8,
       height = 6)
