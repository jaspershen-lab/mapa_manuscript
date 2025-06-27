library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(mclust)
library(tidygraph)
library(igraph)

load("3_data_analysis/02_control_data/04_clustering/bc_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/hc_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/louvain_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/louvain_graph_data.rda")

setwd("3_data_analysis/02_control_data/04_clustering/")

louvain_node_data <- louvain_graph_data |> activate(what = "nodes") |> as_tibble()

hc_node_data <- hc_graph_data |> activate(what = "nodes") |> as_tibble()
hc_cluster_re <- hc_node_data |> select(node, hc_result)

bc_node_data <- bc_graph_data |> activate(what = "nodes") |> as_tibble()
bc_cluster_re <- bc_node_data |> select(node, binary_cut_result)

node_data <-
  louvain_node_data |>
  dplyr::mutate(true_label = as.numeric(sub(pattern = "Functional_module_", replacement = "", x = expected_module))) |>
  dplyr::left_join(hc_cluster_re, by = "node") |>
  dplyr::left_join(bc_cluster_re, by = "node")

# Calculate Adjusted Rand Index ====
library(mclust)
ari_gn <- mclust::adjustedRandIndex(node_data$true_label, node_data$louvain_result)
ari_hc <- mclust::adjustedRandIndex(node_data$true_label, node_data$hc_result)
ari_bc <- mclust::adjustedRandIndex(node_data$true_label, node_data$binary_cut_result)

node_data %>%
  dplyr::select(node, true_label, louvain_result, hc_result, binary_cut_result)

x <- c(1,1,1,2,2,2,3,3,3)
y <- c(3,3,2,2,2,4,2,1,1)
mclust::adjustedRandIndex(x,y)

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
p
ggsave(plot = p, filename = "clustering_evaluation_plot_ARI.pdf",
       width = 8,
       height = 6)
