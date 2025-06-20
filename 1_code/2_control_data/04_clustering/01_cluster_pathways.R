library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(tidygraph)
library(igraph)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_matrix.rda")

# Girvan-Newman clustering ====
sim.cutoff <- 0.5

edge_data <-
  embedding_sim_df %>%
  dplyr::filter(sim > sim.cutoff)

node_data <-
  control_dt %>%
  # dplyr::filter(id %in% edge_data$from | id %in% edge_data$to) %>%
  dplyr::rename(node = id) %>%
  dplyr::select(node, expected_module, expected_count, database, name)

graph_data <-
  tbl_graph(nodes = node_data,
            edges = edge_data,
            directed = FALSE,
            node_key = "node") %>%
  mutate(degree = centrality_degree())

subnetwork <-
  suppressWarnings(cluster_edge_betweenness(graph = graph_data,
                                            weights = abs(edge_attr(graph_data, "sim"))))

# cluster <-
#   paste("Functional_module", as.character(membership(subnetwork)), sep = "_")
edge_btw_clusters <- as.numeric(membership(subnetwork))

gn_graph_data <-
  graph_data %>%
  igraph::upgrade_graph() %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::mutate(gn_result = edge_btw_clusters)

save(gn_graph_data, file = "3_data_analysis/02_control_data/04_clustering/gn_graph_data.rda")

# Hierarchical clustering ====

## 1. Clustering ====
cosine_dist <- 1 - embedding_sim_matrix
cosine_dist_obj <- as.dist(cosine_dist)

hc <- hclust(cosine_dist_obj, method = "complete")
hc_cluster <- cutree(hc, k = 9)
hc_cluster_res <- data.frame(node = names(hc_cluster), hc_result = unname(hc_cluster))

# hc_clusters <-
#   paste("Functional_module", as.character(hc_cluster[node_data$node]), sep = "_")

# clustering <-
# dynamicTreeCut::cutreeDynamic(hc, distM = as.matrix(cosine_dist_obj), verbose = 0, minClusterSize = 1)
# #### Plot the hierarchical tree
# plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "pathways")
# rect.hclust(hc, k = 12, border = "red")  # Highlight 3 clusters
#### Cut the tree to assign clusters (number of clusters(k) or height(h))

## 2. Create graph object ====
edge_data <- embedding_sim_df
node_data <-
  control_dt %>%
  # dplyr::filter(id %in% edge_data$from | id %in% edge_data$to) %>%
  dplyr::rename(node = id) %>%
  dplyr::select(node, expected_module, expected_count, database, name)

hc_graph_data <-
  tbl_graph(nodes = node_data,
            edges = edge_data,
            directed = FALSE,
            node_key = "node") |>
  mutate(degree = centrality_degree())

hc_graph_data <-
  hc_graph_data %>%
  igraph::upgrade_graph() %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::left_join(hc_cluster_res, by = "node")

hc_graph_data <-
  hc_graph_data %>%
  activate(what = "edges") %>%
  filter(
    # Get the module attribute for the from node
    .N()$hc_result[from] ==
      .N()$hc_result[to]
  )

save(hc_graph_data, file = "3_data_analysis/02_control_data/04_clustering/hc_graph_data.rda")

# Markov chain clustering ====
# library(MCL)
# mcl_res <- MCL::mcl(embedding_sim_matrix,
#                     addLoops = FALSE,
#                     max.iter = 500,
#                     expansion = 2,
#                     inflation = 2.5,
#                     allow1 = TRUE)
# mcl_res$Cluster

# Spectral clustering ====
# library(Spectrum)
# spectrum_clusters <- Spectrum::Spectrum(embedding_sim_matrix,
#                                         maxk = 100,
#                                         showres = FALSE,
#                                         silent = TRUE,
#                                         clusteralg = 'km')
# spectrum_clusters$assignments

# Binary cut clustering ====
## 1. Clustering ====
library(simplifyEnrichment)

binary_cut_clusters <- simplifyEnrichment::binary_cut(mat = embedding_sim_matrix,
                                                      cutoff = 0.7)
levels(as.factor(binary_cut_clusters))
simplifyEnrichment::plot_binary_cut(mat = embedding_sim_matrix, cutoff = 0.7)

# binary_cut_clusters <- paste("Functional_module", as.character(binary_cut_clusters), sep = "_")
binary_cut_res <- data.frame(node = rownames(embedding_sim_matrix), binary_cut_result = binary_cut_clusters)

## 2. Create graph object
edge_data <- embedding_sim_df
node_data <-
  control_dt %>%
  # dplyr::filter(id %in% edge_data$from | id %in% edge_data$to) %>%
  dplyr::rename(node = id) %>%
  dplyr::select(node, expected_module, expected_count, database, name)

bc_graph_data <-
  tbl_graph(nodes = node_data,
            edges = edge_data,
            directed = FALSE,
            node_key = "node") |>
  mutate(degree = centrality_degree())

bc_graph_data <-
  bc_graph_data %>%
  igraph::upgrade_graph() %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::left_join(binary_cut_res, by = "node")

bc_graph_data <-
  bc_graph_data %>%
  activate(what = "edges") %>%
  filter(
    # Get the module attribute for the from node
    .N()$binary_cut_result[from] ==
      .N()$binary_cut_result[to]
  )

save(bc_graph_data, file = "3_data_analysis/02_control_data/04_clustering/bc_graph_data.rda")
