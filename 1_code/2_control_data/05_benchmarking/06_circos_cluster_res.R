library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

setwd("3_data_analysis/02_control_data/05_benchmarking")

load("4_methods_all_best_ari.rda")

# Best enrichplot ====
load("enriched_result.rda")
library(enrichplot)

edo <- enrichplot::pairwise_termsim(enriched_result)
set.seed(123)
# Generate emapplot with current nCluster
p <- enrichplot::emapplot(
  x = edo,
  showCategory = 44,
  layout = "fr",
  min_edge = 0.2,
  color_edge = "grey",
  node_label = "group",
  group = TRUE,
  clusterFunction = stats::kmeans,
  nCluster = 17
)

# Extract clustering results
clustering_data <- p$data
cluster_results <- clustering_data[, c("name", "color2")]
colnames(cluster_results) <- c("Description", "cluster_label")

# Join with enriched result to get IDs
cluster_results <- cluster_results |>
  dplyr::left_join(edo@result[, c("ID", "Description")],
                   by = "Description")
# Create numerical cluster IDs
unique_clusters <- unique(cluster_results$cluster_label)
cluster_mapping <- data.frame(
  cluster_label = unique_clusters,
  cluster_id = 1:length(unique_clusters)
)

cluster_results <- merge(cluster_results, cluster_mapping, by = "cluster_label")

enrichplot_cluster_res <- cluster_results |>
  select(ID, cluster_id) |>
  rename(enrichplot_cluster = cluster_id)

# Best aPEAR ====
library(aPEAR)
apear_res_best <- enrichmentNetwork(
  enrichment = enriched_result@result,
  simMethod = "jaccard",
  innerCutoff = 0.1,
  outerCutoff = 0.5,
  clustMethod = "hier",
  clustNameMethod = "pagerank",
  colorBy = "pvalue",
  colorType = "pval",
  nodeSize = "Count",
  plotOnly = FALSE,
  drawEllipses = TRUE,
  minClusterSize = 2,
  fontSize = 3
)

p <- apear_res_best$plot
p

# Extract cluster assignments
apear_cluster <- apear_res_best$clusters |>
  dplyr::rename(Description = ID)

# Create pathway ID to description mapping
pathway_id_desc <- enriched_result@result |>
  dplyr::select(ID, Description) |>
  dplyr::left_join(apear_cluster, by = "Description")

# Create cluster mapper (assign numeric labels to clusters)
cluster_mapper <- data.frame(
  Cluster = unique(apear_cluster$Cluster),
  cluster_label = 1:length(unique(apear_cluster$Cluster))
)

# Join with cluster mapper
all_cluster_label <- pathway_id_desc |>
  dplyr::left_join(cluster_mapper, by = "Cluster") |>
  dplyr::arrange(cluster_label)

# Assign sequential labels to singletons (NA values)
na_indices <- which(is.na(all_cluster_label$cluster_label))
if (length(na_indices) > 0) {
  max_label <- max(all_cluster_label$cluster_label, na.rm = TRUE)
  singleton_labels <- seq(max_label + 1, max_label + length(na_indices))
  all_cluster_label$cluster_label[na_indices] <- singleton_labels
}

apear_cluster_res <- all_cluster_label |>
  select(ID, cluster_label) |>
  rename(apear_cluster = cluster_label)

# Best PAVER ====
library(PAVER)
load("PAVER_01_input.rda")
load("PAVER_02_embedding_matrix.rda")
load("PAVER_03_term2name.rda")

PAVER_result <- prepare_data(input, embedding_matrix, term2name)

minClusterSize <- 1
maxCoreScatter <- 0.22
minGap <- (1 - maxCoreScatter) * 3 / 4
PAVER_result <- generate_themes(
  PAVER_result,
  maxCoreScatter = maxCoreScatter,
  minGap = minGap,
  minClusterSize = minClusterSize
)

# Export clustered data
clustered_input <- PAVER_export(PAVER_result)

# Check if clustering was successful (at least 2 clusters)
unique_clusters <- unique(clustered_input$Cluster)
nClust <- length(unique_clusters)

# Create cluster mapping
cluster_mapping <- data.frame(
  Cluster = unique_clusters,
  paver_cluster_label = seq_len(nClust)
)

# Join with cluster mapping
clustered_input <- clustered_input |>
  dplyr::left_join(cluster_mapping, by = "Cluster")

paver_cluster_res <- clustered_input |>
  select(GOID, paver_cluster_label) |>
  rename(ID = GOID,
         paver_label = paver_cluster_label)


# Best MAPA ====
setwd(get_project_wd())
load("3_data_analysis/02_control_data/04_clustering/louvain_graph_data.rda")
mapa_cluster_res <- louvain_graph_data |>
  tidygraph::activate(what = "nodes") |>
  as_tibble()

mapa_cluster_res <-
  mapa_cluster_res |>
  select(node, louvain_result) |>
  rename(ID = node, mapa_label = louvain_result)

load("3_data_analysis/02_control_data/05_benchmarking/ground_truth_dt.rda")
ground_truth_dt <- ground_truth_dt |> rename(ID = id)

all_clustering_res <-
  ground_truth_dt |>
  left_join(mapa_cluster_res, by = "ID") |>
  left_join(paver_cluster_res, by = "ID") |>
  left_join(apear_cluster_res, by = "ID") |>
  left_join(enrichplot_cluster_res, by = "ID")

save(all_clustering_res, file = "3_data_analysis/02_control_data/05_benchmarking/comparison_result/all_clustering_res.rda")

