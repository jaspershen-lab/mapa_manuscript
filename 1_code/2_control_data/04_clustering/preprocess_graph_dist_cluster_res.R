library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Data preprocessing ====
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_graph_based_result/results_0.01.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/kmeans_results_df.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/hdb_result_df.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/hcluster_results.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/gmm_result_df.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/binary_cut_result.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/ap_results_df.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/ms_result_df.rda")

# Make sure colnames include: Algorithm, Num_clusters, ARI score
colnames(results_0.01)

graph_based_result <- results_0.01 |>
  rename(ari_score = ARI) |>
  select(Algorithm, Num_Clusters, ari_score)

colnames(kmeans_results_df)
kmeans_result <- kmeans_results_df |>
  mutate(Algorithm = "kmeans") |>
  select(Algorithm, Num_Clusters, ari_score)

colnames(hdb_result_df)
hdb_result <- hdb_result_df |>
  rename(ari_score = ari_score_with_noise,
         Num_Clusters = NumClusters) |>
  mutate(Algorithm = "HDBSCAN") |>
  select(Algorithm, Num_Clusters, ari_score)

colnames(hcluster_results)
hcluster_result <- hcluster_results |>
  rename(ari_score = ARI,
         Num_Clusters = K) |>
  mutate(Algorithm = paste("hierarchical", Method, sep = "_")) |>
  select(Algorithm, Num_Clusters, ari_score)

colnames(gmm_result_df)
gmm_result <- gmm_result_df |>
  mutate(Algorithm = "GMM") |>
  select(Algorithm, Num_Clusters, ari_score)

colnames(binary_cut_result_df)
binary_cut_result <- binary_cut_result_df |>
  rename(Num_Clusters = num_cluster) |>
  mutate(Algorithm = "Binary_cut") |>
  select(Algorithm, Num_Clusters, ari_score)

colnames(ap_results_df)
ap_result <- ap_results_df |>
  mutate(Algorithm = "Affinity_propagation") |>
  select(Algorithm, Num_Clusters, ari_score)

colnames(ms_result_df)
ms_result <- ms_result_df |>
  mutate(Algorithm = "Mean_shift") |>
  select(Algorithm, Num_Clusters, ari_score)

dist_based_result <- data.frame(
  Algorithm = character(),
  Num_Clusters = integer(),
  ari_score = numeric()
)
dist_based_result <- rbind(dist_based_result,
                           kmeans_result,
                           hdb_result,
                           hcluster_result,
                           gmm_result,
                           binary_cut_result,
                           ap_result,
                           ms_result)
head(dist_based_result)
dist_based_result <- dist_based_result |>
  group_by(Algorithm, Num_Clusters) |>
  slice_max(order_by = ari_score, n = 1, with_ties = FALSE) |>
  ungroup()

graph_based_result <- graph_based_result |>
  group_by(Algorithm, Num_Clusters) |>
  slice_max(order_by = ari_score, n = 1, with_ties = FALSE) |>
  ungroup()

colnames(dist_based_result)
colnames(graph_based_result)

dist_based_result <- dist_based_result |> mutate(Class = "Distance_based")
graph_based_result <- graph_based_result |> mutate(Class = "Graph_based")
all_results <- rbind(dist_based_result, graph_based_result)

save(all_results, file = "3_data_analysis/02_control_data/04_clustering/all_results.rda")

load("3_data_analysis/02_control_data/04_clustering/all_results.rda")
all_results <- all_results |> filter(!is.na(Num_Clusters))

all_results_2_43 <- all_results |>
  filter((Num_Clusters >= 2 & Num_Clusters <= 43)) |>
  filter(!(Algorithm %in% c("fluid_community", "spinglass")))

heatmap_data <- all_results_2_43 |>
  select(Algorithm, Num_Clusters, ari_score) |>
  pivot_wider(names_from = Num_Clusters,
              values_from = ari_score) |>
  column_to_rownames("Algorithm")
save(heatmap_data, file = "3_data_analysis/02_control_data/04_clustering/heatmap_data_all_clustering_ari.rda")

# # Visualization ====
# # Reshape the data from long to wide format for the heatmap
# load("3_data_analysis/02_control_data/04_clustering/all_results.rda")
# all_results <- all_results |> filter(!is.na(Num_Clusters))
#
# all_results_2_43 <- all_results |>
#   filter((Num_Clusters >= 2 & Num_Clusters <= 43)) |>
#   filter(!(Algorithm %in% c("fluid_community", "spinglass")))
#
# heatmap_data <- all_results_2_43 |>
#   select(Algorithm, Num_Clusters, ari_score) |>
#   pivot_wider(names_from = Num_Clusters,
#               values_from = ari_score) |>
#   column_to_rownames("Algorithm")
# save(heatmap_data, file = "3_data_analysis/02_control_data/04_clustering/heatmap_data_all_clustering_ari.rda")
#
# library(pheatmap)
#
# cluster_nums <- sort(as.numeric(colnames(heatmap_data)))
# heatmap_data <- heatmap_data[, as.character(cluster_nums)]
#
# # Convert to matrix (required for pheatmap)
# heatmap_matrix <- as.matrix(heatmap_data)
#
# # Create row annotation for algorithm classes
# annotation_row <- all_results |>
#   select(Algorithm, Class) |>
#   distinct(.keep_all = TRUE) |>
#   column_to_rownames("Algorithm")
#
# # Make sure annotation matches the algorithms in our matrix
# annotation_row <- annotation_row[rownames(heatmap_matrix), , drop = FALSE]
#
# annotation_colors <- list(Class = class_colors)

# # Create the heatmap single plot
# heatmap_plot <-
#   pheatmap(heatmap_matrix,
#            cluster_rows = TRUE,           # Cluster algorithms by similarity
#            cluster_cols = TRUE,          # Keep cluster numbers in order
#            display_numbers = FALSE,        # Show ARI values in cells
#            number_color = "white",        # Color of the numbers
#            fontsize_number = 8,           # Size of numbers in cells
#            main = "ARI Scores by Algorithm and Number of Clusters",
#            annotation_row = annotation_row,
#            annotation_colors = annotation_colors,
#            fontsize = 10,
#            angle_col = "90",
#            border_color = "white",
#            silent = TRUE)
# heatmap_plot
# ggsave(plot = heatmap_plot, filename = "3_data_analysis/02_control_data/04_clustering/heatmap_plot.pdf", width = 8, height = 6)
#
# # Extract the clustered row order from heatmap
# row_order <- heatmap_plot$tree_row$order
# algorithm_order <- rownames(heatmap_matrix)[row_order]
# all_results_2_43$Algorithm <- factor(all_results_2_43$Algorithm,
#                                      levels = rev(algorithm_order))
#
# combined_scatter_box <- all_results_2_43 |>
#   ggplot(aes(x = ari_score, y = Algorithm)) +
#   # Add box plots with transparency
#   geom_boxplot(aes(fill = Class), alpha = 0.5, outlier.shape = NA, width = 0.6) +
#   # Add scatter points on top
#   geom_point(aes(color = Class), size = 2, alpha = 0.8,
#              position = position_jitter(height = 0.2, width = 0)) +
#   scale_fill_manual(values = class_colors) +
#   scale_color_manual(values = class_colors) +
#   # labs(x = "ARI Score", y = "", title = "ARI Score Distribution") +
#   theme_minimal() +
#   theme(
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(size = 8),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     legend.position = "none",
#     # plot.title = element_text(size = 10, hjust = 0.5),
#     panel.grid.minor = element_blank()
#   ) +
#   xlim(range(all_results_2_43$ari_score, na.rm = TRUE))
#
#
# combined_scatter_box
# heatmap_plot
