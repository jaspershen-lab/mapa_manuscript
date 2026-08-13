## Sensitivity analysis on the informative GO set benchmark (k = 20)
## Step 5: ARI heatmap (algorithm x number of clusters) across all scanned
## clustering methods, in the style of Fig. 3a of the control-data analysis.
## Adapted from 1_code/2_control_data/04_clustering/preprocess_graph_dist_cluster_res.R
## and 1_code/2_control_data/04_clustering/02_fig3a_cluster_res_heatmap.R
## (cluster number range 2-150 instead of 2-43; 313 terms / 57 true modules)

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

res_dir <- "3_data_analysis/10_sensitivity_res_GO"

# Data preprocessing ====
load(file.path(res_dir, "graph_based_results.rda"))
load(file.path(res_dir, "kmeans_results_df.rda"))
load(file.path(res_dir, "hdb_result_df.rda"))
load(file.path(res_dir, "hcluster_results.rda"))
load(file.path(res_dir, "gmm_result_df.rda"))
load(file.path(res_dir, "binary_cut_result_df.rda"))
load(file.path(res_dir, "ap_results_df.rda"))
load(file.path(res_dir, "ms_result_df.rda"))

graph_based_result <- graph_based_results |>
  rename(ari_score = ARI) |>
  select(Algorithm, Num_Clusters, ari_score)

kmeans_result <- kmeans_results_df |>
  mutate(Algorithm = "kmeans") |>
  select(Algorithm, Num_Clusters, ari_score)

hdb_result <- hdb_result_df |>
  rename(ari_score = ari_score_with_noise,
         Num_Clusters = NumClusters) |>
  mutate(Algorithm = "HDBSCAN") |>
  select(Algorithm, Num_Clusters, ari_score)

hcluster_result <- hcluster_results |>
  rename(ari_score = ARI,
         Num_Clusters = K) |>
  mutate(Algorithm = paste("hierarchical", Method, sep = "_")) |>
  select(Algorithm, Num_Clusters, ari_score)

gmm_result <- gmm_result_df |>
  mutate(Algorithm = "GMM") |>
  select(Algorithm, Num_Clusters, ari_score)

binary_cut_result <- binary_cut_result_df |>
  mutate(Algorithm = "Binary_cut") |>
  select(Algorithm, Num_Clusters, ari_score)

ap_result <- ap_results_df |>
  mutate(Algorithm = "Affinity_propagation") |>
  select(Algorithm, Num_Clusters, ari_score)

ms_result <- ms_result_df |>
  mutate(Algorithm = "Mean_shift") |>
  select(Algorithm, Num_Clusters, ari_score)

dist_based_result <- rbind(kmeans_result,
                           hdb_result,
                           hcluster_result,
                           gmm_result,
                           binary_cut_result,
                           ap_result,
                           ms_result)

## Best ARI per algorithm and cluster number
dist_based_result <- dist_based_result |>
  group_by(Algorithm, Num_Clusters) |>
  slice_max(order_by = ari_score, n = 1, with_ties = FALSE) |>
  ungroup()

graph_based_result <- graph_based_result |>
  group_by(Algorithm, Num_Clusters) |>
  slice_max(order_by = ari_score, n = 1, with_ties = FALSE) |>
  ungroup()

dist_based_result <- dist_based_result |> mutate(Class = "Distance_based")
graph_based_result <- graph_based_result |> mutate(Class = "Graph_based")
all_results <- rbind(dist_based_result, graph_based_result)

save(all_results, file = file.path(res_dir, "all_results.rda"))

all_results <- all_results |> filter(!is.na(Num_Clusters))

all_results_2_150 <- all_results |>
  filter((Num_Clusters >= 2 & Num_Clusters <= 150)) |>
  filter(!(Algorithm %in% c("fluid_community", "spinglass")))

heatmap_data <- all_results_2_150 |>
  select(Algorithm, Num_Clusters, ari_score) |>
  pivot_wider(names_from = Num_Clusters,
              values_from = ari_score) |>
  column_to_rownames("Algorithm")
save(heatmap_data,
     file = file.path(res_dir, "heatmap_data_all_clustering_ari.rda"))

# Visualization ====
library(pheatmap)

cluster_nums <- sort(as.numeric(colnames(heatmap_data)))
heatmap_data <- heatmap_data[, as.character(cluster_nums)]

heatmap_matrix <- as.matrix(heatmap_data)

# Row annotation for algorithm classes
annotation_row <- all_results |>
  select(Algorithm, Class) |>
  distinct(.keep_all = TRUE) |>
  column_to_rownames("Algorithm")

annotation_row <- annotation_row[rownames(heatmap_matrix), , drop = FALSE]

annotation_colors <- list(Class = class_colors)

library(ComplexHeatmap)
box_colors <- class_colors[annotation_row$Class]
right_anno <- rowAnnotation(
  "ARI Score" = anno_boxplot(heatmap_matrix,
                             gp = gpar(col = box_colors)),
  width = unit(3, "cm")
)

## With cluster numbers 2-150 the matrix is sparse (each algorithm only
## reaches part of the range), so compute the row clustering with an
## NA-tolerant distance and keep columns in numeric order
row_dist <- dist(heatmap_matrix)   # stats::dist excludes NAs pairwise
row_dist[is.na(row_dist)] <- max(row_dist, na.rm = TRUE)
row_hclust <- hclust(row_dist)

## With 149 columns, only label every 10th cluster number for readability
col_labels <- ifelse(cluster_nums %% 10 == 0 | cluster_nums == 2,
                     cluster_nums, "")

heatmap_plot_2 <-
  pheatmap(heatmap_matrix,
           cluster_rows = row_hclust,     # Cluster algorithms by similarity
           cluster_cols = FALSE,
           labels_col = col_labels,
           display_numbers = FALSE,
           number_color = "white",
           fontsize_number = 8,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           right_annotation = right_anno,
           fontsize = 10,
           angle_col = "90",
           border_color = "white",
           border = TRUE
  )

heatmap_plot_2

library(Cairo)
CairoPDF(file.path(res_dir, "heatmap_plot_ARI_all_algorithms.pdf"),
         width = 14, height = 6)
heatmap_plot_2
dev.off()

# Significance analysis ====
kruskal_result <- kruskal.test(ari_score ~ Algorithm, data = all_results_2_150)
print(kruskal_result)

library(FSA)
dunn_test <- dunnTest(ari_score ~ Algorithm,
                      data = all_results_2_150,
                      method = "bonferroni")
save(dunn_test,
     file = file.path(res_dir, "significance_analysis_dunn_test_all_cluster_ari.rda"))

library(rcompanion)
cld_result <- cldList(P.adj ~ Comparison,
                      data = dunn_test$res,
                      threshold = 0.05)
print(cld_result)
save(cld_result,
     file = file.path(res_dir, "significance_analysis_cld_letters_all_cluster_ari.rda"))
