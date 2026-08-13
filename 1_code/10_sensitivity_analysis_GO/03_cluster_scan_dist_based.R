## Sensitivity analysis on the informative GO set benchmark (k = 20)
## Step 3: Scan distance/similarity-based clustering algorithms and their
## hyper-parameters, evaluated against the ground truth modules.
## Adapted from 1_code/2_control_data/04_clustering/01_cluster_pathways_dist_based.R
## (k / parameter ranges widened for 313 terms; ground truth has 57 modules)

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(mclust)

load("3_data_analysis/10_sensitivity_res_GO/embedding_sim_matrix.rda")

ground_truth_dt <- readr::read_csv("2_data/informative_go_set/k_20/ground_truth.csv",
                                   show_col_types = FALSE)
ground_truth_dt <-
  ground_truth_dt |>
  dplyr::mutate(ground_truth_label = as.numeric(factor(module_id)))
ground_truth_dt <- ground_truth_dt[match(rownames(embedding_sim_matrix),
                                         ground_truth_dt$go_id), ]
stopifnot(!anyNA(ground_truth_dt$ground_truth_label))
ground_truth <- ground_truth_dt$ground_truth_label
names(ground_truth) <- ground_truth_dt$go_id

k_values <- seq(2, 150, 1)

# k-means =====
## L2 normalize rows of the similarity matrix
norms <- sqrt(apply(embedding_sim_matrix^2, 1, sum))
normalized_sim <- embedding_sim_matrix / norms

kmeans_results_df <- data.frame(n_centers = integer(),
                                ari_score = numeric(),
                                Num_Clusters = integer())

for (i in k_values) {
  set.seed(23)
  kmeans_result <- kmeans(normalized_sim,
                          centers = i,
                          nstart = 50)
  cluster_assignments <- kmeans_result$cluster

  current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)
  kmeans_results_df <- rbind(kmeans_results_df,
                             data.frame(n_centers = i,
                                        ari_score = current_ari,
                                        Num_Clusters = length(unique(cluster_assignments))))
}

save(kmeans_results_df,
     file = "3_data_analysis/10_sensitivity_res_GO/kmeans_results_df.rda")
cat("kmeans best:\n")
print(kmeans_results_df[which.max(kmeans_results_df$ari_score), ])

# Hierarchical clustering =====
cosine_dist_matrix <- 1 - embedding_sim_matrix
cosine_dist_object <- as.dist(cosine_dist_matrix)

hclust_methods <- c("ward.D", "ward.D2", "single", "complete",
                    "average", "mcquitty", "median", "centroid")

calculate_ari <- function(dist_obj, method, k, true_labels) {
  tryCatch({
    hc_result <- hclust(dist_obj, method = method)
    clusters <- cutree(hc_result, k = k)
    adjustedRandIndex(clusters, true_labels)
  }, error = function(e) {
    warning(paste("Error with method", method, "and k =", k, ":", e$message))
    NA
  })
}

hcluster_results <- data.frame()
for (method in hclust_methods) {
  cat(paste("\nTesting method:", method, "\n"))
  for (k in k_values) {
    ari_value <- calculate_ari(cosine_dist_object, method, k, ground_truth)
    hcluster_results <- rbind(hcluster_results,
                              data.frame(Method = method,
                                         K = k,
                                         ARI = ari_value))
  }
}

save(hcluster_results,
     file = "3_data_analysis/10_sensitivity_res_GO/hcluster_results.rda")
cat("hclust best:\n")
print(hcluster_results[which.max(hcluster_results$ARI), ])

# Binary cut clustering ====
library(simplifyEnrichment)

bc_cutoffs <- seq(0.5, 0.95, 0.01)
binary_cut_result_df <- data.frame(cutoff = numeric(),
                                   ari_score = numeric(),
                                   Num_Clusters = integer())

for (co in bc_cutoffs) {
  set.seed(123)
  cl <- tryCatch(
    simplifyEnrichment::binary_cut(mat = embedding_sim_matrix,
                                   cutoff = co,
                                   try_all_partition_fun = TRUE),
    error = function(e) NULL
  )
  if (is.null(cl)) next
  binary_cut_result_df <- rbind(binary_cut_result_df,
                                data.frame(cutoff = co,
                                           ari_score = adjustedRandIndex(cl, ground_truth),
                                           Num_Clusters = length(unique(cl))))
}

save(binary_cut_result_df,
     file = "3_data_analysis/10_sensitivity_res_GO/binary_cut_result_df.rda")
cat("binary cut best:\n")
print(binary_cut_result_df[which.max(binary_cut_result_df$ari_score), ])

# HDBSCAN =====
library(dbscan)

hdb_result_df <- data.frame(minpts = integer(),
                            ari_without_noise = numeric(),
                            ari_score_with_noise = numeric(),
                            noise_percentage = numeric(),
                            NumClusters = integer())

for (i in 2:100) {
  hdb_result <- tryCatch(hdbscan(cosine_dist_object, minPts = i),
                         error = function(e) NULL)
  if (is.null(hdb_result)) next
  cluster_assignments <- hdb_result$cluster

  non_noise_indices <- which(cluster_assignments != 0)
  filtered_assignments <- cluster_assignments[non_noise_indices]
  filtered_ground_truth <- ground_truth[non_noise_indices]

  noise_percentage <- (length(cluster_assignments) - length(filtered_assignments)) /
    length(cluster_assignments) * 100

  ari_without_noise <- adjustedRandIndex(filtered_assignments, filtered_ground_truth)
  ari_score_with_noise <- adjustedRandIndex(cluster_assignments, ground_truth)
  hdb_result_df <- rbind(hdb_result_df,
                         data.frame(minpts = i,
                                    ari_without_noise = ari_without_noise,
                                    ari_score_with_noise = ari_score_with_noise,
                                    noise_percentage = noise_percentage,
                                    NumClusters = sum(cluster_assignments == 0) +
                                      length(unique(cluster_assignments[non_noise_indices]))))
}

save(hdb_result_df,
     file = "3_data_analysis/10_sensitivity_res_GO/hdb_result_df.rda")
cat("HDBSCAN best (with noise):\n")
print(hdb_result_df[which.max(hdb_result_df$ari_score_with_noise), ])

# Affinity Propagation ====
library(apcluster)

p_values <- seq(0.2, 0.9, 0.01)

ap_results_df <- data.frame(p = numeric(),
                            ari_score = numeric(),
                            Num_Clusters = integer())

for (p in p_values) {
  set.seed(23)
  ap_result <- tryCatch(apcluster(s = embedding_sim_matrix, p = p),
                        error = function(e) NULL)
  if (is.null(ap_result) || length(ap_result@clusters) == 0) next

  cluster_assignments <- as.integer(labels(ap_result, type = "enum"))

  current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)
  ap_results_df <- rbind(ap_results_df,
                         data.frame(p = p,
                                    ari_score = current_ari,
                                    Num_Clusters = length(ap_result@clusters)))
}

save(ap_results_df,
     file = "3_data_analysis/10_sensitivity_res_GO/ap_results_df.rda")
cat("Affinity propagation best:\n")
print(ap_results_df[which.max(ap_results_df$ari_score), ])

# Mean-Shift ====
library(LPCM)
bandwidths <- seq(0.05, 1, 0.01)

ms_result_df <- data.frame(bandwidth = numeric(),
                           ari_score = numeric(),
                           Num_Clusters = integer())

for (h in bandwidths) {
  set.seed(42)
  ms_result <- tryCatch(ms(embedding_sim_matrix, h = h, plot = FALSE),
                        error = function(e) NULL)
  if (is.null(ms_result)) next

  cluster_assignments <- ms_result$cluster.label

  current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)
  ms_result_df <- rbind(ms_result_df,
                        data.frame(bandwidth = h,
                                   ari_score = current_ari,
                                   Num_Clusters = length(unique(cluster_assignments))))
}

save(ms_result_df,
     file = "3_data_analysis/10_sensitivity_res_GO/ms_result_df.rda")
cat("Mean-shift best:\n")
print(ms_result_df[which.max(ms_result_df$ari_score), ])

# Gaussian Mixtures ====
Gs <- seq(2, 100, 1)

gmm_result_df <- data.frame(G = integer(),
                            ari_score = numeric(),
                            Num_Clusters = integer())
for (g in Gs) {
  set.seed(42)
  gmm_result <- tryCatch(Mclust(embedding_sim_matrix, G = g, verbose = FALSE),
                         error = function(e) NULL)
  if (is.null(gmm_result)) next

  cluster_assignments <- gmm_result$classification

  current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)
  gmm_result_df <- rbind(gmm_result_df,
                         data.frame(G = g,
                                    ari_score = current_ari,
                                    Num_Clusters = length(unique(cluster_assignments))))
}

save(gmm_result_df,
     file = "3_data_analysis/10_sensitivity_res_GO/gmm_result_df.rda")
if (nrow(gmm_result_df) > 0) {
  cat("GMM best:\n")
  print(gmm_result_df[which.max(gmm_result_df$ari_score), ])
}
