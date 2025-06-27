library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_matrix.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_matrix.rda")

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)
ground_truth_dt <-
  control_dt |>
  dplyr::mutate(ground_truth_label = as.numeric(sub(pattern = "Functional_module_", replacement = "", x = expected_module)))
ground_truth_dt <- ground_truth_dt[match(rownames(embedding_matrix), ground_truth_dt$id),]
ground_truth <- ground_truth_dt$ground_truth_label
names(ground_truth) <- ground_truth_dt$id
ground_truth

library(mclust)

# k-means =====
# --- Step 1: L2 Normalize similarity matrix ---
norms <- sqrt(apply(embedding_sim_matrix^2, 1, sum)) # all the lengths of vectors are already 1
normalized_sim <- embedding_sim_matrix / norms

# # --- Step 2: Reduce dimensionality on the NORMALIZED data using UMAP ---
# library(uwot)      # For running UMAP
# library(mclust)   # For calculating Adjusted Rand Index (ARI)
#
# component_values <- seq(2, 40, 1)  #strictly less than min(nrow(A), ncol(A))
#
# # Create an empty data frame to store the results
# results_df <- data.frame(n_components = integer(),
#                          ari_score = numeric())
#
# # Loop through each n_components value
# for (n_comp in component_values) {
#
#   cat(paste("--- Testing n_components =", n_comp, "---\n"))
#
#   # Run UMAP on the *normalized* high-dimensional data
#   set.seed(42) # Use the same seed for UMAP in each loop for consistency
#   umap_output <- umap(embedding_matrix,
#                       n_components = n_comp,
#                       metric = "cosine") # Explicitly use cosine metric
#
#   # Run K-Means on the low-dimensional UMAP output
#   set.seed(42) # Use the same seed for k-means for consistency
#   kmeans_result <- kmeans(umap_output,
#                           centers = 12, # Use the known number of clusters
#                           nstart = 50)  # Run multiple times to find a good solution
#
#   cluster_assignments <- kmeans_result$cluster
#
#   # Calculate the Adjusted Rand Index (ARI)
#   current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)
#
#   cat(paste("ARI Score:", round(current_ari, 4), "\n"))
#
#   # Store the results
#   results_df <- rbind(results_df, data.frame(n_components = n_comp,
#                                              ari_score = current_ari))
# }
#
# # --- Step 4: Analyze and Visualize the Results ---
# # Find the best n_components value
# best_result <- results_df[which.max(results_df$ari_score), ]
# cat(paste("\nBest performance found with n_components =", best_result$n_components,
#           "which yielded an ARI of", round(best_result$ari_score, 4), "\n"))

# Plot the results to visualize the performance curve
# tuning_plot <- ggplot(results_df, aes(x = n_components, y = ari_score)) +
#   geom_line(color = "steelblue", linewidth = 1) +
#   geom_point(color = "steelblue", size = 3) +
#   geom_point(data = best_result, aes(x = n_components, y = ari_score), color = "red", size = 5) +
#   geom_text(data = best_result,
#             aes(label = paste("Best ARI:", round(ari_score, 3))),
#             vjust = -1.5, color = "red") +
#   labs(title = "UMAP n_components vs. Clustering Performance (ARI)",
#        x = "Number of UMAP Components",
#        y = "Adjusted Rand Index (ARI)") +
#   theme_minimal(base_size = 14) +
#   scale_x_continuous(breaks = component_values)
#
# tuning_plot
# ggsave(plot = tuning_plot, filename = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/kmeans_umap_tuning_plot.pdf", width = 8, height = 6)

# set.seed(42)
# umap_output <- umap(embedding_matrix,
#                     n_components = best_result$n_components,
#                     metric = "cosine") # Explicitly use cosine metric

kmeans_results_df <- data.frame(n_centers = integer(),
                                ari_score = numeric(),
                                Num_Clusters = integer())

for (i in 2:43) {
  set.seed(23)
  kmeans_result <- kmeans(normalized_sim,
                          centers = i, # Use the known number of clusters
                          nstart = 50)  # Run multiple times to find a good solution (the lowest total within-cluster sum of squares)
  cluster_assignments <- kmeans_result$cluster

  # Calculate the Adjusted Rand Index (ARI)
  current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)
  # Store the results
  kmeans_results_df <- rbind(kmeans_results_df, data.frame(n_centers = i,
                                                           ari_score = current_ari,
                                                           Num_Clusters = length(unique(cluster_assignments))))
}

best_ncenters <- kmeans_results_df[which.max(kmeans_results_df$ari_score), ]
best_ncenters$ari_score
set.seed(12)
kmeans_result <- kmeans(normalized_sim,
                        centers = best_ncenters$n_centers,
                        nstart = 50)
kmeans_cluster_assignments <- kmeans_result$cluster
save(kmeans_results_df, file = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/kmeans_results_df.rda")

# hierarchical clustering =====
## Convert similarity to distance
cosine_dist_matrix <- 1 - embedding_sim_matrix

## Coerce into a 'dist' object for R's clustering functions
cosine_dist_object <- as.dist(cosine_dist_matrix)

# Define all available methods in hclust
hclust_methods <- c("ward.D", "ward.D2", "single", "complete",
                    "average", "mcquitty", "median", "centroid")

# Define range of clusters to test
k_values <- seq(2, 43, 1)

# Initialize results dataframe
hcluster_results <- data.frame()

## calculate ARI for given method and k
calculate_ari <- function(dist_obj, method, k, true_labels) {
  tryCatch({
    # Perform hierarchical clustering
    hc_result <- hclust(dist_obj, method = method)

    # Cut tree to get k clusters
    clusters <- cutree(hc_result, k = k)

    # Calculate ARI
    ari <- adjustedRandIndex(clusters, true_labels)

    return(ari)
  }, error = function(e) {
    warning(paste("Error with method", method, "and k =", k, ":", e$message))
    return(NA)
  })
}

for (method in hclust_methods) {
  cat(paste("\nTesting method:", method, "\n"))

  for (k in k_values) {
    # Calculate ARI
    ari_value <- calculate_ari(cosine_dist_object, method, k, ground_truth)

    # Store results
    hcluster_results <- rbind(hcluster_results,
                              data.frame(Method = method,
                                         K = k,
                                         ARI = ari_value))
  }
}

best_combination_method_k <- hcluster_results[which.max(hcluster_results$ARI), ]
best_combination_method_k

save(hcluster_results, file = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/hcluster_results.rda")


# Binary cut clustering ====
# select the best cutoff
library(simplifyEnrichment)
result <- select_cutoff(mat = embedding_sim_matrix,
                        cutoff = seq(0.5, 0.95, 0.01),
                        try_all_partition_fun = TRUE)

# k = 13
# cut.off = 0.69
select_cutoff = function(mat, cutoff = seq(0.6, 0.98, by = 0.01), verbose = se_opt$verbose, ...) {

  cutoff = cutoff[cutoff >= 0.5 & cutoff <= 1]
  if(length(cutoff) == 0) {
    stop("`cutoff` should be within [0.5, 1].")
  }

  s1 = s2 = s3 = s4 = ari_score = NULL
  for(i in seq_along(cutoff)) {
    # if(verbose) message(GetoptLong::qq("@{i}/@{length(cutoff)}, cutoff = @{cutoff[i]}..."))
    set.seed(123)
    cl = binary_cut(mat, cutoff = cutoff[i], try_all_partition_fun = TRUE)
    ari_score[i] = adjustedRandIndex(cl, ground_truth)
    s1[i] = difference_score(mat, cl)
    tb = table(cl)
    s2[i] = length(tb)
    s3[i] = sum(tb >= 5)
    s4[i] = block_mean(mat, cl)
  }

  check_pkg("cowplot", bioc = FALSE)
  check_pkg("ggplot2", bioc = FALSE)

  suppressWarnings(
    p1 <- ggplot2::ggplot(data = NULL, ggplot2::aes(x = cutoff, y = s1)) +
      ggplot2::geom_point() + ggplot2::ylab("Difference score") +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
  )

  df1 = data.frame(cutoff, s2); colnames(df1) = c("method", "value")
  df2 = data.frame(cutoff, s3); colnames(df2) = c("method", "value")
  df1$type = "All sizes"
  df2$type = "Size >= 5"
  df = rbind(df1, df2)
  suppressWarnings(
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = df$method, y = df$value, col = df$type, fill = df$type)) +
      ggplot2::geom_point() + ggplot2::ylab("Number of clusters") + ggplot2::labs(col = "Type", fill = "Type") +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
  )

  suppressWarnings(
    p3 <- ggplot2::ggplot(data = NULL, ggplot2::aes(x = cutoff, y = s4)) +
      ggplot2::geom_point() + ggplot2::ylab("Block mean") + ggplot2::xlab("Cutoff") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  )

  suppressWarnings(print(cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v", axis = "lr", rel_heights = c(1, 1, 1.3))))
}

block_mean = function(mat, cl) {
  n = nrow(mat)
  l_block = matrix(FALSE, nrow = nrow(mat), ncol = ncol(mat))

  for(le in unique(cl)) {
    l = cl == le
    l_block[l, l] = TRUE
  }
  l_block2 = l_block
  l_block2[upper.tri(mat, diag = TRUE)] = FALSE
  x1 = mat[l_block2]
  mean(x1)
}

binary_cut_result_df <- data.frame(cutoff = cutoff, difference = s1, num_cluster = s2, block_mean = s4, ari_score = ari_score)
save(binary_cut_result_df, file = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/binary_cut_result.rda")
## best cutoff
binary_cut_result <- binary_cut(mat = embedding_sim_matrix,
                                cutoff = 0.69)


# DBSCAN ======
# library(dbscan)
#
# ## Convert similarity to distance - DBSCAN needs a distance matrix.
# cosine_dist_matrix <- 1 - embedding_sim_matrix
# cosine_dist_object <- as.dist(cosine_dist_matrix)
#
# # Set MinPts
# min_pts_value <- 7
# # Create the k-distance plot using k = MinPts - 1
# kNNdistplot(cosine_dist_object, k = min_pts_value - 1)
# abline(h = 0.59, lty = 2) # Look for the "knee" in the plot
#
# # Perform DBSCAN clustering
# db_result <- dbscan(cosine_dist_object, eps = 0.59, minPts = min_pts_value)
#
# # Print the results (cluster 0 represents noise points)
# print(db_result)

# HDBSCAN =====
library(dbscan)

cosine_dist_matrix <- 1 - embedding_sim_matrix
cosine_dist_object <- as.dist(cosine_dist_matrix)

hdb_result_df <- data.frame(minpts = integer(),
                            ari_without_noise = numeric(),
                            ari_score_with_noise = numeric(),
                            noise_percentage = numeric(),
                            NumClusters = integer())

for (i in 2:43) {
  hdb_result <- hdbscan(cosine_dist_object, minPts = i)
  cluster_assignments <- hdb_result$cluster

  # We find the indices of all points that are NOT noise (i.e., cluster != 0)
  non_noise_indices <- which(cluster_assignments != 0)

  filtered_assignments <- cluster_assignments[non_noise_indices]
  filtered_ground_truth <- ground_truth[non_noise_indices]

  noise_percentage <- (length(cluster_assignments) - length(filtered_assignments)) / length(cluster_assignments) * 100

  # Calculate the Adjusted Rand Index (ARI)
  ari_without_noise <- adjustedRandIndex(filtered_assignments, filtered_ground_truth)
  ari_score_with_noise <- adjustedRandIndex(cluster_assignments, ground_truth)
  # Store the results
  hdb_result_df <- rbind(hdb_result_df, data.frame(minpts = i,
                                                   ari_without_noise = ari_without_noise,
                                                   ari_score_with_noise = ari_score_with_noise,
                                                   noise_percentage = noise_percentage,
                                                   NumClusters = sum(hdb_result$cluster == 0) + length(unique(cluster_assignments)) - 1
                                                   )
                         )
}

save(hdb_result_df, file = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/hdb_result_df.rda")

# best_minpts
best_min_pts <- hdb_result_df[which.max(hdb_result_df$ari_without_noise),]
best_min_pts$minpts

# Affinity Propagation ====
library(apcluster)

# --- Run Affinity Propagation ---
p_values <- seq(0.2, 0.9, 0.01)

ap_results_df <- data.frame(p = integer(),
                            ari_score = numeric(),
                            Num_Clusters = integer())

for (p in p_values) {
  set.seed(23)
  ap_result <- apcluster(s = embedding_sim_matrix, p = p)

  cluster_assignments <- as.integer(labels(ap_result, type="enum"))

  # Calculate the Adjusted Rand Index (ARI)
  current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)
  # Store the results
  ap_results_df <- rbind(ap_results_df, data.frame(p = p,
                                                   ari_score = current_ari,
                                                   Num_Clusters = length(ap_result@clusters)))
}

# Find the best p value
best_p <- ap_results_df[which.max(ap_results_df$ari_score), ]
ap_result <- apcluster(s = embedding_sim_matrix, p = best_p$p)
# To get a simple integer vector of cluster assignments (like with kmeans/dbscan),
cluster_assignments <- as.integer(labels(ap_result, type="enum"))
table(cluster_assignments)

# You can see which points are the exemplars
exemplars <- ap_result@exemplars
print(exemplars)

save(ap_results_df, file = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/ap_results_df.rda")

# Mean-Shift ====
library(LPCM)
bandwidths <- seq(0.05, 1, 0.01)

ms_result_df <- data.frame(bandwidth = numeric(),
                           ari_score = numeric(),
                           Num_Clusters = integer())

for (h in bandwidths) {
  set.seed(42)
  ms_result <- ms(embedding_sim_matrix, h = h)

  # Get the cluster assignments
  cluster_assignments <- ms_result$cluster.label

  # Calculate the Adjusted Rand Index (ARI)
  current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)

  # Store the results
  ms_result_df <- rbind(ms_result_df, data.frame(bandwidth = h,
                                                  ari_score = current_ari,
                                                  Num_Clusters = length(unique(cluster_assignments))))
}

save(ms_result_df, file = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/ms_result_df.rda")

# Gaussian Mixtures ====
library(mclust)

Gs <- seq(2, 43, 1)

gmm_result_df <- data.frame(G = integer(),
                            ari_score = numeric(),
                            Num_Clusters = integer())
for (g in Gs) {
  set.seed(42)
  gmm_result <- Mclust(embedding_sim_matrix, G = g)

  # Get the cluster assignments
  cluster_assignments <- gmm_result$classification

  # Calculate the Adjusted Rand Index (ARI)
  current_ari <- adjustedRandIndex(cluster_assignments, ground_truth)

  # Store the results
  gmm_result_df <- rbind(gmm_result_df, data.frame(G = g,
                                                    ari_score = current_ari,
                                                    Num_Clusters = length(unique(cluster_assignments))))
}

save(gmm_result_df, file = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_dist_based_result/gmm_result_df.rda")
