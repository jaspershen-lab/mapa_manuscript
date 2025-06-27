library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(igraph)
library(tidygraph)

# Permutation test for hierarchical_ward.D2 and louvain
# load("3_data_analysis/02_control_data/04_clustering/all_results.rda")
# best_params <- all_results |>
#   filter(Algorithm %in% c("hierarchical_ward.D2", "louvain"))
# best_result <- best_params |> group_by(Algorithm) |> slice_max(order_by = ari_score) |> ungroup()
#
# load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_graph_based_result/results_0.01.rda")
# best_louvain_params <- results_0.01 |> filter(Algorithm == "louvain" & Num_Clusters == 14)
# best cutoff = 0.55 for louvain
# best num_clusters = 8 for hierarchical_ward.D2

# louvain clustering
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_matrix.rda")
control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

setwd("3_data_analysis/02_control_data/04_clustering")

edge_data_with_cutoff <-
  embedding_sim_df |>
  rename(weight = sim) |>
  dplyr::filter(weight > 0.55)

node_data <-
  control_dt |>
  dplyr::rename(node = id) |>
  dplyr::select(node, expected_module, expected_count, database, name)

graph_data_with_cutoff <-
  tbl_graph(
    nodes = node_data,
    edges = edge_data_with_cutoff,
    directed = FALSE,
    node_key = "node"
  ) |>
  mutate(degree = centrality_degree())

louvain_res <- tryCatch(cluster_louvain(graph_data_with_cutoff), error = function(e) NA)
membership(louvain_res)
louvain_res$membership
louvain_graph_data <-
  graph_data_with_cutoff |>
  igraph::upgrade_graph() |>
  tidygraph::activate(what = "nodes") |>
  dplyr::mutate(louvain_result = louvain_res$membership)
# save(louvain_graph_data, file = "louvain_graph_data.rda")

# hierarchical clustering with ward.D2 linkage method
cosine_dist_matrix <- 1 - embedding_sim_matrix
cosine_dist_object <- as.dist(cosine_dist_matrix)
hc_result <- hclust(cosine_dist_object, method = "ward.D2")
clusters <- cutree(hc_result, k = 8)
clusters

edge_data <- embedding_sim_df
node_data <-
  control_dt %>%
  # dplyr::filter(id %in% edge_data$from | id %in% edge_data$to) %>%
  dplyr::rename(node = id) %>%
  dplyr::select(node, expected_module, expected_count, database, name)

hc_graph_data <-
  tbl_graph(
    nodes = node_data,
    edges = edge_data,
    directed = FALSE,
    node_key = "node"
  ) |>
  mutate(degree = centrality_degree())

hc_cluster_res <- data.frame(node = names(clusters),
                             hc_result = unname(clusters))
hc_graph_data <-
  hc_graph_data %>%
  igraph::upgrade_graph() %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::left_join(hc_cluster_res, by = "node")

# permutation test
ground_truth_dt <-
  control_dt |>
  dplyr::mutate(ground_truth_label = as.numeric(sub(pattern = "Functional_module_", replacement = "", x = expected_module)))

louvain_node_data <- louvain_graph_data |> activate(what = "nodes") |> as_tibble()
louvain_cluster_re <- louvain_node_data |> select(node, louvain_result)

node_data <-
  ground_truth_dt |>
  rename(node = id) |>
  dplyr::left_join(hc_cluster_res, by = "node") |>
  dplyr::left_join(louvain_cluster_re, by = "node")

library(mclust)

compare_ari_test <- function(ground_truth, clusters1, clusters2, n_permutations = 9999) {
  # Calculate the observed ARIs for both methods and their actual difference.
  ari1 <- adjustedRandIndex(ground_truth, clusters1)
  ari2 <- adjustedRandIndex(ground_truth, clusters2)
  observed_diff <- ari1 - ari2

  n_points <- length(ground_truth)
  permuted_diffs <- numeric(n_permutations)

  # The Null Hypothesis (H0) is that there is no difference between the two methods.
  # Under H0, for any data point, its label from method 1 and method 2 are interchangeable.
  # We simulate this by randomly swapping the labels for each data point.
  for (i in 1:n_permutations) {
    # Create empty vectors for the new, permuted label sets
    permuted_labels1 <- integer(n_points)
    permuted_labels2 <- integer(n_points)

    # For each data point, randomly decide whether to swap its labels between the two methods.
    # rbinom creates a vector of 0s and 1s (like a coin flip for each point).
    swap_indices <- rbinom(n_points, 1, 0.5) == 1

    # Where swap_indices is FALSE, keep the original label assignment
    permuted_labels1[!swap_indices] <- clusters1[!swap_indices]
    permuted_labels2[!swap_indices] <- clusters2[!swap_indices]

    # Where swap_indices is TRUE, swap the labels between the two methods
    permuted_labels1[swap_indices] <- clusters2[swap_indices]
    permuted_labels2[swap_indices] <- clusters1[swap_indices]

    # Calculate the ARI scores for these new, jumbled label sets
    perm_ari1 <- adjustedRandIndex(ground_truth, permuted_labels1)
    perm_ari2 <- adjustedRandIndex(ground_truth, permuted_labels2)

    # Store the difference in ARI from this permutation
    permuted_diffs[i] <- perm_ari1 - perm_ari2
  }

  # Calculate the one-sided p-value.
  # This is the probability of seeing a difference as large or larger than the one we observed.
  count_extreme <- sum(permuted_diffs >= observed_diff)
  p_value <- (count_extreme + 1) / (n_permutations + 1)

  return(list(
    method1_ari = ari1,
    method2_ari = ari2,
    observed_ari_difference = observed_diff,
    p_value = p_value,
    permuted_scores = permuted_diffs
  ))
}

comparison_result <- compare_ari_test(
  ground_truth = node_data$ground_truth_label,
  clusters2 = node_data$hc_result,  # hierarchical_ward.D2
  clusters1 = node_data$louvain_result,      # Method 2
  n_permutations = 9999
)
df_comp <- data.frame(score = comparison_result$permuted_scores)
cutoff_comp <- quantile(df_comp$score, 0.95)

plot3 <- ggplot(df_comp, aes(x = score)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  geom_vline(aes(xintercept = comparison_result$observed_ari_difference, color = "Observed Difference"), linetype = "solid", linewidth = 1) +
  geom_vline(aes(xintercept = cutoff_comp, color = "p = 0.05 Cutoff"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Legend", values = c("Observed Difference" = "red", "p = 0.05 Cutoff" = "grey")) +
  scale_x_continuous(breaks = seq(-0.25, 0.25, by = 0.05)) +
  labs(
    x = "Difference in Adjusted Rand Index (ARI(louvain) - ARI(ward.D2))",
    y = "Frequency"
  ) +
  theme_bw()
plot3
ggsave(plot = plot3,
       filename = "3_data_analysis/02_control_data/04_clustering/permutation_test_louvain_hc_ward.pdf",
       width = 8, height = 6)
