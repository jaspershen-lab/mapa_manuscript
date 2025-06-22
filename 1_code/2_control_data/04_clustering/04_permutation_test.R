library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# -----------------------------------------------------------------------------
# Permutation Test for Adjusted Rand Index (ARI)
# -----------------------------------------------------------------------------
# This script calculates the statistical significance of clustering results
# by comparing the observed ARI to a distribution of ARIs generated from
# random permutations of the cluster labels.

# 1. Load data and packages ====
library(mclust)

load("3_data_analysis/02_control_data/04_clustering/bc_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/hc_graph_data.rda")
load("3_data_analysis/02_control_data/04_clustering/gn_graph_data.rda")

gn_node_data <- gn_graph_data |> activate(what = "nodes") |> as_tibble()

hc_node_data <- hc_graph_data |> activate(what = "nodes") |> as_tibble()
hc_cluster_re <- hc_node_data |> select(node, hc_result)

bc_node_data <- bc_graph_data |> activate(what = "nodes") |> as_tibble()
bc_cluster_re <- bc_node_data |> select(node, binary_cut_result)

node_data <-
  gn_node_data |>
  dplyr::mutate(true_label = as.numeric(sub(pattern = "Functional_module_", replacement = "", x = expected_module))) |>
  dplyr::left_join(hc_cluster_re, by = "node") |>
  dplyr::left_join(bc_cluster_re, by = "node")

# 2. DEFINE THE PERMUTATION TEST FUNCTION =====
# -----------------------------------------------------------------------------
#' Performs a permutation test to calculate the p-value for an observed ARI.
#'
#' @param ground_truth A numeric or character vector of the true class labels.
#' @param predicted_clusters A numeric or character vector of the predicted cluster labels.
#' @param n_permutations An integer specifying the number of permutations to run.
#'        A value of 9999 is common for generating p-values up to 0.0001.
#'
#' @return A list containing the observed_ari and the calculated p_value.

ari_permutation_test <- function(ground_truth, predicted_clusters, n_permutations = 9999) {

  # --- Step A: Calculate the observed ARI from the actual data ---
  observed_ari <- adjustedRandIndex(ground_truth, predicted_clusters)

  # --- Step B: Generate the null distribution ---
  #Null Hypothesis: The observed clustering is no better than a random assignment of data points to clusters of the same sizes.
  #                 The observed ARI score could have easily been obtained by chance.
  # We create a null distribution of ARI scores by repeatedly shuffling the
  # predicted labels and calculating the ARI against the true labels. This
  # simulates the scores we would get under the null hypothesis (i.e., by random chance).

  # Create an empty vector to store the ARI from each permutation
  permuted_aris <- numeric(n_permutations)

  # Loop to perform permutations
  for (i in 1:n_permutations) {
    # Randomly shuffle (permute) the predicted cluster labels.
    # This breaks the true association between the prediction and the ground truth
    # while keeping the number and size of clusters identical.
    shuffled_labels <- sample(predicted_clusters)

    # Calculate ARI for the permuted data and store it
    permuted_aris[i] <- adjustedRandIndex(ground_truth, shuffled_labels)
  }

  # --- Step C: Calculate the p-value ---
  # The p-value is the probability of observing an ARI as high as or higher than
  # the actual observed ARI, assuming the null hypothesis is true.

  # Count how many times the permuted ARI was greater than or equal to the observed ARI
  count_extreme_or_more <- sum(permuted_aris >= observed_ari)

  # The p-value is (number of "more extreme" results + 1) / (total permutations + 1).
  # We add 1 to the numerator and denominator to account for the observed data point itself.
  p_value <- (count_extreme_or_more + 1) / (n_permutations + 1)

  # Return the results in a named list for clarity
  return(list(
    observed_ari = observed_ari,
    p_value = p_value,
    permutations = n_permutations,
    permuted_scores = permuted_aris
  ))
}

# 3. RUN THE SIGNIFICANCE TEST ON `node_data` =====
# -----------------------------------------------------------------------------
# Set a seed for random number generation to ensure the permutation results are reproducible.
set.seed(123)

# --- Test 1: For `binary_cut_result` (ARI1) ---
ari1_test <- ari_permutation_test(
  ground_truth = node_data$true_label,
  predicted_clusters = node_data$binary_cut_result,
  n_permutations = 9999 # Using 9999 permutations for a robust p-value
)

# --- Test 2: For `gn_result` (ARI2) ---
ari2_test <- ari_permutation_test(
  ground_truth = node_data$true_label,
  predicted_clusters = node_data$gn_result,
  n_permutations = 9999
)


# 4. VIEW AND INTERPRET THE SIGNIFICANCE RESULTS =====
# -----------------------------------------------------------------------------

cat("--- Permutation Test for ARI1 ('binary_cut_result') vs. Random ---\n")
print(ari1_test)

cat("\n--- Permutation Test for ARI2 ('gn_result') vs. Random ---\n")
print(ari2_test)


# 5. COMPARE THE TWO CLUSTERING METHODS ====
# -----------------------------------------------------------------------------
#' Performs a permutation test to determine if one clustering result is
#' statistically significantly better than another. This is for dependent ARIs
#' (i.e., two methods run on the same data).
#'
#' @param ground_truth A vector of the true class labels.
#' @param clusters1 The first set of predicted labels (the one you hypothesize is better).
#' @param clusters2 The second set of predicted labels.
#' @param n_permutations The number of permutations to run.
#' @return A list with the observed ARI difference and the p-value for that difference.

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

# --- Run the comparison test ---
# We hypothesize that binary_cut_result (ARI1) is better than gn_result (ARI2)
comparison_result <- compare_ari_test(
  ground_truth = node_data$true_label,
  clusters1 = node_data$binary_cut_result, # Method 1
  clusters2 = node_data$gn_result,      # Method 2
  n_permutations = 9999
)

# --- View the comparison result ---
cat("\n--- Statistical Comparison of ARI1 vs ARI2 ---\n")
print(comparison_result)


# 6. VISUALIZE THE PERMUTATION TEST RESULTS ====
# -----------------------------------------------------------------------------
# Convert permuted scores to data frames for ggplot
df_ari1 <- data.frame(score = ari1_test$permuted_scores)
df_ari2 <- data.frame(score = ari2_test$permuted_scores)
df_comp <- data.frame(score = comparison_result$permuted_scores)

# --- Plot 1: Significance of ARI1 ('bc' method) ---
# Calculate the 0.001 cutoff (99.9th percentile of the random distribution)
cutoff1 <- quantile(df_ari1$score, 0.999)

plot1 <- ggplot(df_ari1, aes(x = score)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  geom_vline(aes(xintercept = ari1_test$observed_ari, color = "Observed ARI"), linetype = "solid", linewidth = 1) +
  geom_vline(aes(xintercept = cutoff1, color = "p = 0.001 Cutoff"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Legend", values = c("Observed ARI" = "red", "p = 0.001 Cutoff" = "blue")) +
  labs(
    title = "Permutation Test for ARI1 ('bc' Method)",
    subtitle = "Is the result significantly better than random?",
    x = "Adjusted Rand Index (ARI)",
    y = "Frequency"
  ) +
  theme_minimal()

plot1
# --- Plot 2: Significance of ARI2 ('gn' method) ---
cutoff2 <- quantile(df_ari2$score, 0.999)
plot2 <- ggplot(df_ari2, aes(x = score)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  geom_vline(aes(xintercept = ari2_test$observed_ari, color = "Observed ARI"), linetype = "solid", linewidth = 1) +
  geom_vline(aes(xintercept = cutoff2, color = "p = 0.001 Cutoff"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Legend", values = c("Observed ARI" = "red", "p = 0.001 Cutoff" = "blue")) +
  labs(
    title = "Permutation Test for ARI2 ('gn' Method)",
    subtitle = "Is the result significantly better than random?",
    x = "Adjusted Rand Index (ARI)",
    y = "Frequency"
  ) +
  theme_minimal()
plot2
# --- Plot 3: Comparison of ARI1 vs ARI2 ---
# Calculate the 0.05 cutoff (95th percentile of the random distribution)
cutoff_comp <- quantile(df_comp$score, 0.95)

plot3 <- ggplot(df_comp, aes(x = score)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  geom_vline(aes(xintercept = comparison_result$observed_ari_difference, color = "Observed Difference"), linetype = "solid", linewidth = 1) +
  geom_vline(aes(xintercept = cutoff_comp, color = "p = 0.05 Cutoff"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "Legend", values = c("Observed Difference" = "red", "p = 0.05 Cutoff" = "grey")) +
  labs(
    title = "Comparison of Clustering Methods (ARI(bc) - ARI(gn))",
    subtitle = "Is the difference between methods statistically significant?",
    x = "Difference in Adjusted Rand Index (ARI(bc) - ARI(gn))",
    y = "Frequency"
  ) +
  theme_minimal()
plot3

ggsave(plot = plot1, filename = "3_data_analysis/02_control_data/04_clustering/permutation_test/ARI_bc_test.pdf", width = 8, height = 6)
ggsave(plot = plot2, filename = "3_data_analysis/02_control_data/04_clustering/permutation_test/ARI_gn_test.pdf", width = 8, height = 6)
ggsave(plot = plot3, filename = "3_data_analysis/02_control_data/04_clustering/permutation_test/ARI_comparison_bc_vs_gn_test.pdf", width = 8, height = 6)
