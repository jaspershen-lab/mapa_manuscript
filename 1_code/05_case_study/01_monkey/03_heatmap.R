library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/result_with_module.rda")

library(pheatmap)

# mixed_cases <- result_with_module |>
#   # Group by tissue and fm_module
#   group_by(tissue, fm_module) %>%
#   summarise(
#     unique_clusters = list(unique(cluster)),
#     cluster_types = paste(sort(unique(cluster)), collapse = " & "),
#     has_cluster_u = "Cluster U" %in% cluster,
#     has_cluster_d = "Cluster D" %in% cluster,
#     total_entries = n(),
#     .groups = 'drop'
#   ) %>%
#   # Filter for cases where BOTH clusters are present in the same tissue-module combination
#   filter(has_cluster_u & has_cluster_d)

# Step 1: Create signed values based on cluster type
heatmap_data <- result_with_module |>
  mutate(
    # Create signed molecule count: positive for cluster U, negative for cluster D
    signed_count = case_when(
      cluster == "Cluster U" ~ mapped_molecule_count,
      cluster == "Cluster D" ~ -mapped_molecule_count,
      TRUE ~ mapped_molecule_count  # fallback for any other cluster types
    )
  ) |>
  # Step 2: Group by tissue and fm_module, sum the signed counts
  # (in case there are multiple entries for the same tissue-module combination)
  group_by(tissue, fm_module) |>
  summarise(total_signed_count = sum(signed_count, na.rm = TRUE), .groups = 'drop')

# Step 3: Convert to wide format for heatmap matrix
heatmap_matrix <- heatmap_data %>%
  pivot_wider(
    names_from = fm_module,
    values_from = total_signed_count,
    values_fill = NA
  ) %>%
  tibble::column_to_rownames("tissue")  # Set tissue as row names

# Convert to matrix format (required for pheatmap)
heatmap_matrix <- as.matrix(heatmap_matrix)
any(!is.finite(heatmap_matrix))
heatmap_matrix[is.infinite(heatmap_matrix)] <- NA


# Step 4: Create the heatmap
pheatmap::pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,           # Cluster tissues
  cluster_cols = FALSE,           # Cluster modules
  display_numbers = FALSE,
  fontsize = 8,                  # Adjust font size as needed
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45,                # Rotate column labels
  na_col = "grey"
)

summary(as.vector(heatmap_matrix))
anyNA(heatmap_matrix)
any(!is.finite(as.matrix(heatmap_matrix)))   # TRUE if Inf/NaN present

heatmap_plot_2 <-
  pheatmap(heatmap_matrix,
           cluster_rows = TRUE,           # Cluster algorithms by similarity
           cluster_cols = TRUE,          # Keep cluster numbers in order
           display_numbers = FALSE,        # Show ARI values in cells
           number_color = "white",        # Color of the numbers
           fontsize_number = 8,           # Size of numbers in cells
           # main = "ARI Scores by Algorithm and Number of Clusters",
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           right_annotation = right_anno,
           fontsize = 10,
           angle_col = "90",
           border_color = "white",
           border = TRUE
  )

save(heatmap_matrix, file = "2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/heatmap_matrix.rda")

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/heatmap_matrix.rda")
