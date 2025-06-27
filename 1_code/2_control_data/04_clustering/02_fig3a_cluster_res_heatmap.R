library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/02_control_data/04_clustering/heatmap_data_all_clustering_ari.rda")
load("3_data_analysis/02_control_data/04_clustering/all_results.rda")

setwd("3_data_analysis/02_control_data/04_clustering/")

library(pheatmap)

cluster_nums <- sort(as.numeric(colnames(heatmap_data)))
heatmap_data <- heatmap_data[, as.character(cluster_nums)]

# Convert to matrix (required for pheatmap)
heatmap_matrix <- as.matrix(heatmap_data)

# Create row annotation for algorithm classes
annotation_row <- all_results |>
  select(Algorithm, Class) |>
  distinct(.keep_all = TRUE) |>
  column_to_rownames("Algorithm")

# Make sure annotation matches the algorithms in our matrix
annotation_row <- annotation_row[rownames(heatmap_matrix), , drop = FALSE]

annotation_colors <- list(Class = class_colors)

## combine plot
library(ComplexHeatmap)
# After attach ComplexHeatmap, this will generate different object with before
box_colors <- class_colors[annotation_row$Class]
right_anno <- rowAnnotation(
  "ARI Score" = anno_boxplot(heatmap_matrix,
                             gp = gpar(col = box_colors)),
  width = unit(3, "cm")
)
right_anno

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

heatmap_plot_2

library(Cairo)
CairoPDF("heatmap_plot_ARI_all_algorithms.pdf", width = 10, height = 6)
heatmap_plot_2
dev.off()

# Significance analysis ====
# Hypothesis test: Perform the Kruskal-Wallis rank sum test on your data
kruskal_result <- kruskal.test(ari_score ~ Algorithm, data = all_results_2_43)

# Print the results
print(kruskal_result)
kruskal_result$p.value < 0.05

# Post-Hoc Testing: Dunn's test
library(FSA)
dunn_test <- dunnTest(ari_score ~ Algorithm,
                      data = all_results_2_43,
                      method = "bonferroni")
save(dunn_test, file = "significance_analysis_dunn_test_all_cluster_ari.rda")
# CLD: Compact Letter Display
library(rcompanion)
cld_result <- cldList(P.adj ~ Comparison,
                      data = dunn_test$res,
                      threshold = 0.05)

save(cld_result, file = "significance_analysis_cld_letters_all_cluster_ari.rda")
