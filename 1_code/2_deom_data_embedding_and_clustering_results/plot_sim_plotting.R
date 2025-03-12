# Script to compare two pairwise similarity matrices

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(corrplot)

gemini_sim_matrix <- semantic_sim_matrix
openai_sim_matrix <- semantic_sim_matrix
openai_sim_matrix_with_genename <- sim_matrix

sim_matrix1 <- gemini_sim_matrix
sim_matrix2 <- openai_sim_matrix_with_genename

# Function to extract pairwise values from matrices
extract_pairwise_values <- function(matrix1, matrix2) {
  # Ensure both matrices have the same dimensions
  if (!all(dim(matrix1) == dim(matrix2)) ||
      !all(rownames(matrix1) == rownames(matrix2)) ||
      !all(colnames(matrix1) == colnames(matrix2))) {
    stop("Matrices must have the same dimensions and row/column names")
  }

  # Get the lower triangle indices (exclude diagonal)
  indices <- which(lower.tri(matrix1), arr.ind = TRUE)

  # Extract values
  values1 <- sapply(1:nrow(indices), function(i) {
    matrix1[indices[i, 1], indices[i, 2]]
  })

  values2 <- sapply(1:nrow(indices), function(i) {
    matrix2[indices[i, 1], indices[i, 2]]
  })

  # Create pair labels
  pairs <- sapply(1:nrow(indices), function(i) {
    paste(rownames(matrix1)[indices[i, 1]], colnames(matrix1)[indices[i, 2]], sep = "-")
  })

  return(data.frame(pair = pairs, sim1 = values1, sim2 = values2))
}

# Extract pairwise values
pair_data <- extract_pairwise_values(sim_matrix1, sim_matrix2)

# Calculate correlations
pearson_cor <- cor(pair_data$sim1, pair_data$sim2, method = "pearson")
spearman_cor <- cor(pair_data$sim1, pair_data$sim2, method = "spearman")
kendall_cor <- cor(pair_data$sim1, pair_data$sim2, method = "kendall")

# Print correlation results
cat("Correlation between similarity matrices:\n")
cat("Pearson correlation:", round(pearson_cor, 4), "\n")
cat("Spearman correlation:", round(spearman_cor, 4), "\n")
cat("Kendall correlation:", round(kendall_cor, 4), "\n")

# Visualize the correlation with a scatter plot
p1 <- ggplot(pair_data, aes(x = sim1, y = sim2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(
    title = "Correlation between Similarity Matrices",
    subtitle = paste("Pearson correlation:", round(pearson_cor, 4)),
    x = "Similarity values (Matrix 1)",
    y = "Similarity values (Matrix 2)"
  ) +
  theme_minimal()

# Create a heatmap of differences
diff_matrix <- sim_matrix1 - sim_matrix2
melted_diff <- melt(diff_matrix)
colnames(melted_diff) <- c("Row", "Column", "Difference")

p2 <- ggplot(melted_diff, aes(x = Column, y = Row, fill = Difference)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Difference between Similarity Matrices") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot both matrices side by side for visual comparison
par(mfrow = c(1, 2))
corrplot(sim_matrix1, method = "color", type = "lower",
         title = "Matrix 1", mar = c(0,0,1,0))
corrplot(sim_matrix2, method = "color", type = "lower",
         title = "Matrix 2", mar = c(0,0,1,0))

# Reset par
par(mfrow = c(1, 1))

# Print the scatter plot
print(p1)

# Print the difference heatmap
print(p2)

# Optional: Save the visualizations
# ggsave("similarity_correlation_scatter.png", p1, width = 8, height = 6)
# ggsave("similarity_difference_heatmap.png", p2, width = 8, height = 6)


sim_matrix1 <- openai_sim_matrix
sim_matrix2 <- openai_sim_matrix_with_genename
sim_matrix3 <- gemini_sim_matrix

# Function to extract similarity values from a matrix (lower triangle)
extract_similarity_values <- function(matrix) {
  # Get the lower triangle indices (exclude diagonal)
  indices <- which(lower.tri(matrix), arr.ind = TRUE)

  # Extract values
  values <- sapply(1:nrow(indices), function(i) {
    matrix[indices[i, 1], indices[i, 2]]
  })

  return(values)
}

# Extract similarity values from each matrix
values1 <- extract_similarity_values(sim_matrix1)
values2 <- extract_similarity_values(sim_matrix2)
values3 <- extract_similarity_values(sim_matrix3)

# Combine the values into a data frame for plotting
similarity_data <- data.frame(
  value = c(values1, values2, values3),
  matrix = factor(c(
    rep("openai", length(values1)),
    rep("openai_with_gene_name", length(values2)),
    rep("gemini", length(values3))
  ))
)

# Calculate summary statistics
summary_stats <- similarity_data %>%
  group_by(matrix) %>%
  summarize(
    min = min(value),
    q1 = quantile(value, 0.25),
    median = median(value),
    mean = mean(value),
    q3 = quantile(value, 0.75),
    max = max(value),
    sd = sd(value)
  )

# Print summary statistics
print(summary_stats)

# Create a single boxplot with custom colors
custom_colors <- c("#1E88E5", "#FFC107", "#D81B60")  # Blue, Amber, Pink

p <- ggplot(similarity_data, aes(x = matrix, y = value, fill = matrix)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5) +
  labs(
    title = "Comparison of Similarity Values Across Three Matrices",
    x = "Similairty Metrics",
    y = "Similarity Value"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Display the plot
print(p)

# Optional: Add notches to show confidence intervals around median
p_notched <- ggplot(similarity_data, aes(x = matrix, y = value, fill = matrix)) +
  geom_boxplot(width = 0.7, outlier.size = 1.5, notch = TRUE) +
  labs(
    title = "Comparison of Similarity Values Across Three Matrices",
    subtitle = "Notches show 95% confidence interval around median",
    x = "Matrix",
    y = "Similarity Value"
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Display the notched plot
print(p_notched)

# Optional: Save the plots
# ggsave("similarity_boxplot_colored.png", p, width = 8, height = 6, dpi = 300)
# ggsave("similarity_boxplot_notched.png", p_notched, width = 8, height = 6, dpi = 300)
