library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(reshape2)

load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/embedding_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/jaccard_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/kappa_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/op_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/dice_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/wang_sim_df.rda")

df_lists <- list(wang_sim_df, embedding_sim_df, jaccard_sim_df, op_sim_df, kappa_sim_df, dice_sim_df)
suffixes <- c("_wang", "_embedding", "_jaccard", "_op", "_kappa", "_dice")

# Rename the sim columns
renamed_dfs <- map2(dataframes, suffixes, function(df, suffix) {
  df %>% rename_with(~paste0("sim", suffix), .cols = "sim")
})
# Reduce to merge all dataframes
go_combined_sim_df <- purrr::reduce(renamed_dfs, function(x, y) {
  left_join(x, y, by = c("from", "to"))
})

save(go_combined_sim_df, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/go_combined_sim_df.rda")

# Calculate Pearson correlations
cor_semantic <- cor(go_combined_sim_df$sim_embedding, go_combined_sim_df$sim_wang, use = "complete.obs", method = "pearson")
cor_op <- cor(go_combined_sim_df$sim_embedding, go_combined_sim_df$sim_op, use = "complete.obs", method = "pearson")
cor_kappa <- cor(go_combined_sim_df$sim_embedding, go_combined_sim_df$sim_kappa, use = "complete.obs", method = "pearson")
cor_dice <- cor(go_combined_sim_df$sim_embedding, go_combined_sim_df$sim_dice, use = "complete.obs", method = "pearson")
cor_jaccard <- cor(go_combined_sim_df$sim_embedding, go_combined_sim_df$sim_jaccard, use = "complete.obs", method = "pearson")

# Create a data frame for correlation results
cor_results <- data.frame(
  metric = c("Wang", "Overlap", "Kappa", "Dice", "Jaccard"),
  correlation = c(cor_semantic, cor_op, cor_kappa, cor_dice, cor_jaccard)
)

# Print the correlation results
print(cor_results)

metric_colors <- c("Wang" = "#f2a361", "Jaccard" = "#e9c46b", "Overlap" = "#e66f51", "Kappa" = "#264653", "Dice" = "#299e8c")

plot_data <- data.frame(
  sim_embedding = rep(go_combined_sim_df$sim_embedding, 5),
  sim_value = c(go_combined_sim_df$sim_wang, go_combined_sim_df$sim_jaccard, go_combined_sim_df$sim_op, go_combined_sim_df$sim_kappa, go_combined_sim_df$sim_dice),
  metric = factor(rep(c("Wang", "Jaccard", "Overlap", "Kappa", "Dice"), each = nrow(go_combined_sim_df)))
)

# Create the plot with all regression lines in one panel
p <- ggplot(plot_data, aes(x = sim_embedding, y = sim_value, color = metric)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_color_manual(values = metric_colors) +
  labs(
    title = "Correlation between Pathway Biotext Embedding Similarity and Other Metrics (GO terms)",
    x = "Pathway Biotext Embedding Similarity",
    y = "Similarity Value",
    color = "Similarity Metric"
  ) +
  theme_minimal() +
  # Add annotation with correlation values
  annotate(
    "text",
    x = rep(min(go_combined_sim_df$sim_embedding, na.rm = TRUE) + 0.05, 5),
    y = seq(
      from = max(c(go_combined_sim_df$sim_wang, go_combined_sim_df$sim_jaccard, go_combined_sim_df$sim_op, go_combined_sim_df$sim_kappa, go_combined_sim_df$sim_dice), na.rm = TRUE) - 0.05,
      to = max(c(go_combined_sim_df$sim_wang, go_combined_sim_df$sim_jaccard, go_combined_sim_df$sim_op, go_combined_sim_df$sim_kappa, go_combined_sim_df$sim_dice), na.rm = TRUE) - 0.35,
      length.out = 5
    ),
    label = sprintf(
      "%s: r = %.3f",
      c("Wang", "Jaccard", "Overlap", "Kappa", "Dice"),
      c(cor_semantic, cor_jaccard, cor_op, cor_kappa, cor_dice)
    ),
    color = metric_colors,
    hjust = 0,
    fontface = "bold",
    size = 4
  )

p

ggsave(plot = p, filename = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/f02_bioembedsim_vs_wang_go_terms.pdf", width = 8, height = 8)
