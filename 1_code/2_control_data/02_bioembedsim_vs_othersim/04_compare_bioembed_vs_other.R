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

# Compare embedding sim with overlap-based sim ====
## 1. Combine all the sim
df_lists <- list(embedding_sim_df, jaccard_sim_df, op_sim_df, kappa_sim_df, dice_sim_df)
suffixes <- c("_embedding", "_jaccard", "_op", "_kappa", "_dice")

renamed_dfs <- map2(df_lists, suffixes, function(df, suffix) {
  df %>% rename_with(~paste0("sim", suffix), .cols = "sim")
})

# Merge all dataframes
combined_sim_df <- purrr::reduce(renamed_dfs, function(x, y) {
  full_join(x, y, by = c("from", "to"))
})

save(combined_sim_df, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/combined_sim_df.rda")

## 2. Calculate Pearson correlations
cor_jaccard <- cor(combined_sim_df$sim_embedding, combined_sim_df$sim_jaccard, use = "complete.obs", method = "pearson")
cor_op <- cor(combined_sim_df$sim_embedding, combined_sim_df$sim_op, use = "complete.obs", method = "pearson")
cor_kappa <- cor(combined_sim_df$sim_embedding, combined_sim_df$sim_kappa, use = "complete.obs", method = "pearson")
cor_dice <- cor(combined_sim_df$sim_embedding, combined_sim_df$sim_dice, use = "complete.obs", method = "pearson")

# Create a data frame for correlation results
cor_results <- data.frame(
  metric = c("Jaccard", "Overlap", "Kappa", "Dice"),
  correlation = c(cor_jaccard, cor_op, cor_kappa, cor_dice)
)

# Print the correlation results
print(cor_results)

# Reshape data for plotting
plot_data <- data.frame(
  sim_embedding = rep(combined_sim_df$sim_embedding, 4),
  sim_value = c(combined_sim_df$sim_jaccard, combined_sim_df$sim_op, combined_sim_df$sim_kappa, combined_sim_df$sim_dice),
  metric = factor(rep(c("Jaccard", "Overlap", "Kappa", "Dice"), each = nrow(combined_sim_df)))
)

metric_colors <-
  c("Jaccard" = "#e9c46b",
    "Overlap" = "#e66f51",
    "Kappa" = "#264653",
    "Dice" = "#299e8c")

# Create the plot with all regression lines in one panel
p <- ggplot(plot_data,
            aes(x = sim_embedding, y = sim_value, color = metric)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_color_manual(values = metric_colors) +
  labs(
    title = "Correlations between Pathway Biotext Embedding Similarity and Other Metrics",
    x = "Pathway Biotext Embedding Similarity",
    y = "Similarity Value",
    color = "Similarity Metric"
  ) +
  theme_minimal() +
  # Add annotation with correlation values
  annotate(
    "text",
    x = rep(min(combined_sim_df$sim_embedding, na.rm = TRUE) + 0.05, 4),
    y = seq(
      from = max(c(combined_sim_df$sim_jaccard, combined_sim_df$sim_op, combined_sim_df$sim_kappa, combined_sim_df$sim_dice), na.rm = TRUE) - 0.05,
      to = max(c(combined_sim_df$sim_jaccard, combined_sim_df$sim_op, combined_sim_df$sim_kappa, combined_sim_df$sim_dice), na.rm = TRUE) - 0.35,
      length.out = 4
    ),
    label = sprintf(
      "%s: r = %.3f",
      c("Jaccard", "Overlap", "Kappa", "Dice"),
      c(cor_jaccard, cor_op, cor_kappa, cor_dice)
    ),
    color = metric_colors,
    hjust = 0,
    fontface = "bold",
    size = 4
  )

p

ggsave(plot = p, filename = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/f01_bioembedsim_vs_othersim.pdf", width = 8, height = 8)
