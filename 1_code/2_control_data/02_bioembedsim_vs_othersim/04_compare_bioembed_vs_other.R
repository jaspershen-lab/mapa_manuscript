library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(reshape2)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/embedding_sim_df.rda"
)
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/jaccard_sim_df.rda"
)
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/kappa_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/op_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/dice_sim_df.rda")

setwd("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim")

# Compare embedding sim with overlap-based sim ====
## 1. Combine all the sim
df_lists <- list(embedding_sim_df,
                 jaccard_sim_df,
                 op_sim_df,
                 kappa_sim_df,
                 dice_sim_df)
suffixes <- c("_embedding", "_jaccard", "_op", "_kappa", "_dice")

renamed_dfs <- map2(df_lists, suffixes, function(df, suffix) {
  df %>% rename_with( ~ paste0("sim", suffix), .cols = "sim")
})

# Merge all dataframes
combined_sim_df <- purrr::reduce(renamed_dfs, function(x, y) {
  full_join(x, y, by = c("from", "to"))
})

# save(combined_sim_df, file = "combined_sim_df.rda")
load("combined_sim_df.rda")
export(combined_sim_df, file = "all_pathways_sim_res.xlsx")
## 2. Calculate Pearson correlations
cor_jaccard <- cor(
  combined_sim_df$sim_embedding,
  combined_sim_df$sim_jaccard,
  use = "complete.obs",
  method = "pearson"
)
cor_op <- cor(
  combined_sim_df$sim_embedding,
  combined_sim_df$sim_op,
  use = "complete.obs",
  method = "pearson"
)
cor_kappa <- cor(
  combined_sim_df$sim_embedding,
  combined_sim_df$sim_kappa,
  use = "complete.obs",
  method = "pearson"
)
cor_dice <- cor(
  combined_sim_df$sim_embedding,
  combined_sim_df$sim_dice,
  use = "complete.obs",
  method = "pearson"
)

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
  sim_value = c(
    combined_sim_df$sim_jaccard,
    combined_sim_df$sim_op,
    combined_sim_df$sim_kappa,
    combined_sim_df$sim_dice
  ),
  metric = factor(rep(
    c("Jaccard", "Overlap", "Kappa", "Dice"),
    each = nrow(combined_sim_df)
  ))
)

metric_colors <-
  c(
    "Jaccard" = "#e9c46b",
    "Overlap" = "#e66f51",
    "Kappa" = "#264653",
    "Dice" = "#299e8c"
  )

# Create the plot with all regression lines in one panel
p <- ggplot(plot_data, aes(x = sim_embedding, y = sim_value, color = metric)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              se = FALSE) +
  scale_color_manual(values = metric_colors) +
  labs(
    title = "Correlations between Pathway Biotext Embedding Similarity and Other Metrics",
    x = "Pathway Biotext Embedding Similarity",
    y = "Similarity Value",
    color = "Similarity Metric"
  ) +
  theme_bw() +
  # Add annotation with correlation values
  annotate(
    "text",
    x = rep(min(combined_sim_df$sim_embedding, na.rm = TRUE) + 0.05, 4),
    y = seq(
      from = max(
        c(
          combined_sim_df$sim_jaccard,
          combined_sim_df$sim_op,
          combined_sim_df$sim_kappa,
          combined_sim_df$sim_dice
        ),
        na.rm = TRUE
      ) - 0.05,
      to = max(
        c(
          combined_sim_df$sim_jaccard,
          combined_sim_df$sim_op,
          combined_sim_df$sim_kappa,
          combined_sim_df$sim_dice
        ),
        na.rm = TRUE
      ) - 0.35,
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

ggsave(
  plot = p,
  filename = "f01_bioembedsim_vs_othersim.pdf",
  width = 8,
  height = 8
)

control_dt

temp_data <-
  combined_sim_df %>%
  dplyr::left_join(control_dt[, c("id", "database", "name", "expected_module")], by = c("from" = "id")) %>%
  dplyr::left_join(control_dt[, c("id", "database", "name", "expected_module")], by = c("to" = "id")) %>%
  dplyr::mutate(
    class = case_when(
      expected_module.x == expected_module.y ~ "same_module",
      expected_module.x != expected_module.y ~ "different_module"
    )
  ) %>%
  dplyr::mutate(
    class2 = case_when(
      database.x == database.y ~ "intra_database",
      database.x != database.y ~ "inter_database"
    )
  ) %>%
  dplyr::mutate(class = factor(class, levels = c("same_module", "different_module")))

temp_data$class %>%
  table()

library(gghalves)

library(ggsignif)


####calculate the correlation between each two similarity
library(GGally)

temp_data

plot <-
  ggpairs(data = temp_data,
          columns = 3:7,
          ggplot2::aes(color = class)) +
  theme_bw() +
  scale_color_manual(values = same_different_module_color) +
  scale_fill_manual(values = same_different_module_color)

plot

ggsave(
  plot = plot,
  filename = "pathbiotext_vs_othersim_scatterplot.pdf",
  width = 10,
  height = 10
)


####add median values for each group
plot <-
  temp_data %>%
  ggplot(aes(x = class, y = sim_embedding)) +
  geom_half_point(
    aes(fill = class),
    side = "l",
    position = position_nudge(x = -0.1),
    color = "black",
    size = 4,
    shape = 21,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_half_boxplot(
    outlier.shape = NA,
    side = "l",
    color = "black",
    fill = NA,
    show.legend = FALSE
  ) +
  geom_half_violin(
    aes(fill = class),
    side = "r",
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_fill_manual(values = same_different_module_color) +
  scale_color_manual(values = same_different_module_color) +
  labs(x = "", y = "PathBiotext") +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.5,
    size = 3,
    color = "black"
  ) +
  ggsignif::geom_signif(
    test = "wilcox.test",
    comparisons = list(c("same_module", "different_module")),
    map_signif_level = TRUE,
    textsize = 3.5,
    y_position = 0.95,
    vjust = 0.5
  )
plot

ggsave(
  plot = plot,
  filename = "pathbiotext_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)




####add median values for each group
plot <-
  temp_data %>%
  ggplot(aes(x = class, y = sim_jaccard)) +
  geom_half_point(
    aes(fill = class),
    side = "l",
    position = position_nudge(x = -0.1),
    color = "black",
    size = 4,
    shape = 21,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_half_boxplot(
    outlier.shape = NA,
    side = "l",
    color = "black",
    fill = NA,
    show.legend = FALSE
  ) +
  geom_half_violin(
    aes(fill = class),
    side = "r",
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_fill_manual(values = same_different_module_color) +
  scale_color_manual(values = same_different_module_color) +
  labs(x = "", y = "Jaccard index") +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.5,
    size = 3,
    color = "black"
  ) +
  ggsignif::geom_signif(
    test = "wilcox.test",
    comparisons = list(c("same_module", "different_module")),
    map_signif_level = TRUE,
    textsize = 3.5,
    y_position = 0.95,
    vjust = 0.5
  )
plot

ggsave(
  plot = plot,
  filename = "jaccard_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)


####add median values for each group
plot <-
  temp_data %>%
  ggplot(aes(x = class, y = sim_op)) +
  geom_half_point(
    aes(fill = class),
    side = "l",
    position = position_nudge(x = -0.1),
    color = "black",
    size = 4,
    shape = 21,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_half_boxplot(
    outlier.shape = NA,
    side = "l",
    color = "black",
    fill = NA,
    show.legend = FALSE
  ) +
  geom_half_violin(
    aes(fill = class),
    side = "r",
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_fill_manual(values = same_different_module_color) +
  scale_color_manual(values = same_different_module_color) +
  labs(x = "", y = "Overlap cofficient") +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.5,
    size = 3,
    color = "black"
  ) +
  ggsignif::geom_signif(
    test = "wilcox.test",
    comparisons = list(c("same_module", "different_module")),
    map_signif_level = TRUE,
    textsize = 3.5,
    y_position = 0.95,
    vjust = 0.5
  )
plot

ggsave(
  plot = plot,
  filename = "overlap_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)


####sim_kappa
plot <-
  temp_data %>%
  ggplot(aes(x = class, y = sim_kappa)) +
  geom_half_point(
    aes(fill = class),
    side = "l",
    position = position_nudge(x = -0.1),
    color = "black",
    size = 4,
    shape = 21,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_half_boxplot(
    outlier.shape = NA,
    side = "l",
    color = "black",
    fill = NA,
    show.legend = FALSE
  ) +
  geom_half_violin(
    aes(fill = class),
    side = "r",
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_fill_manual(values = same_different_module_color) +
  scale_color_manual(values = same_different_module_color) +
  labs(x = "", y = "Kappa") +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.5,
    size = 3,
    color = "black"
  ) +
  ggsignif::geom_signif(
    test = "wilcox.test",
    comparisons = list(c("same_module", "different_module")),
    map_signif_level = TRUE,
    textsize = 3.5,
    y_position = 0.95,
    vjust = 0.5
  )
plot

ggsave(
  plot = plot,
  filename = "kappa_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)



###sim_kappa
plot <-
  temp_data %>%
  ggplot(aes(x = class, y = sim_dice)) +
  geom_half_point(
    aes(fill = class),
    side = "l",
    position = position_nudge(x = -0.1),
    color = "black",
    size = 4,
    shape = 21,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_half_boxplot(
    outlier.shape = NA,
    side = "l",
    color = "black",
    fill = NA,
    show.legend = FALSE
  ) +
  geom_half_violin(
    aes(fill = class),
    side = "r",
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw() +
  scale_fill_manual(values = same_different_module_color) +
  scale_color_manual(values = same_different_module_color) +
  labs(x = "", y = "Dice") +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.5,
    size = 3,
    color = "black"
  ) +
  ggsignif::geom_signif(
    test = "wilcox.test",
    comparisons = list(c("same_module", "different_module")),
    map_signif_level = TRUE,
    textsize = 3.5,
    y_position = 0.95,
    vjust = 0.5
  )
plot

ggsave(
  plot = plot,
  filename = "dice_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)



