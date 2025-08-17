library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(reshape2)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/embedding_sim_df.rda"
)

###GO similarity
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/wang_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/resnik_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/rel_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/jiang_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/lin_sim_df.rda")
###component based
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/op_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/kappa_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/dice_sim_df.rda")
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/jaccard_sim_df.rda"
)

dir.create("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/GO_terms")
setwd("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/GO_terms")

embedding_sim_df <-
  embedding_sim_df %>%
  dplyr::filter(stringr::str_detect(from, "GO:") &
                  stringr::str_detect(to, "GO:"))

df_lists <- list(
  embedding_sim_df,
  wang_sim_df,
  resnik_sim_df,
  rel_sim_df,
  jiang_sim_df,
  lin_sim_df,
  jaccard_sim_df,
  op_sim_df,
  kappa_sim_df,
  dice_sim_df
)

suffixes <- c(
  "_embedding",
  "_wang",
  "_resnik",
  "_rel",
  "_jiang",
  "_lin",
  "_jaccard",
  "_op",
  "_kappa",
  "_dice"
)

# Rename the sim columns
renamed_dfs <- map2(df_lists, suffixes, function(df, suffix) {
  df %>% rename_with(~ paste0("sim", suffix), .cols = "sim")
})

# Reduce to merge all dataframes
go_combined_sim_df <- purrr::reduce(renamed_dfs, function(x, y) {
  left_join(x, y, by = c("from", "to"))
})

save(go_combined_sim_df, file = "go_combined_sim_df.rda")

load("go_combined_sim_df.rda")
rio::export(go_combined_sim_df, file = "go_term_sim_res.xlsx")

###remove NA values
go_combined_sim_df <- go_combined_sim_df %>%
  dplyr::filter(
    !is.na(sim_embedding) &
      !is.na(sim_wang) &
      !is.na(sim_resnik) &
      !is.na(sim_rel) &
      !is.na(sim_jiang) &
      !is.na(sim_lin) &
      !is.na(sim_jaccard) &
      !is.na(sim_op) &
      !is.na(sim_kappa) &
      !is.na(sim_dice)
  )

# # Calculate Pearson correlations
# cor_semantic <- cor(
#   go_combined_sim_df$sim_embedding,
#   go_combined_sim_df$sim_wang,
#   use = "complete.obs",
#   method = "pearson"
# )
#
# cor_op <- cor(
#   go_combined_sim_df$sim_embedding,
#   go_combined_sim_df$sim_op,
#   use = "complete.obs",
#   method = "pearson"
# )
# cor_kappa <- cor(
#   go_combined_sim_df$sim_embedding,
#   go_combined_sim_df$sim_kappa,
#   use = "complete.obs",
#   method = "pearson"
# )
# cor_dice <- cor(
#   go_combined_sim_df$sim_embedding,
#   go_combined_sim_df$sim_dice,
#   use = "complete.obs",
#   method = "pearson"
# )
# cor_jaccard <- cor(
#   go_combined_sim_df$sim_embedding,
#   go_combined_sim_df$sim_jaccard,
#   use = "complete.obs",
#   method = "pearson"
# )
#
# # Create a data frame for correlation results
# cor_results <- data.frame(
#   metric = c("Wang", "Overlap", "Kappa", "Dice", "Jaccard"),
#   correlation = c(cor_semantic, cor_op, cor_kappa, cor_dice, cor_jaccard)
# )
#
# # Print the correlation results
# print(cor_results)
#
# metric_colors <- c(
#   "Wang" = "#f2a361",
#   "Jaccard" = "#e9c46b",
#   "Overlap" = "#e66f51",
#   "Kappa" = "#264653",
#   "Dice" = "#299e8c"
# )
#
# plot_data <- data.frame(
#   sim_embedding = rep(go_combined_sim_df$sim_embedding, 5),
#   sim_value = c(
#     go_combined_sim_df$sim_wang,
#     go_combined_sim_df$sim_jaccard,
#     go_combined_sim_df$sim_op,
#     go_combined_sim_df$sim_kappa,
#     go_combined_sim_df$sim_dice
#   ),
#   metric = factor(rep(
#     c("Wang", "Jaccard", "Overlap", "Kappa", "Dice"),
#     each = nrow(go_combined_sim_df)
#   ))
# )
#
# # Create the plot with all regression lines in one panel
# p <- ggplot(plot_data, aes(x = sim_embedding, y = sim_value, color = metric)) +
#   geom_point(alpha = 0.3) +
#   geom_smooth(method = "lm",
#               formula = y ~ x,
#               se = FALSE) +
#   scale_color_manual(values = metric_colors) +
#   labs(
#     title = "Correlation between Pathway Biotext Embedding Similarity and Other Metrics (GO terms)",
#     x = "Pathway Biotext Embedding Similarity",
#     y = "Similarity Value",
#     color = "Similarity Metric"
#   ) +
#   theme_minimal() +
#   # Add annotation with correlation values
#   annotate(
#     "text",
#     x = rep(min(
#       go_combined_sim_df$sim_embedding, na.rm = TRUE
#     ) + 0.05, 5),
#     y = seq(
#       from = max(
#         c(
#           go_combined_sim_df$sim_wang,
#           go_combined_sim_df$sim_jaccard,
#           go_combined_sim_df$sim_op,
#           go_combined_sim_df$sim_kappa,
#           go_combined_sim_df$sim_dice
#         ),
#         na.rm = TRUE
#       ) - 0.05,
#       to = max(
#         c(
#           go_combined_sim_df$sim_wang,
#           go_combined_sim_df$sim_jaccard,
#           go_combined_sim_df$sim_op,
#           go_combined_sim_df$sim_kappa,
#           go_combined_sim_df$sim_dice
#         ),
#         na.rm = TRUE
#       ) - 0.35,
#       length.out = 5
#     ),
#     label = sprintf(
#       "%s: r = %.3f",
#       c("Wang", "Jaccard", "Overlap", "Kappa", "Dice"),
#       c(cor_semantic, cor_jaccard, cor_op, cor_kappa, cor_dice)
#     ),
#     color = metric_colors,
#     hjust = 0,
#     fontface = "bold",
#     size = 4
#   )
#
# p
#
# ggsave(
#   plot = p,
#   filename = "f02_bioembedsim_vs_wang_go_terms.pdf",
#   width = 8,
#   height = 8
# )



temp_data <-
  go_combined_sim_df %>%
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
          columns = 3:12,
          ggplot2::aes(color = class)) +
  theme_bw() +
  scale_color_manual(values = same_different_module_color) +
  scale_fill_manual(values = same_different_module_color)

plot

ggsave(
  plot = plot,
  filename = "pathbiotext_vs_othersim_scatterplot.pdf",
  width = 15,
  height = 15
)


####wang vs resnik
cor.test(temp_data$sim_wang, temp_data$sim_resnik, method = "pearson")
cor.test(temp_data$sim_wang, temp_data$sim_rel, method = "pearson")
cor.test(temp_data$sim_wang, temp_data$sim_jiang, method = "pearson")





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
  ggplot(aes(x = class, y = sim_wang)) +
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
  labs(x = "", y = "Wang") +
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
  filename = "wang_same_different_module_comparison.pdf",
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
  labs(x = "", y = "Overlap coefficient") +
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
  filename = "op_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)


temp_data %>%
  dplyr::mutate(variable_id = paste0(from, "_", to)) %>%
  dplyr::select(
    variable_id,
    expected_module.x,
    expected_module.y,
    sim_embedding,
    sim_wang,
    sim_jaccard,
    sim_op,
    class
  ) %>%
  tidyr::pivot_longer(
    cols = c(sim_embedding, sim_wang, sim_jaccard, sim_op),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(metric = factor(
    metric,
    levels = c("sim_embedding", "sim_wang", "sim_jaccard", "sim_op")
  )) %>%
  ggplot(aes(metric, value)) +
  geom_point(aes(color = class))

####UMAP
