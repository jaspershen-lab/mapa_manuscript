library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

setwd("3_data_analysis/02_control_data/05_benchmarking/comparison_result/")

load("combined_similarity_results.rda")

# combined_similarity_results_ordered <- combined_similarity_results |>
#   mutate(tool_name = factor(tool_name,
#                             levels = c("mapa_cluster_label", "apear_cluster_label", "paver_cluster_label", "enrich_plot_cluster_label")))
#
# combined_p <- ggplot(combined_similarity_results_ordered, aes(x = factor(tool_cluster), y = cosine_sim, fill = tool_name)) +
#   geom_boxplot(alpha = 0.7, position = "dodge") +
#   geom_point(aes(fill = tool_name),
#              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2),
#              size = 1,
#              alpha = 0.6,
#              shape = 21,
#              color = "black") +
#   labs(
#     x = "Functional Module",
#     y = "Cosine Similarity"
#   ) +
#   scale_fill_manual(values = methods) +
#   theme_bw() +
#   theme(
#     axis.text = element_text(size = 8),
#     axis.title = element_text(size = 8, face = "bold"),
#     legend.position = "bottom",
#     legend.text = element_text(size = 9),
#     legend.title = element_blank()
#   ) +
#   scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))
#

combined_similarity_results_ordered <- combined_similarity_results |>
  mutate(tool_name = factor(tool_name,
                            levels = c("mapa_cluster_label", "apear_cluster_label", "paver_cluster_label", "enrich_plot_cluster_label")))

combined_similarity_results_ordered_remove_jc_5 <-
  combined_similarity_results_ordered |>
  filter(!(expert_name == "JiangChao" & expert_cluster == 5))

library(ggsignif)

cluster_labeller <- function(variable, value) {
  return(paste("Functional Module", value))
}

combined_p <- ggplot(combined_similarity_results_ordered_remove_jc_5, aes(x = tool_name, y = cosine_sim, fill = tool_name)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(aes(fill = tool_name),
             position = position_jitter(width = 0.2),
             size = 1,
             alpha = 0.6,
             shape = 21,
             color = "black") +
  geom_signif(
    comparisons = list(
      c("mapa_cluster_label", "apear_cluster_label"),
      c("mapa_cluster_label", "paver_cluster_label"),
      c("mapa_cluster_label", "enrich_plot_cluster_label")
    ),
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    step_increase = 0.1,
    tip_length = 0.01,
    size = 0.5,
    textsize = 2.5
  ) +
  facet_wrap(~ tool_cluster, labeller = cluster_labeller, scales = "free_x", nrow = 2) +
  labs(
    x = NULL,  # Remove x-axis title
    y = "Cosine Similarity",
    fill = "Methods"  # Legend title
  ) +
  scale_fill_manual(
    values = methods,
    labels = c("MAPA", "aPEAR", "PAVER", "enrichplot")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    # axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    legend.position = "bottom",  # Show legend at bottom
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    strip.text = element_text(size = 8)
  ) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  guides(fill = guide_legend(nrow = 1))

combined_p
ggsave(plot = combined_p, filename = "remove_na_fm_5_jc_combined_p.pdf", height = 6, width = 8)

fm5 <- combined_similarity_results_ordered_remove_jc_5 |>
  filter(tool_cluster == 5)

