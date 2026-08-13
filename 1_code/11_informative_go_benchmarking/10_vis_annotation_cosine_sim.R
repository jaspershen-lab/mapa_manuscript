## Benchmarking on the informative GO set benchmark (k = 20)
## Step 10: Visualize the cosine similarity between each tool's module
## annotation and the ground-truth anchor annotation.
## Style follows 1_code/2_control_data/05_benchmarking/05_vis_annotation_cosine_sim.R.
## The control data had multiple human experts per module (-> per-module
## facets); here each module has a single ground-truth annotation, so the
## distribution across the 57 modules is shown in one panel per tool.

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(ggsignif)

load("3_data_analysis/11_informative_go_benchmarking/annotation_similarity_results.rda")

annotation_similarity_ordered <-
  annotation_similarity_results |>
  dplyr::mutate(method = factor(
    method,
    levels = c("mapa_cluster_label", "apear_cluster_label",
               "paver_cluster_label", "enrich_plot_cluster_label")
  ))

combined_p <-
  ggplot(annotation_similarity_ordered,
         aes(x = method, y = cosine_sim, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(fill = method),
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
  labs(
    x = NULL,
    y = "Cosine Similarity",
    fill = "Methods"
  ) +
  scale_fill_manual(
    values = methods,
    labels = c("MAPA", "aPEAR", "PAVER", "enrichplot")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2)) +
  guides(fill = guide_legend(nrow = 1))

combined_p

ggsave(plot = combined_p,
       filename = "3_data_analysis/11_informative_go_benchmarking/annotation_cosine_sim_vs_ground_truth.pdf",
       height = 5, width = 6)

## Wilcoxon tests (exact statistics)
mapa_sim <- annotation_similarity_results |>
  dplyr::filter(method == "mapa_cluster_label") |> dplyr::pull(cosine_sim)
for (m in c("apear_cluster_label", "paver_cluster_label", "enrich_plot_cluster_label")) {
  other_sim <- annotation_similarity_results |>
    dplyr::filter(method == m) |> dplyr::pull(cosine_sim)
  cat("MAPA vs", m, ": p =",
      format(wilcox.test(mapa_sim, other_sim)$p.value, digits = 3), "\n")
}
