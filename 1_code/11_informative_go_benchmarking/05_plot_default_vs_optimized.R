## Benchmarking on the informative GO set benchmark (k = 20)
## Step 5: Plot default vs optimized ARI for the four methods.
## For each method two points are shown (open = default, filled = optimized)
## connected by an arrow pointing from the default to the optimized result.

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/11_informative_go_benchmarking/default_vs_optimized_ari.rda")

names(methods) <- c("MAPA", "aPEAR", "PAVER", "enrichplot")

plot_data <-
  default_vs_optimized_ari |>
  ## Discrete y-axis levels are drawn from bottom to top, so reverse the
  ## requested top-to-bottom order here.
  dplyr::mutate(methods = factor(
    methods,
    levels = rev(c("MAPA", "PAVER", "aPEAR", "enrichplot"))
  ))

point_data <-
  plot_data |>
  tidyr::pivot_longer(cols = c(default_ari, optimized_ari),
                      names_to = "type",
                      values_to = "ARI") |>
  dplyr::mutate(type = dplyr::recode(type,
                                     default_ari = "Default",
                                     optimized_ari = "Optimized"))

p <-
  ggplot(plot_data, aes(y = methods)) +
  geom_segment(
    ## Trim both ends so the arrow does not pierce the two points
    aes(yend = methods,
        x = default_ari + 0.02 * sign(optimized_ari - default_ari),
        xend = optimized_ari - 0.018 * sign(optimized_ari - default_ari),
        color = methods),
    arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_point(
    data = point_data,
    aes(x = ARI, color = methods, fill = methods, shape = type),
    size = 3,
    stroke = 0.5
  ) +
  scale_shape_manual(values = c(Default = 21, Optimized = 16)) +
  scale_color_manual(values = methods, guide = "none") +
  scale_fill_manual(values = scales::alpha(methods, 0), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  geom_text(
    aes(x = default_ari, label = sprintf("%.2f", default_ari)),
    hjust = 1.4, size = 2.5
  ) +
  geom_text(
    aes(x = optimized_ari, label = sprintf("%.2f", optimized_ari)),
    hjust = -0.4, size = 2.5
  ) +
  guides(shape = guide_legend(
    title = NULL,
    override.aes = list(color = "black", fill = NA)
  )) +
  labs(x = "Adjusted Rand Index (ARI)",
       y = "Methods") +
  theme_bw() +
  theme(legend.position = "top")
p

ggsave(plot = p,
       filename = "3_data_analysis/11_informative_go_benchmarking/compare_best_ari.pdf",
       width = 5, height = 6)
