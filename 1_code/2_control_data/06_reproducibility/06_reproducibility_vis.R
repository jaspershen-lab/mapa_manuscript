library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

setwd("3_data_analysis/08_reproducibility")

load("pairwise_annotation_comparisons_sim_res.rda")

all_node_sets <- unique(pairwise_annotation_comparisons_sim_res$node_set_1)

library(gghalves)
library(ggsignif)

# Create a mapping for node set names to more readable labels
node_set_labels <- paste("Node Set", 1:8)
names(node_set_labels) <- all_node_sets

# Prepare the data with facet labels
all_node_set_sim <- data.frame()

for (node_set in all_node_sets) {
  sim_1 <- pairwise_annotation_comparisons_sim_res |>
    filter(node_set_1 == node_set) |>
    mutate(node_set_label = unname(node_set_labels[node_set])) |>
    select(node_set_label, sim, class)
  sim_2 <- pairwise_annotation_comparisons_sim_res |>
    filter(node_set_2 == node_set) |>
    mutate(node_set_label = unname(node_set_labels[node_set])) |>
    select(node_set_label, sim, class)
  node_set_all_info <- rbind(sim_1, sim_2)

  all_node_set_sim <- rbind(all_node_set_sim,
                            node_set_all_info)
}

plot_data <- all_node_set_sim
save(plot_data, file = "plot_data.rda")

same_different_module_color <-
  c(
    same_fm = "#ff0000",
    diff_fm = "#65684c"
  )

actual_classes <- unique(all_node_set_sim$class)

plot_data <- plot_data |>
  mutate(class = factor(class, levels = c("same_fm", "diff_fm")))

plot_data_same <- plot_data |>
  filter(class == "same_fm") |>
  distinct(node_set_label, sim, .keep_all = TRUE)

plot_data_remove_duplicate <- plot_data |>
  filter(class == "diff_fm") |>
  rbind(plot_data_same)

# Create the faceted plot
faceted_plot <- plot_data_remove_duplicate |>
  ggplot(aes(x = class, y = sim)) +
  geom_half_point(
    aes(fill = class),
    side = "l",
    position = position_nudge(x = -0.1),
    color = "black",
    size = 3,
    shape = 21,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_half_boxplot(
    outlier.shape = NA,
    side = "l",
    color = "black",
    fill = NA,
    show.legend = FALSE,
    width = 0.4
  ) +
  geom_half_violin(
    aes(fill = class),
    side = "r",
    alpha = 0.5,
    show.legend = FALSE,
    width = 0.8
  ) +
  scale_y_continuous(limits = c(0.2, 1.1)) +  # Increased upper limit for significance
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = same_different_module_color) +
  labs(
    x = "",
    y = "Cosine Similarity"
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.5,
    size = 2.5,
    color = "black"
  ) +
  facet_wrap(~ node_set_label, ncol = 4, scales = "free") +
  geom_signif(
    comparisons = list(actual_classes),  # Use actual class names from data
    test = "wilcox.test",
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    step_increase = 0.1,
    tip_length = 0.01,
    size = 0.5,
    textsize = 2.5,
    y_position = 1.05  # Position significance at top
  )

faceted_plot

ggsave(filename = "annotation_producibility_plot.pdf",
       plot = faceted_plot,
       width = 12, height = 7)

