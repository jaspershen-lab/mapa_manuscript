## Benchmarking on the informative GO set benchmark (k = 20)
## Step 8: Wang semantic similarity (GO-graph-based, GOSemSim) between GO
## terms in the same vs different ground-truth modules.
## Calculation follows 1_code/2_control_data/02_bioembedsim_vs_othersim/
## 03_go_term_sim_calculation.R; the plot follows the style of
## 05_compare_bioembed_vs_go_terms.R (same adaptation as script 07:
## geom_half_point replaced by manually jittered geom_point)

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(GOSemSim)
library(org.Hs.eg.db)
library(gghalves)
library(ggsignif)

ground_truth_dt <- readr::read_csv("2_data/informative_go_set/k_20/ground_truth.csv",
                                   show_col_types = FALSE)
go_ids <- ground_truth_dt$go_id

# 1. Wang similarity ====
## Wang similarity only uses the GO DAG structure, no IC needed
bp_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "BP", computeIC = FALSE)

wang_sim_matrix <-
  termSim(
    t1 = go_ids,
    t2 = go_ids,
    semData = bp_semgodata,
    method = "Wang"
  )

## Terms missing from the installed GO.db release (if any) give NA rows
n_na_terms <- sum(rowSums(!is.na(wang_sim_matrix)) == 0)
message("Terms with no Wang similarity (not in installed GO.db): ", n_na_terms)

wang_sim_df <-
  wang_sim_matrix |>
  as.data.frame.table(responseName = "sim") |>
  dplyr::filter(Var1 != Var2) |>
  dplyr::rename(from = Var1, to = Var2) |>
  dplyr::mutate(across(c(from, to), as.character)) |>
  dplyr::filter(from < to)

save(wang_sim_matrix,
     file = "3_data_analysis/11_informative_go_benchmarking/wang_sim_matrix.rda")
save(wang_sim_df,
     file = "3_data_analysis/11_informative_go_benchmarking/wang_sim_df.rda")

# 2. Add to the combined similarity table ====
load("3_data_analysis/11_informative_go_benchmarking/go_combined_sim_df.rda")

go_combined_sim_df <-
  go_combined_sim_df |>
  dplyr::select(-dplyr::any_of("sim_wang")) |>
  dplyr::left_join(wang_sim_df |> dplyr::rename(sim_wang = sim),
                   by = c("from", "to"))

message("Pairs with NA Wang similarity: ", sum(is.na(go_combined_sim_df$sim_wang)))
save(go_combined_sim_df,
     file = "3_data_analysis/11_informative_go_benchmarking/go_combined_sim_df.rda")

# 3. Plot ====
plot_same_diff <- function(data, sim_col, y_lab) {
  set.seed(123)
  data$x_jit <- as.numeric(data$class) - 0.1 - runif(nrow(data), 0, 0.2)

  data |>
    ggplot(aes(x = class, y = .data[[sim_col]])) +
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
    geom_point(
      data = data,
      aes(x = x_jit, y = .data[[sim_col]], fill = class),
      inherit.aes = FALSE,
      color = "black",
      size = 1.5,
      shape = 21,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    scale_fill_manual(values = same_different_module_color) +
    scale_color_manual(values = same_different_module_color) +
    labs(x = "", y = y_lab) +
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
}

plot_data <- go_combined_sim_df |> dplyr::filter(!is.na(sim_wang))

plot <- plot_same_diff(plot_data, "sim_wang", "Wang")
plot
ggsave(
  plot = plot,
  filename = "3_data_analysis/11_informative_go_benchmarking/wang_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)

# Wilcoxon test and group medians
print(wilcox.test(sim_wang ~ class, data = plot_data))

plot_data |>
  dplyr::group_by(class) |>
  dplyr::summarise(
    median_wang = median(sim_wang),
    n = dplyr::n()
  ) |>
  print()
