## Benchmarking on the informative GO set benchmark (k = 20)
## Step 7: Jaccard / overlap-coefficient similarity between GO terms in the
## same vs different ground-truth modules.
## Similarity calculation follows
## 1_code/2_control_data/02_bioembedsim_vs_othersim/02_overlap_sim_calculation.R
## (term_similarity_internal from utils.R, gene sets = UniProt accessions from
## module_gene_sets.tsv) and the plots follow the style of
## 1_code/2_control_data/02_bioembedsim_vs_othersim/05_compare_bioembed_vs_go_terms.R
## (point size reduced: 48,828 term pairs vs 231 in the control data)

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

library(gghalves)
library(ggsignif)

# 1. Gene sets and ground truth ====
gene_sets <- readr::read_tsv("2_data/informative_go_set/k_20/module_gene_sets.tsv",
                             show_col_types = FALSE)
ground_truth_dt <- readr::read_csv("2_data/informative_go_set/k_20/ground_truth.csv",
                                   show_col_types = FALSE)

gene_lists <- strsplit(gene_sets$uniprot_accessions, ";")
names(gene_lists) <- gene_sets$go_id

# 2. Similarity based on gene overlap ====
jaccard_sim_matrix <- term_similarity_internal(gl = gene_lists,
                                               measure.method = "jaccard")
op_sim_matrix <- term_similarity_internal(gl = gene_lists,
                                          measure.method = "overlap")

matrix_to_df <- function(m) {
  as.data.frame.table(m, responseName = "sim") |>
    dplyr::filter(Var1 != Var2) |>
    dplyr::rename(from = Var1, to = Var2) |>
    dplyr::mutate(across(c(from, to), as.character)) |>
    dplyr::filter(from < to)
}

jaccard_sim_df <- matrix_to_df(jaccard_sim_matrix)
op_sim_df <- matrix_to_df(op_sim_matrix)

save(jaccard_sim_df,
     file = "3_data_analysis/11_informative_go_benchmarking/jaccard_sim_df.rda")
save(op_sim_df,
     file = "3_data_analysis/11_informative_go_benchmarking/op_sim_df.rda")

# 3. Combine with the biotext embedding similarity and label pairs ====
## Biotext embedding similarity computed by mapa::get_bioembedsim()
## (1_code/10_sensitivity_analysis_GO/01_build_object_and_embedding.R)
load("3_data_analysis/10_sensitivity_res_GO/embedding_sim_df.rda")

go_combined_sim_df <-
  embedding_sim_df |>
  dplyr::rename(sim_embedding = sim) |>
  dplyr::left_join(jaccard_sim_df |> dplyr::rename(sim_jaccard = sim),
                   by = c("from", "to")) |>
  dplyr::left_join(op_sim_df |> dplyr::rename(sim_op = sim),
                   by = c("from", "to")) |>
  dplyr::left_join(ground_truth_dt[, c("go_id", "module_id")],
                   by = c("from" = "go_id")) |>
  dplyr::left_join(ground_truth_dt[, c("go_id", "module_id")],
                   by = c("to" = "go_id")) |>
  dplyr::mutate(
    class = case_when(
      module_id.x == module_id.y ~ "same_module",
      module_id.x != module_id.y ~ "different_module"
    )
  ) |>
  dplyr::mutate(class = factor(class, levels = c("same_module", "different_module")))

stopifnot(!anyNA(go_combined_sim_df$sim_jaccard),
          !anyNA(go_combined_sim_df$sim_op))
print(table(go_combined_sim_df$class))
print(summary(go_combined_sim_df$sim_embedding))
save(go_combined_sim_df,
     file = "3_data_analysis/11_informative_go_benchmarking/go_combined_sim_df.rda")

# 4. Plots ====
## Note: geom_half_point() is broken with the currently installed ggplot2
## (error "argument 'layout' is missing"), so the left-side jittered points
## are drawn with geom_point() on manually jittered x positions instead
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
    ## Jittered points added after the discrete-x layers so that the numeric
    ## x positions do not change the x-scale type
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

plot <- plot_same_diff(go_combined_sim_df, "sim_jaccard", "Jaccard index")
plot
ggsave(
  plot = plot,
  filename = "3_data_analysis/11_informative_go_benchmarking/jaccard_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)

plot <- plot_same_diff(go_combined_sim_df, "sim_op", "Overlap coefficient")
plot
ggsave(
  plot = plot,
  filename = "3_data_analysis/11_informative_go_benchmarking/op_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)

plot <- plot_same_diff(go_combined_sim_df, "sim_embedding", "PathBiotext")
plot
ggsave(
  plot = plot,
  filename = "3_data_analysis/11_informative_go_benchmarking/pathbiotext_same_different_module_comparison.pdf",
  width = 5,
  height = 5
)

# Wilcoxon tests (exact statistics)
print(wilcox.test(sim_jaccard ~ class, data = go_combined_sim_df))
print(wilcox.test(sim_op ~ class, data = go_combined_sim_df))
print(wilcox.test(sim_embedding ~ class, data = go_combined_sim_df))

## Group medians
go_combined_sim_df |>
  dplyr::group_by(class) |>
  dplyr::summarise(
    median_jaccard = median(sim_jaccard),
    median_op = median(sim_op),
    median_embedding = median(sim_embedding),
    n = dplyr::n()
  ) |>
  print()
