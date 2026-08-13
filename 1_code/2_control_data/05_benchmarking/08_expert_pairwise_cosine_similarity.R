## Pairwise cosine similarity among expert annotations within each cluster

input_file <- file.path(
  "3_data_analysis", "02_control_data", "05_benchmarking",
  "comparison_result", "all_result_long_data.rda"
)
output_dir <- dirname(input_file)

if (!file.exists(input_file)) {
  stop(
    "Input file not found. Run this script from the project root: ",
    input_file
  )
}

if (!requireNamespace("dplyr", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE)) {
  stop("This script requires the dplyr and ggplot2 packages.")
}

input_env <- new.env(parent = emptyenv())
load(input_file, envir = input_env)

if (!exists("all_result_long_data", envir = input_env, inherits = FALSE)) {
  stop("all_result_long_data was not found in: ", input_file)
}

all_result_long_data <- input_env$all_result_long_data
required_columns <- c("cluster", "method", "label", "combined_emb")
missing_columns <- setdiff(required_columns, names(all_result_long_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

tool_methods <- c(
  "mapa_cluster_label",
  "apear_cluster_label",
  "paver_cluster_label",
  "enrich_plot_cluster_label"
)

expert_data <- all_result_long_data |>
  dplyr::filter(!method %in% tool_methods) |>
  dplyr::select(cluster, expert_id = method, label, combined_emb) |>
  dplyr::arrange(cluster, expert_id)

if (anyDuplicated(expert_data[c("cluster", "expert_id")])) {
  stop("Each cluster-expert combination must occur exactly once.")
}

embedding_lengths <- lengths(expert_data$combined_emb)
if (any(embedding_lengths == 0L) || length(unique(embedding_lengths)) != 1L) {
  stop("All combined_emb values must be non-empty vectors of equal length.")
}

cosine_similarity <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  denominator <- sqrt(sum(x^2)) * sqrt(sum(y^2))

  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }

  sum(x * y) / denominator
}

calculate_cluster_pairs <- function(cluster_data) {
  if (nrow(cluster_data) < 2L) {
    return(NULL)
  }

  pair_indices <- utils::combn(seq_len(nrow(cluster_data)), 2L)

  data.frame(
    cluster = cluster_data$cluster[1L],
    expert_id_1 = cluster_data$expert_id[pair_indices[1L, ]],
    expert_id_2 = cluster_data$expert_id[pair_indices[2L, ]],
    label_1 = cluster_data$label[pair_indices[1L, ]],
    label_2 = cluster_data$label[pair_indices[2L, ]],
    cosine_similarity = vapply(
      seq_len(ncol(pair_indices)),
      function(i) cosine_similarity(
        cluster_data$combined_emb[[pair_indices[1L, i]]],
        cluster_data$combined_emb[[pair_indices[2L, i]]]
      ),
      numeric(1L)
    ),
    stringsAsFactors = FALSE
  )
}

expert_pairwise_cosine <- split(expert_data, expert_data$cluster) |>
  lapply(calculate_cluster_pairs) |>
  do.call(what = rbind) |>
  dplyr::as_tibble() |>
  dplyr::mutate(
    involves_not_available =
      trimws(label_1) == "Not Available" |
      trimws(label_2) == "Not Available",
    reporting_status = dplyr::if_else(
      involves_not_available,
      "Involves 'Not Available' (excluded from pooled summary)",
      "Included in pooled summary"
    )
  ) |>
  dplyr::arrange(cluster, expert_id_1, expert_id_2)

if (anyNA(expert_pairwise_cosine$cosine_similarity)) {
  stop("At least one cosine similarity could not be calculated.")
}

summarise_similarity <- function(data) {
  data |>
    dplyr::summarise(
      n_pairs = dplyr::n(),
      median = stats::median(cosine_similarity),
      q1 = as.numeric(stats::quantile(cosine_similarity, 0.25)),
      q3 = as.numeric(stats::quantile(cosine_similarity, 0.75)),
      .groups = "drop"
    )
}

cluster_similarity_summary <- expert_pairwise_cosine |>
  dplyr::group_by(cluster) |>
  summarise_similarity() |>
  dplyr::rename(
    n_pairs_all = n_pairs,
    median_all = median,
    q1_all = q1,
    q3_all = q3
  ) |>
  dplyr::left_join(
    expert_pairwise_cosine |>
      dplyr::filter(!involves_not_available) |>
      dplyr::group_by(cluster) |>
      summarise_similarity() |>
      dplyr::rename(
        n_pairs_reporting = n_pairs,
        median_reporting = median,
        q1_reporting = q1,
        q3_reporting = q3
      ),
    by = "cluster"
  )

pooled_similarity_summary <- expert_pairwise_cosine |>
  dplyr::filter(!involves_not_available) |>
  summarise_similarity() |>
  dplyr::mutate(
    excluded_pairs = sum(expert_pairwise_cosine$involves_not_available),
    report = sprintf(
      "a median pairwise cosine similarity of %.3f (IQR: %.3f\u2013%.3f)",
      median, q1, q3
    )
  )

cluster_levels <- sort(unique(expert_pairwise_cosine$cluster))
plot_data <- expert_pairwise_cosine |>
  dplyr::mutate(cluster = factor(cluster, levels = cluster_levels))

set.seed(20260805)
expert_similarity_boxplot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = cluster, y = cosine_similarity)
) +
  ggplot2::geom_boxplot(
    width = 0.62,
    fill = "#DCE6F1",
    colour = "#2F3E4E",
    linewidth = 0.45,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    data = dplyr::filter(plot_data, !involves_not_available),
    width = 0.13,
    height = 0,
    shape = 21,
    size = 2.2,
    stroke = 0.35,
    colour = "#2F3E4E",
    fill = "#4C78A8",
    alpha = 0.85
  ) +
  ggplot2::geom_jitter(
    data = dplyr::filter(plot_data, involves_not_available),
    ggplot2::aes(
      colour = "Involves Expert #1's 'Not Available' annotation"
    ),
    width = 0.13,
    height = 0,
    shape = 4,
    size = 3.1,
    stroke = 1
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c(
      "Involves Expert #1's 'Not Available' annotation" = "#D55E00"
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = ggplot2::expansion(mult = c(0.02, 0.04))
  ) +
  ggplot2::labs(
    x = "Cluster",
    y = "Pairwise cosine similarity",
    caption = paste0(
      "All expert pairs are shown. The four cluster 5 pairs involving Expert #1's ",
      "'Not Available' annotation\nare marked in orange and excluded only from ",
      "the pooled median and IQR."
    )
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.justification = "left",
    plot.caption = ggplot2::element_text(
      size = 7.5,
      hjust = 0,
      colour = "#444444"
    )
  )

plot_summary <- cluster_similarity_summary |>
  dplyr::mutate(cluster = factor(cluster, levels = cluster_levels))

set.seed(20260805)
expert_similarity_barplot <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(x = cluster, y = median_all)
) +
  ggplot2::geom_col(
    width = 0.62,
    fill = "#DCE6F1",
    colour = "#2F3E4E",
    linewidth = 0.45
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = q1_all, ymax = q3_all),
    width = 0.18,
    colour = "#2F3E4E",
    linewidth = 0.55
  ) +
  ggplot2::geom_jitter(
    data = dplyr::filter(plot_data, !involves_not_available),
    ggplot2::aes(x = cluster, y = cosine_similarity),
    inherit.aes = FALSE,
    width = 0.13,
    height = 0,
    shape = 21,
    size = 2.2,
    stroke = 0.35,
    colour = "#2F3E4E",
    fill = "#4C78A8",
    alpha = 0.85
  ) +
  ggplot2::geom_jitter(
    data = dplyr::filter(plot_data, involves_not_available),
    ggplot2::aes(
      x = cluster,
      y = cosine_similarity,
      colour = "Involves Expert #1's 'Not Available' annotation"
    ),
    inherit.aes = FALSE,
    width = 0.13,
    height = 0,
    shape = 4,
    size = 3.1,
    stroke = 1
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c(
      "Involves Expert #1's 'Not Available' annotation" = "#D55E00"
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = ggplot2::expansion(mult = c(0, 0.04))
  ) +
  ggplot2::labs(
    x = "Cluster",
    y = "Pairwise cosine similarity",
    caption = paste0(
      "Bars and error bars show the median and IQR of all plotted pairs, respectively. ",
      "The four cluster 5 pairs involving Expert #1's\n'Not Available' annotation ",
      "are marked in orange and excluded only from the pooled median and IQR."
    )
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.justification = "left",
    plot.caption = ggplot2::element_text(
      size = 7.5,
      hjust = 0,
      colour = "#444444"
    )
  )

utils::write.csv(
  expert_pairwise_cosine,
  file.path(output_dir, "expert_pairwise_cosine_similarity.csv"),
  row.names = FALSE
)
utils::write.csv(
  cluster_similarity_summary,
  file.path(output_dir, "expert_pairwise_cosine_similarity_by_cluster.csv"),
  row.names = FALSE
)
utils::write.csv(
  pooled_similarity_summary,
  file.path(output_dir, "expert_pairwise_cosine_similarity_pooled_summary.csv"),
  row.names = FALSE
)

save(
  expert_pairwise_cosine,
  cluster_similarity_summary,
  pooled_similarity_summary,
  expert_similarity_boxplot,
  expert_similarity_barplot,
  file = file.path(output_dir, "expert_pairwise_cosine_similarity.rda")
)

ggplot2::ggsave(
  filename = file.path(output_dir, "expert_pairwise_cosine_similarity_boxplot.pdf"),
  plot = expert_similarity_boxplot,
  width = 7.2,
  height = 4.8,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(output_dir, "expert_pairwise_cosine_similarity_boxplot.png"),
  plot = expert_similarity_boxplot,
  width = 7.2,
  height = 4.8,
  units = "in",
  dpi = 300
)
ggplot2::ggsave(
  filename = file.path(output_dir, "expert_pairwise_cosine_similarity_barplot.pdf"),
  plot = expert_similarity_barplot,
  width = 7.2,
  height = 4.8,
  units = "in"
)
ggplot2::ggsave(
  filename = file.path(output_dir, "expert_pairwise_cosine_similarity_barplot.png"),
  plot = expert_similarity_barplot,
  width = 7.2,
  height = 4.8,
  units = "in",
  dpi = 300
)

message(pooled_similarity_summary$report)
message(
  "Reporting pairs: ", pooled_similarity_summary$n_pairs,
  "; excluded 'Not Available' pairs: ", pooled_similarity_summary$excluded_pairs
)
