## Revised benchmarking: curated 44-pathway control dataset.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")
script_dir <- dirname(script_file)
options(mapa_manuscript_root = normalizePath(file.path(script_dir, "../../.."),
                                             winslash = "/"))
source(file.path(dirname(script_dir), "benchmark_functions.R"))

get_curated_mapa_results <- function(root, ground_truth_dt) {
  suppressPackageStartupMessages({
    library(tidygraph)
    library(dplyr)
  })
  graph <- load_single_object(
    file.path(root, "3_data_analysis/02_control_data/04_clustering/louvain_graph_data.rda"),
    "louvain_graph_data"
  )
  node_data <- graph |>
    tidygraph::activate("nodes") |>
    tibble::as_tibble()
  labels <- node_data$louvain_result[match(ground_truth_dt$id, node_data$node)]
  if (anyNA(labels)) stop("MAPA curated assignments are incomplete")

  list(
    summary = data.frame(
      method = "MAPA",
      default_ari = mclust::adjustedRandIndex(labels, ground_truth_dt$ground_truth_label),
      default_n_clusters = dplyr::n_distinct(labels),
      default_native_n_clusters = dplyr::n_distinct(labels),
      default_n_unassigned = 0L,
      default_settings = "sim.cutoff=0.55; cluster_method=louvain",
      optimized_ari = mclust::adjustedRandIndex(labels, ground_truth_dt$ground_truth_label),
      optimized_n_clusters = dplyr::n_distinct(labels),
      optimized_native_n_clusters = dplyr::n_distinct(labels),
      optimized_n_unassigned = 0L,
      optimized_settings = "sim.cutoff=0.55; cluster_method=louvain (existing result)",
      default_status = "ok"
    ),
    optimized_labels = labels,
    default_labels = labels
  )
}

run_complete_benchmark(list(
  dataset_name = "Curated pathway dataset",
  input_dir = "3_data_analysis/02_control_data/05_benchmarking",
  output_dir = "3_data_analysis/12_benchmark_clustering/01_curated_pathway_dataset",
  get_mapa_results = get_curated_mapa_results
))
