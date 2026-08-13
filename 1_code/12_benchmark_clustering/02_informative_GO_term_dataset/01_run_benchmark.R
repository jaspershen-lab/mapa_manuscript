## Revised benchmarking: informative GO-term dataset.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")
script_dir <- dirname(script_file)
options(mapa_manuscript_root = normalizePath(file.path(script_dir, "../../.."),
                                             winslash = "/"))
source(file.path(dirname(script_dir), "benchmark_functions.R"))

get_informative_mapa_results <- function(root, ground_truth_dt) {
  summary_path <- file.path(
    root,
    "3_data_analysis/11_informative_go_benchmarking/default_vs_optimized_ari.csv"
  )
  assignment_path <- file.path(
    root,
    "3_data_analysis/11_informative_go_benchmarking/default_optimized_clustering_assignments.csv"
  )
  existing_summary <- readr::read_csv(summary_path, show_col_types = FALSE) |>
    dplyr::filter(methods == "MAPA")
  existing_assignments <- readr::read_csv(assignment_path, show_col_types = FALSE)
  labels <- existing_assignments$mapa_optimized_label[
    match(ground_truth_dt$id, existing_assignments$ID)
  ]
  default_labels <- existing_assignments$mapa_default_label[
    match(ground_truth_dt$id, existing_assignments$ID)
  ]
  if (anyNA(labels) || anyNA(default_labels)) {
    stop("MAPA informative-GO assignments are incomplete")
  }

  list(
    summary = data.frame(
      method = "MAPA",
      default_ari = existing_summary$default_ari,
      default_n_clusters = existing_summary$default_n_clusters,
      default_native_n_clusters = existing_summary$default_n_clusters,
      default_n_unassigned = 0L,
      default_settings = existing_summary$default_params,
      optimized_ari = existing_summary$optimized_ari,
      optimized_n_clusters = existing_summary$optimized_n_clusters,
      optimized_native_n_clusters = existing_summary$optimized_n_clusters,
      optimized_n_unassigned = 0L,
      optimized_settings = paste0(existing_summary$optimized_params,
                                  " (existing result)"),
      default_status = "ok"
    ),
    optimized_labels = labels,
    default_labels = default_labels
  )
}

run_complete_benchmark(list(
  dataset_name = "Informative GO-term dataset",
  input_dir = "3_data_analysis/11_informative_go_benchmarking",
  output_dir = "3_data_analysis/12_benchmark_clustering/02_informative_GO_term_dataset",
  get_mapa_results = get_informative_mapa_results
))
