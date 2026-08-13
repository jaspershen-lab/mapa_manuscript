## Recreate the native network visualizations for the optimized benchmark
## settings on the curated 44-pathway dataset.

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this file with Rscript so the manuscript root can be located.")
}
script_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/")
root <- normalizePath(file.path(dirname(script_file), "../../.."), winslash = "/")

input_dir <- file.path(root, "3_data_analysis/02_control_data/05_benchmarking")
output_dir <- file.path(
  root,
  "3_data_analysis/12_benchmark_clustering/01_curated_pathway_dataset"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(aPEAR)
  library(dplyr)
  library(enrichplot)
  library(ggplot2)
  library(mclust)
  library(patchwork)
  library(PAVER)
  library(readr)
})

load_object <- function(path, object_name) {
  env <- new.env(parent = globalenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Object '", object_name, "' was not found in ", path)
  }
  get(object_name, envir = env, inherits = FALSE)
}

assert_same_partition <- function(observed, expected, tool) {
  if (length(observed) != length(expected) || anyNA(observed) || anyNA(expected)) {
    stop(tool, " plot labels are incomplete or have an unexpected length.")
  }
  agreement <- mclust::adjustedRandIndex(as.character(observed), expected)
  if (!isTRUE(all.equal(agreement, 1, tolerance = 1e-12))) {
    stop(tool, " plot labels do not reproduce the optimized assignments (ARI = ",
         signif(agreement, 6), ").")
  }
  invisible(agreement)
}

assignments <- readr::read_csv(
  file.path(output_dir, "optimized_clustering_assignments.csv"),
  show_col_types = FALSE
)
enriched_result <- load_object(
  file.path(input_dir, "enriched_result.rda"), "enriched_result"
)

## PAVER ---------------------------------------------------------------------
paver <- load_object(
  file.path(output_dir, "paver_sensitivity.rda"), "paver"
)
paver_best <- paver$best
paver_env <- new.env(parent = globalenv())
load(file.path(input_dir, "PAVER_01_input.rda"), envir = paver_env)
load(file.path(input_dir, "PAVER_02_embedding_matrix.rda"), envir = paver_env)
load(file.path(input_dir, "PAVER_03_term2name.rda"), envir = paver_env)

paver_result <- PAVER::prepare_data(
  paver_env$input,
  paver_env$embedding_matrix,
  paver_env$term2name
)
set.seed(234)
paver_result <- PAVER::generate_themes(
  paver_result,
  hclust_method = paver_best$hclust_method[[1]],
  maxCoreScatter = paver_best$maxCoreScatter[[1]],
  minGap = paver_best$minGap[[1]],
  minClusterSize = paver_best$minClusterSize[[1]]
)

paver_export <- PAVER::PAVER_export(paver_result)
paver_expected <- assignments$PAVER[match(paver_export$GOID, assignments$ID)]
assert_same_partition(paver_export$Cluster, paver_expected, "PAVER")

paver_plot <- PAVER::PAVER_interpretation_plot(paver_result) +
  PAVER::PAVER_theme_plot(paver_result)
ggplot2::ggsave(
  filename = file.path(output_dir, "paver_best_clustering_network.pdf"),
  plot = paver_plot,
  width = 12,
  height = 6
)

## aPEAR ---------------------------------------------------------------------
apear <- load_object(
  file.path(output_dir, "apear_sensitivity.rda"), "apear"
)
apear_best <- apear$best
set.seed(135)
apear_result <- aPEAR::enrichmentNetwork(
  enrichment = enriched_result@result,
  simMethod = apear_best$simMethod[[1]],
  innerCutoff = 0.1,
  outerCutoff = 0.5,
  clustMethod = apear_best$clustMethod[[1]],
  clustNameMethod = "pagerank",
  colorBy = "pvalue",
  colorType = "pval",
  nodeSize = "Count",
  plotOnly = FALSE,
  drawEllipses = TRUE,
  minClusterSize = apear_best$minClusterSize[[1]],
  fontSize = 3
)

apear_map <- enriched_result@result |>
  dplyr::select(ID, Description)
apear_validation <- as.data.frame(apear_result$clusters) |>
  dplyr::rename(Description = ID, observed = Cluster) |>
  dplyr::left_join(apear_map, by = "Description")
apear_expected <- assignments$aPEAR[
  match(apear_validation$ID, assignments$ID)
]
assert_same_partition(apear_validation$observed, apear_expected, "aPEAR")

ggplot2::ggsave(
  filename = file.path(output_dir, "apear_best_clustering_network.pdf"),
  plot = apear_result$plot,
  width = 10,
  height = 8
)

## enrichplot ----------------------------------------------------------------
enrichplot_result <- load_object(
  file.path(output_dir, "enrichplot_sensitivity.rda"), "enrichplot_res"
)
enrichplot_best <- enrichplot_result$best
n_terms <- nrow(enriched_result@result)
edo <- enrichplot::pairwise_termsim(
  enriched_result,
  method = "JC",
  showCategory = n_terms
)

set.seed(123)
enrichplot_plot <- enrichplot::emapplot(
  x = edo,
  showCategory = n_terms,
  layout = enrichplot_best$layout[[1]],
  min_edge = enrichplot_best$min_edge[[1]],
  color_edge = "grey",
  node_label = "group",
  group = TRUE,
  clusterFunction = stats::kmeans,
  nCluster = enrichplot_best$nCluster[[1]]
)

enrichplot_map <- edo@result |>
  dplyr::select(ID, Description)
enrichplot_validation <- enrichplot_plot$data |>
  dplyr::transmute(Description = name, observed = color2) |>
  dplyr::left_join(enrichplot_map, by = "Description")
enrichplot_expected <- assignments$enrichplot[
  match(enrichplot_validation$ID, assignments$ID)
]
assert_same_partition(
  enrichplot_validation$observed,
  enrichplot_expected,
  "enrichplot"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "enrichplot_best_clustering_network.pdf"),
  plot = enrichplot_plot,
  width = 10,
  height = 8
)

message("Generated optimized PAVER, aPEAR, and enrichplot network PDFs in: ",
        output_dir)
