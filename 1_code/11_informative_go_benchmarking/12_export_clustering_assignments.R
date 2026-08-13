## Benchmarking on the informative GO set benchmark (k = 20)
## Export term-level cluster assignments for the default and optimized
## configurations of MAPA, PAVER, aPEAR, and enrichplot.

library(r4projects)
setwd(get_project_wd())
rm(list = ls())

suppressPackageStartupMessages({
  library(mapa)
  library(dplyr)
  library(readr)
  library(mclust)
  library(enrichplot)
  library(aPEAR)
  library(PAVER)
})

res_dir <- "3_data_analysis/11_informative_go_benchmarking"

load(file.path(res_dir, "ground_truth_dt.rda"))
load(file.path(res_dir, "enriched_result.rda"))
load(file.path(res_dir, "default_vs_optimized_ari.rda"))

stopifnot(
  nrow(ground_truth_dt) == 313L,
  !anyDuplicated(ground_truth_dt$id),
  !anyNA(ground_truth_dt$ground_truth_label)
)

as_sequential_integer <- function(x) {
  as.integer(match(x, unique(x)))
}

assignment_frame <- function(id, label, label_name) {
  out <- data.frame(
    id = as.character(id),
    label = as_sequential_integer(label),
    stringsAsFactors = FALSE
  )
  names(out)[2] <- label_name
  stopifnot(!anyDuplicated(out$id), !anyNA(out[[label_name]]))
  out
}

join_assignments <- function(...) {
  pieces <- list(...)
  out <- ground_truth_dt
  for (piece in pieces) {
    out <- dplyr::left_join(out, piece, by = "id")
  }
  stopifnot(nrow(out) == nrow(ground_truth_dt), !anyNA(out))
  out |>
    dplyr::rename(ID = id) |>
    dplyr::arrange(mapa_label, ID)
}

check_assignment <- function(table, method, configuration, label_column) {
  expected <- default_vs_optimized_ari |>
    dplyr::filter(methods == method)

  expected_ari <- if (configuration == "default") {
    expected$default_ari
  } else {
    expected$optimized_ari
  }
  expected_n <- if (configuration == "default") {
    expected$default_n_clusters
  } else {
    expected$optimized_n_clusters
  }

  observed_ari <- mclust::adjustedRandIndex(
    table[[label_column]], table$ground_truth_label
  )
  observed_n <- dplyr::n_distinct(table[[label_column]])

  if (!isTRUE(all.equal(observed_ari, expected_ari, tolerance = 1e-10)) ||
      observed_n != expected_n) {
    stop(
      method, " ", configuration, " assignment validation failed: ",
      "observed ARI=", signif(observed_ari, 12),
      ", expected ARI=", signif(expected_ari, 12),
      ", observed clusters=", observed_n,
      ", expected clusters=", expected_n
    )
  }

  data.frame(
    configuration = configuration,
    method = method,
    ari = observed_ari,
    n_clusters = observed_n,
    stringsAsFactors = FALSE
  )
}

# MAPA: default ------------------------------------------------------------
load("3_data_analysis/10_sensitivity_res_GO/bioembed_sim_res.rda")
load(file.path(res_dir, "pathway_genes.rda"))

go_result <- bioembed_sim_res$enriched_pathway@enrichment_go_result@result
go_result$geneID <- vapply(
  pathway_genes[go_result$ID], paste, collapse = "/", FUN.VALUE = character(1)
)
bioembed_sim_res$enriched_pathway@enrichment_go_result@result <- go_result

all_accessions <- unique(unlist(pathway_genes, use.names = FALSE))
bioembed_sim_res$enriched_pathway@variable_info <- data.frame(
  ensembl = all_accessions,
  entrezid = all_accessions,
  uniprot = all_accessions,
  symbol = all_accessions,
  variable_id = paste0("gene_", seq_along(all_accessions))
)
bioembed_sim_res$enriched_pathway@enrichment_kegg_result <- NULL
bioembed_sim_res$enriched_pathway@enrichment_reactome_result <- NULL

set.seed(23)
mapa_default_result <- get_functional_modules(
  object = bioembed_sim_res,
  sim.cutoff = 0.55,
  cluster_method = "louvain"
)
mapa_default_raw <- mapa_default_result@merged_module$result_with_module |>
  dplyr::select(id = node, label = module) |>
  as.data.frame()
mapa_default <- assignment_frame(
  mapa_default_raw$id,
  mapa_default_raw$label,
  "mapa_label"
)

# MAPA: optimized ----------------------------------------------------------
load("3_data_analysis/10_sensitivity_res_GO/embedding_sim_matrix.rda")
norms <- sqrt(rowSums(embedding_sim_matrix^2))
normalized_sim <- embedding_sim_matrix / norms

set.seed(23)
mapa_optimized_fit <- stats::kmeans(
  normalized_sim,
  centers = 77,
  nstart = 50
)
mapa_optimized <- assignment_frame(
  rownames(normalized_sim),
  mapa_optimized_fit$cluster,
  "mapa_label"
)

# enrichplot: default and optimized ---------------------------------------
n_terms <- nrow(enriched_result@result)
edo <- enrichplot::pairwise_termsim(enriched_result, showCategory = n_terms)

extract_enrichplot_assignment <- function(n_cluster = NULL) {
  set.seed(123)
  plot_args <- list(
    x = edo,
    showCategory = n_terms,
    node_label = "group",
    group = TRUE
  )
  if (!is.null(n_cluster)) {
    plot_args$nCluster <- n_cluster
    plot_args$clusterFunction <- stats::kmeans
    plot_args$layout <- "fr"
    plot_args$min_edge <- 0.2
    plot_args$color_edge <- "grey"
  }
  p <- do.call(enrichplot::emapplot, plot_args)

  joined <- p$data[, c("name", "color2")] |>
    dplyr::rename(Description = name, label = color2) |>
    dplyr::left_join(
      edo@result[, c("ID", "Description")],
      by = "Description"
    )
  stopifnot(!anyNA(joined$ID), !anyDuplicated(joined$ID))
  assignment_frame(joined$ID, joined$label, "enrichplot_cluster")
}

enrichplot_default <- extract_enrichplot_assignment()
enrichplot_optimized <- extract_enrichplot_assignment(73L)

# aPEAR: default (markov) and optimized (hier) ----------------------------
extract_apear_assignment <- function(apear_result) {
  apear_clusters <- apear_result$clusters |>
    dplyr::rename(Description = ID)

  pathway_clusters <- enriched_result@result |>
    dplyr::select(ID, Description) |>
    dplyr::left_join(apear_clusters, by = "Description")

  non_singleton_labels <- unique(apear_clusters$Cluster)
  pathway_clusters$cluster_label <- match(
    pathway_clusters$Cluster,
    non_singleton_labels
  )

  singleton_rows <- which(is.na(pathway_clusters$cluster_label))
  if (length(singleton_rows) > 0L) {
    max_label <- max(pathway_clusters$cluster_label, na.rm = TRUE)
    pathway_clusters$cluster_label[singleton_rows] <-
      seq.int(max_label + 1L, length.out = length(singleton_rows))
  }

  assignment_frame(
    pathway_clusters$ID,
    pathway_clusters$cluster_label,
    "apear_cluster"
  )
}

apear_results <- list()
set.seed(135)
for (method in c("markov", "hier")) {
  apear_results[[method]] <- aPEAR::enrichmentNetwork(
    enrichment = enriched_result@result,
    simMethod = "jaccard",
    innerCutoff = 0.1,
    outerCutoff = 0.5,
    clustMethod = method,
    clustNameMethod = "pagerank",
    colorBy = "pvalue",
    colorType = "pval",
    nodeSize = "Count",
    plotOnly = FALSE,
    drawEllipses = TRUE,
    minClusterSize = 2,
    fontSize = 3
  )
}
apear_default <- extract_apear_assignment(apear_results$markov)
apear_optimized <- extract_apear_assignment(apear_results$hier)

# PAVER: default and optimized --------------------------------------------
load(file.path(res_dir, "PAVER_01_input.rda"))
load(file.path(res_dir, "PAVER_02_embedding_matrix.rda"))
load(file.path(res_dir, "PAVER_03_term2name.rda"))
PAVER_result <- PAVER::prepare_data(input, embedding_matrix, term2name)

extract_paver_assignment <- function(paver_result) {
  exported <- PAVER::PAVER_export(paver_result)
  assignment_frame(exported$GOID, exported$Cluster, "paver_label")
}

set.seed(234)
paver_default <- extract_paver_assignment(
  PAVER::generate_themes(PAVER_result)
)

set.seed(234)
paver_optimized <- extract_paver_assignment(
  PAVER::generate_themes(
    PAVER_result,
    maxCoreScatter = 0.88,
    minGap = (1 - 0.88) * 3 / 4,
    minClusterSize = 1
  )
)

# Merge, validate, and export ---------------------------------------------
default_assignments <- join_assignments(
  mapa_default,
  paver_default,
  apear_default,
  enrichplot_default
)
optimized_assignments <- join_assignments(
  mapa_optimized,
  paver_optimized,
  apear_optimized,
  enrichplot_optimized
)

combined_assignments <- default_assignments |>
  dplyr::rename(
    mapa_default_label = mapa_label,
    paver_default_label = paver_label,
    apear_default_cluster = apear_cluster,
    enrichplot_default_cluster = enrichplot_cluster
  ) |>
  dplyr::left_join(
    optimized_assignments |>
      dplyr::select(-ground_truth_label) |>
      dplyr::rename(
        mapa_optimized_label = mapa_label,
        paver_optimized_label = paver_label,
        apear_optimized_cluster = apear_cluster,
        enrichplot_optimized_cluster = enrichplot_cluster
      ),
    by = "ID"
  ) |>
  dplyr::select(
    ID,
    ground_truth_label,
    mapa_default_label,
    mapa_optimized_label,
    paver_default_label,
    paver_optimized_label,
    apear_default_cluster,
    apear_optimized_cluster,
    enrichplot_default_cluster,
    enrichplot_optimized_cluster
  ) |>
  dplyr::arrange(mapa_optimized_label, ID)

stopifnot(
  nrow(combined_assignments) == nrow(ground_truth_dt),
  !anyDuplicated(combined_assignments$ID),
  !anyNA(combined_assignments)
)

validation <- dplyr::bind_rows(
  check_assignment(default_assignments, "MAPA", "default", "mapa_label"),
  check_assignment(default_assignments, "PAVER", "default", "paver_label"),
  check_assignment(default_assignments, "aPEAR", "default", "apear_cluster"),
  check_assignment(default_assignments, "enrichplot", "default", "enrichplot_cluster"),
  check_assignment(optimized_assignments, "MAPA", "optimized", "mapa_label"),
  check_assignment(optimized_assignments, "PAVER", "optimized", "paver_label"),
  check_assignment(optimized_assignments, "aPEAR", "optimized", "apear_cluster"),
  check_assignment(optimized_assignments, "enrichplot", "optimized", "enrichplot_cluster")
)

readr::write_csv(
  default_assignments,
  file.path(res_dir, "default_clustering_assignments.csv")
)
readr::write_csv(
  optimized_assignments,
  file.path(res_dir, "optimized_clustering_assignments.csv")
)
readr::write_csv(
  combined_assignments,
  file.path(res_dir, "default_optimized_clustering_assignments.csv")
)
readr::write_csv(
  validation,
  file.path(res_dir, "clustering_assignments_validation.csv")
)

save(
  default_assignments,
  optimized_assignments,
  combined_assignments,
  validation,
  file = file.path(res_dir, "default_optimized_clustering_assignments.rda")
)

cat("\nValidated term-level assignments:\n")
print(validation, row.names = FALSE)
cat("\nRows per table:", nrow(default_assignments), "\n")
