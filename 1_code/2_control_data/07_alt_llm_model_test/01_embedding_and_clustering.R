source("1_code/2_control_data/07_alt_llm_model_test/00_config.R")

required_packages <- c(
  "dplyr", "ggplot2", "igraph", "mclust", "readxl", "stringr", "tibble", "uwot"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

load_one <- function(path, expected_name = NULL) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  if (!is.null(expected_name) && expected_name %in% object_names) {
    return(env[[expected_name]])
  }
  if (length(object_names) != 1L) {
    stop("Expected one object in ", path, "; found: ", paste(object_names, collapse = ", "))
  }
  env[[object_names]]
}

cosine_similarity_matrix <- function(x) {
  x <- as.matrix(x)
  row_norms <- sqrt(rowSums(x^2))
  if (any(!is.finite(row_norms) | row_norms == 0)) {
    stop("Embedding matrix contains a non-finite or zero-length row.")
  }
  x_normalized <- x / row_norms
  result <- tcrossprod(x_normalized)
  dimnames(result) <- list(rownames(x), rownames(x))
  diag(result) <- 1
  result
}

embed_pathway_texts <- function(pathway_text, model_name) {
  if (!requireNamespace("mapa", quietly = TRUE)) {
    stop("The mapa package is required to generate embeddings.")
  }
  api_key <- first_nonempty_env(c("SILICONFLOW_API_KEY", "QWEN_API_KEY"))
  if (!nzchar(api_key)) {
    stop(
      "No cached Qwen embedding matrix was found. Set SILICONFLOW_API_KEY ",
      "(or QWEN_API_KEY) to generate it."
    )
  }
  rows <- lapply(pathway_text, function(item) {
    value <- mapa::get_embedding(
      chunk = item$text_info,
      api_key = api_key,
      model_name = model_name,
      api_provider = "siliconflow"
    )
    if (is.null(value)) stop("Embedding failed for pathway: ", item$id)
    value
  })
  if (length(unique(lengths(rows))) != 1L) {
    stop("Qwen embedding rows have inconsistent dimensions.")
  }
  result <- do.call(rbind, rows)
  rownames(result) <- vapply(pathway_text, `[[`, character(1L), "id")
  result
}

run_default_louvain <- function(similarity_matrix, cutoff, seed) {
  edge_data <- as.data.frame(as.table(similarity_matrix), stringsAsFactors = FALSE)
  names(edge_data) <- c("from", "to", "similarity")
  edge_data <- edge_data |>
    dplyr::mutate(
      from = as.character(from),
      to = as.character(to),
      similarity = as.numeric(similarity)
    ) |>
    dplyr::filter(from < to, similarity >= cutoff)

  vertices <- data.frame(name = rownames(similarity_matrix))
  graph <- igraph::graph_from_data_frame(edge_data, directed = FALSE, vertices = vertices)
  if (nrow(edge_data) == 0L) {
    membership <- seq_len(nrow(vertices))
    names(membership) <- vertices$name
  } else {
    set.seed(seed)
    membership <- igraph::membership(
      igraph::cluster_louvain(graph, weights = igraph::E(graph)$similarity)
    )
  }
  list(graph = graph, edge_data = edge_data, membership = membership[vertices$name])
}

control_data <- readxl::read_excel(control_data_file, sheet = 1)
original_embedding_matrix <- load_one(original_embedding_file, "embedding_matrix")

qwen_embedding_source <- "generated_in_this_run"
if (file.exists(qwen_embedding_cache_file)) {
  qwen_embedding_matrix <- load_one(
    qwen_embedding_cache_file, "siliconflow_embedding_matrix"
  )
  qwen_embedding_source <- paste0(
    "reused_verified_cache:",
    file.path("3_data_analysis", "09_select_model_param", "siliconflow_embedding_matrix.rda")
  )
} else {
  pathway_text <- load_one(pathway_text_file, "all_combined_info")
  qwen_embedding_matrix <- embed_pathway_texts(pathway_text, qwen_embedding_model)
}

if (is.null(rownames(qwen_embedding_matrix)) || is.null(rownames(original_embedding_matrix))) {
  stop("Both embedding matrices must have pathway IDs as row names.")
}
expected_ids <- as.character(control_data$id)
if (!setequal(rownames(qwen_embedding_matrix), expected_ids)) {
  stop("Qwen embedding IDs do not exactly match control_data.xlsx.")
}
if (!setequal(rownames(original_embedding_matrix), expected_ids)) {
  stop("Original embedding IDs do not exactly match control_data.xlsx.")
}
qwen_embedding_matrix <- qwen_embedding_matrix[expected_ids, , drop = FALSE]
original_embedding_matrix <- original_embedding_matrix[expected_ids, , drop = FALSE]

qwen_embedding_sim_matrix <- cosine_similarity_matrix(qwen_embedding_matrix)
original_embedding_sim_matrix <- cosine_similarity_matrix(original_embedding_matrix)
qwen_embedding_sim_df <- as.data.frame(
  as.table(qwen_embedding_sim_matrix), stringsAsFactors = FALSE
) |>
  dplyr::rename(from = Var1, to = Var2, sim = Freq) |>
  dplyr::mutate(from = as.character(from), to = as.character(to)) |>
  dplyr::filter(from < to)

qwen_cluster <- run_default_louvain(
  qwen_embedding_sim_matrix, mapa_similarity_cutoff, analysis_seed
)
original_cluster <- run_default_louvain(
  original_embedding_sim_matrix, mapa_similarity_cutoff, analysis_seed
)

ground_truth <- stats::setNames(
  as.character(control_data$expected_module), as.character(control_data$id)
)[expected_ids]
qwen_membership <- qwen_cluster$membership[expected_ids]
original_membership <- original_cluster$membership[expected_ids]

qwen_ari_ground_truth <- mclust::adjustedRandIndex(ground_truth, qwen_membership)
original_ari_ground_truth <- mclust::adjustedRandIndex(
  ground_truth, original_membership
)
cluster_ari_qwen_vs_original <- mclust::adjustedRandIndex(
  qwen_membership, original_membership
)

upper_index <- upper.tri(qwen_embedding_sim_matrix, diag = FALSE)
qwen_pairwise <- qwen_embedding_sim_matrix[upper_index]
original_pairwise <- original_embedding_sim_matrix[upper_index]
pearson_test <- stats::cor.test(qwen_pairwise, original_pairwise, method = "pearson")
spearman_test <- suppressWarnings(
  stats::cor.test(qwen_pairwise, original_pairwise, method = "spearman", exact = FALSE)
)

embedding_consistency_metrics <- data.frame(
  metric = c(
    "Qwen ARI vs ground truth",
    "Original ARI vs ground truth",
    "Qwen vs original clustering ARI",
    "Pairwise similarity Pearson r",
    "Pairwise similarity Spearman rho",
    "Qwen number of clusters",
    "Original number of clusters",
    "Qwen edges at cutoff",
    "Original edges at cutoff"
  ),
  value = c(
    qwen_ari_ground_truth,
    original_ari_ground_truth,
    cluster_ari_qwen_vs_original,
    unname(pearson_test$estimate),
    unname(spearman_test$estimate),
    length(unique(qwen_membership)),
    length(unique(original_membership)),
    nrow(qwen_cluster$edge_data),
    nrow(original_cluster$edge_data)
  ),
  stringsAsFactors = FALSE
)

default_louvain_clustering_result <- data.frame(
  pathway_id = expected_ids,
  pathway_name = control_data$name,
  database = control_data$database,
  ground_truth = ground_truth,
  qwen_default_louvain = unname(qwen_membership),
  original_default_louvain = unname(original_membership),
  stringsAsFactors = FALSE
)

similarity_comparison <- qwen_embedding_sim_df |>
  dplyr::rename(qwen_similarity = sim) |>
  dplyr::mutate(
    original_similarity = original_embedding_sim_matrix[
      cbind(match(from, expected_ids), match(to, expected_ids))
    ]
  )

embedding_metadata <- data.frame(
  field = c(
    "qwen_embedding_model", "original_embedding_model", "qwen_dimensions",
    "original_dimensions", "pathway_count", "similarity_measure",
    "clustering_method", "similarity_cutoff", "seed", "qwen_embedding_source"
  ),
  value = c(
    qwen_embedding_model, annotation_embedding_model,
    paste(dim(qwen_embedding_matrix), collapse = " x "),
    paste(dim(original_embedding_matrix), collapse = " x "),
    nrow(qwen_embedding_matrix), "cosine", mapa_clustering_method,
    mapa_similarity_cutoff, analysis_seed, qwen_embedding_source
  ),
  stringsAsFactors = FALSE
)

save(qwen_embedding_matrix, file = file.path(output_dir, "qwen_embedding_matrix.rda"))
save(
  qwen_embedding_sim_matrix,
  file = file.path(output_dir, "qwen_embedding_sim_matrix.rda")
)
save(qwen_embedding_sim_df, file = file.path(output_dir, "qwen_embedding_sim_df.rda"))
save(
  default_louvain_clustering_result,
  file = file.path(output_dir, "default_louvain_clustering_result.rda")
)
save(
  embedding_consistency_metrics,
  file = file.path(output_dir, "embedding_consistency_metrics.rda")
)
utils::write.csv(
  default_louvain_clustering_result,
  file.path(output_dir, "default_louvain_clustering_result.csv"),
  row.names = FALSE
)
utils::write.csv(
  embedding_consistency_metrics,
  file.path(output_dir, "embedding_consistency_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  similarity_comparison,
  file.path(output_dir, "pairwise_similarity_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  embedding_metadata,
  file.path(output_dir, "embedding_metadata.csv"),
  row.names = FALSE
)

similarity_plot <- ggplot2::ggplot(
  similarity_comparison,
  ggplot2::aes(x = original_similarity, y = qwen_similarity)
) +
  ggplot2::geom_point(size = 1.7, alpha = 0.55, colour = "#4C78A8") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                       colour = "#D55E00", linewidth = 0.8) +
  ggplot2::annotate(
    "text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.25,
    label = sprintf(
      "Pearson r = %.3f\nSpearman rho = %.3f",
      unname(pearson_test$estimate), unname(spearman_test$estimate)
    ),
    size = 4
  ) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = "Original embedding cosine similarity",
    y = "Qwen3-Embedding-8B cosine similarity"
  )
ggplot2::ggsave(
  file.path(output_dir, "embedding_similarity_consistency.pdf"),
  similarity_plot, width = 7, height = 6
)
ggplot2::ggsave(
  file.path(output_dir, "embedding_similarity_consistency.png"),
  similarity_plot, width = 7, height = 6, dpi = 300
)

ari_plot_data <- data.frame(
  embedding = factor(c("Original", "Qwen3-Embedding-8B"),
                     levels = c("Original", "Qwen3-Embedding-8B")),
  ARI = c(original_ari_ground_truth, qwen_ari_ground_truth)
)
ari_plot <- ggplot2::ggplot(
  ari_plot_data, ggplot2::aes(x = embedding, y = ARI, fill = embedding)
) +
  ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.3f", ARI)), vjust = -0.5, size = 4
  ) +
  ggplot2::scale_fill_manual(values = c("#9ECAE1", "#4C78A8")) +
  ggplot2::scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = NULL, y = "Adjusted Rand Index")
ggplot2::ggsave(
  file.path(output_dir, "clustering_ari_comparison.pdf"),
  ari_plot, width = 7, height = 6
)

umap_input <- 1 - qwen_embedding_sim_matrix
set.seed(analysis_seed)
umap_coordinates <- uwot::umap(
  umap_input,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean",
  verbose = FALSE
)
umap_data <- data.frame(
  pathway_id = expected_ids,
  UMAP1 = umap_coordinates[, 1L],
  UMAP2 = umap_coordinates[, 2L]
) |>
  dplyr::left_join(
    control_data |>
      dplyr::transmute(
        pathway_id = as.character(id), expected_module, database
      ),
    by = "pathway_id"
  ) |>
  dplyr::mutate(
    expected_module = factor(
      expected_module,
      levels = stringr::str_sort(unique(expected_module), numeric = TRUE)
    ),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )
ellipse_modules <- umap_data |>
  dplyr::count(expected_module) |>
  dplyr::filter(n > 3L) |>
  dplyr::pull(expected_module)
umap_plot <- ggplot2::ggplot(
  umap_data, ggplot2::aes(x = UMAP1, y = UMAP2)
) +
  ggplot2::geom_point(
    ggplot2::aes(colour = expected_module, shape = database),
    size = 5, alpha = 0.8
  ) +
  ggplot2::stat_ellipse(
    data = dplyr::filter(umap_data, expected_module %in% ellipse_modules),
    ggplot2::aes(colour = expected_module),
    type = "t", geom = "polygon", fill = NA, show.legend = FALSE
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "right") +
  ggplot2::labs(x = "UMAP1", y = "UMAP2", colour = "Expected module")
ggplot2::ggsave(
  file.path(output_dir, "qwen_embedding_umap.pdf"),
  umap_plot, width = 9, height = 7
)

summary_lines <- c(
  sprintf("Qwen embedding model: %s", qwen_embedding_model),
  sprintf("MAPA default clustering: %s, similarity cutoff %.2f", mapa_clustering_method, mapa_similarity_cutoff),
  sprintf("Pathways: %d", length(expected_ids)),
  sprintf("Qwen ARI versus ground truth: %.6f", qwen_ari_ground_truth),
  sprintf("Original ARI versus ground truth: %.6f", original_ari_ground_truth),
  sprintf("Qwen versus original clustering ARI: %.6f", cluster_ari_qwen_vs_original),
  sprintf("Pairwise cosine-similarity Pearson r: %.6f", unname(pearson_test$estimate)),
  sprintf("Pairwise cosine-similarity Spearman rho: %.6f", unname(spearman_test$estimate)),
  sprintf("Qwen clusters/edges: %d/%d", length(unique(qwen_membership)), nrow(qwen_cluster$edge_data)),
  sprintf("Original clusters/edges: %d/%d", length(unique(original_membership)), nrow(original_cluster$edge_data)),
  sprintf("Qwen embedding source: %s", qwen_embedding_source)
)
writeLines(summary_lines, file.path(output_dir, "embedding_clustering_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
