## Shared functions for the revised clustering benchmark.
##
## The three tools are evaluated on the same complete set of pathway IDs.
## If a method leaves a pathway unassigned, that pathway is represented as a
## singleton cluster for ARI calculation. Native coverage and the number of
## unassigned pathways are retained in the result tables.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mclust)
  library(tibble)
})

find_manuscript_root <- function() {
  configured_root <- getOption("mapa_manuscript_root", NULL)
  candidates <- unique(normalizePath(
    c(
      configured_root,
      getwd(),
      file.path(getwd(), "mapa_manuscript"),
      dirname(getwd()),
      file.path(dirname(getwd()), "mapa_manuscript")
    ),
    winslash = "/",
    mustWork = FALSE
  ))
  hit <- candidates[file.exists(file.path(candidates, "1_code", "100_tools.R"))]
  if (length(hit) == 0) {
    stop("Could not locate the mapa_manuscript project root.")
  }
  hit[[1]]
}

load_single_object <- function(path, object_name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Object '", object_name, "' was not found in ", path)
  }
  get(object_name, envir = env, inherits = FALSE)
}

relabel_complete_partition <- function(all_ids, assigned_ids, assigned_labels) {
  if (anyDuplicated(all_ids)) stop("all_ids must be unique")
  if (length(assigned_ids) != length(assigned_labels)) {
    stop("assigned_ids and assigned_labels must have equal lengths")
  }

  keep <- !is.na(assigned_ids) & !is.na(assigned_labels)
  assigned_ids <- as.character(assigned_ids[keep])
  assigned_labels <- as.character(assigned_labels[keep])
  if (anyDuplicated(assigned_ids)) {
    stop("A method returned duplicate assignments for at least one pathway ID")
  }

  out <- rep(NA_integer_, length(all_ids))
  names(out) <- all_ids
  if (length(assigned_ids) > 0) {
    numeric_labels <- match(assigned_labels, unique(assigned_labels))
    idx <- match(assigned_ids, all_ids)
    valid <- !is.na(idx)
    out[idx[valid]] <- numeric_labels[valid]
  }

  missing <- which(is.na(out))
  if (length(missing) > 0) {
    start <- if (all(is.na(out))) 0L else max(out, na.rm = TRUE)
    out[missing] <- seq.int(start + 1L, start + length(missing))
  }
  unname(out)
}

score_partition <- function(labels, ground_truth) {
  if (length(labels) != length(ground_truth) || anyNA(labels) || anyNA(ground_truth)) {
    stop("ARI requires complete, equally sized label vectors")
  }
  list(
    ari = mclust::adjustedRandIndex(labels, ground_truth),
    n_clusters = dplyr::n_distinct(labels)
  )
}

run_paver_grid <- function(input_dir, ground_truth_dt) {
  suppressPackageStartupMessages({
    library(PAVER)
    library(dynamicTreeCut)
  })

  env <- new.env(parent = globalenv())
  load(file.path(input_dir, "PAVER_01_input.rda"), envir = env)
  load(file.path(input_dir, "PAVER_02_embedding_matrix.rda"), envir = env)
  load(file.path(input_dir, "PAVER_03_term2name.rda"), envir = env)

  prepared <- PAVER::prepare_data(env$input, env$embedding_matrix, env$term2name)
  embedding_df <- prepared$embedding_mat
  embedding_mat <- embedding_df |>
    tibble::column_to_rownames("UniqueID") |>
    as.matrix()
  paver_ids <- prepared$prepared_data$GOID[
    match(rownames(embedding_mat), prepared$prepared_data$UniqueID)
  ]
  if (anyNA(paver_ids) || anyDuplicated(paver_ids)) {
    stop("Unable to construct a unique PAVER UniqueID-to-pathway mapping")
  }

  all_ids <- ground_truth_dt$id
  gt <- ground_truth_dt$ground_truth_label
  dist_obj <- PAVER:::cosine_dissimilarity(embedding_mat)
  dist_mat <- as.matrix(dist_obj)

  hclust_methods <- c("ward.D2", "average", "complete")
  min_cluster_sizes <- c(1L, 2L, 5L, 10L, 15L, 20L)
  max_core_values <- seq(0.01, 1, by = 0.01)
  trees <- setNames(lapply(hclust_methods, function(method) {
    stats::hclust(dist_obj, method = method)
  }), hclust_methods)

  evaluate_cut <- function(hclust_method, min_cluster_size,
                           max_core_scatter = NULL, min_gap = NULL) {
    args <- list(
      dendro = trees[[hclust_method]],
      distM = dist_mat,
      minClusterSize = min_cluster_size,
      verbose = 0
    )
    if (!is.null(max_core_scatter)) args$maxCoreScatter <- max_core_scatter
    if (!is.null(min_gap)) args$minGap <- min_gap

    raw_labels <- do.call(dynamicTreeCut::cutreeDynamic, args)
    if (length(raw_labels) != length(paver_ids)) {
      stop("cutreeDynamic returned an unexpected number of labels")
    }
    native_assigned <- raw_labels != 0 & !is.na(raw_labels)
    complete_labels <- relabel_complete_partition(
      all_ids = all_ids,
      assigned_ids = paver_ids[native_assigned],
      assigned_labels = raw_labels[native_assigned]
    )
    scored <- score_partition(complete_labels, gt)

    list(
      labels = complete_labels,
      ari = scored$ari,
      n_clusters = scored$n_clusters,
      native_n_clusters = dplyr::n_distinct(raw_labels[native_assigned]),
      native_n_assigned = sum(native_assigned),
      n_unassigned = sum(!native_assigned)
    )
  }

  result_list <- vector("list", length(hclust_methods) *
                          length(min_cluster_sizes) * length(max_core_values))
  counter <- 0L
  for (h_method in hclust_methods) {
    message("PAVER: hclust_method = ", h_method)
    for (min_size in min_cluster_sizes) {
      for (max_core in max_core_values) {
        counter <- counter + 1L
        min_gap <- (1 - max_core) * 3 / 4
        evaluated <- tryCatch(
          evaluate_cut(h_method, min_size, max_core, min_gap),
          error = function(e) e
        )
        if (inherits(evaluated, "error")) {
          result_list[[counter]] <- data.frame(
            hclust_method = h_method,
            minClusterSize = min_size,
            maxCoreScatter = max_core,
            minGap = min_gap,
            ARI = NA_real_,
            n_clusters = NA_integer_,
            native_n_clusters = NA_integer_,
            native_n_assigned = NA_integer_,
            n_unassigned = NA_integer_,
            status = conditionMessage(evaluated)
          )
        } else {
          result_list[[counter]] <- data.frame(
            hclust_method = h_method,
            minClusterSize = min_size,
            maxCoreScatter = max_core,
            minGap = min_gap,
            ARI = evaluated$ari,
            n_clusters = evaluated$n_clusters,
            native_n_clusters = evaluated$native_n_clusters,
            native_n_assigned = evaluated$native_n_assigned,
            n_unassigned = evaluated$n_unassigned,
            status = "ok"
          )
        }
      }
    }
  }
  results <- dplyr::bind_rows(result_list)

  best_row <- results |>
    dplyr::filter(status == "ok", !is.na(ARI)) |>
    dplyr::arrange(dplyr::desc(ARI), hclust_method, minClusterSize,
                   maxCoreScatter) |>
    dplyr::slice(1)
  best_eval <- evaluate_cut(
    best_row$hclust_method,
    best_row$minClusterSize,
    best_row$maxCoreScatter,
    best_row$minGap
  )

  default_eval <- evaluate_cut("ward.D2", 20L)
  public_default <- tryCatch({
    set.seed(234)
    PAVER::generate_themes(prepared)
    "ok"
  }, error = function(e) paste0("public generate_themes error: ", conditionMessage(e)))

  list(
    results = results,
    best = best_row,
    best_labels = best_eval$labels,
    default = data.frame(
      hclust_method = "ward.D2",
      minClusterSize = 20L,
      maxCoreScatter = NA_real_,
      minGap = NA_real_,
      ARI = default_eval$ari,
      n_clusters = default_eval$n_clusters,
      native_n_clusters = default_eval$native_n_clusters,
      native_n_assigned = default_eval$native_n_assigned,
      n_unassigned = default_eval$n_unassigned,
      status = public_default
    ),
    default_labels = default_eval$labels
  )
}

run_apear_grid <- function(enriched_result, ground_truth_dt) {
  suppressPackageStartupMessages(library(aPEAR))

  enrichment_df <- enriched_result@result
  if (!"geneID" %in% names(enrichment_df)) {
    stop("The enrichment result does not contain a geneID column")
  }
  if (anyDuplicated(enrichment_df$Description)) {
    stop("aPEAR requires unique pathway descriptions for unambiguous mapping")
  }

  all_ids <- ground_truth_dt$id
  gt <- ground_truth_dt$ground_truth_label
  id_by_description <- setNames(enrichment_df$ID, enrichment_df$Description)
  sim_methods <- c("jaccard", "cosine", "cor")
  cluster_methods <- c("markov", "hier")
  min_cluster_sizes <- c(2L, 3L, 5L, 10L)

  similarity <- setNames(lapply(sim_methods, function(method) {
    tryCatch(
      aPEAR::pathwaySimilarity(
        enrichment = enrichment_df,
        geneCol = "geneID",
        method = method
      ),
      error = function(e) e
    )
  }), sim_methods)

  evaluate_apear <- function(sim_method, cluster_method, min_cluster_size) {
    sim <- similarity[[sim_method]]
    if (inherits(sim, "error")) stop(conditionMessage(sim))
    set.seed(135)
    clusters <- aPEAR:::findClusters(
      sim = sim,
      minClusterSize = min_cluster_size,
      method = cluster_method,
      nameMethod = "pagerank",
      verbose = FALSE
    )
    assigned_ids <- unname(id_by_description[names(clusters)])
    complete_labels <- relabel_complete_partition(
      all_ids = all_ids,
      assigned_ids = assigned_ids,
      assigned_labels = clusters
    )
    scored <- score_partition(complete_labels, gt)
    list(
      labels = complete_labels,
      ari = scored$ari,
      n_clusters = scored$n_clusters,
      native_n_clusters = dplyr::n_distinct(clusters),
      native_n_assigned = length(clusters),
      n_unassigned = length(all_ids) - length(unique(assigned_ids))
    )
  }

  grid <- expand.grid(
    simMethod = sim_methods,
    clustMethod = cluster_methods,
    minClusterSize = min_cluster_sizes,
    stringsAsFactors = FALSE
  )
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    evaluated <- tryCatch(
      evaluate_apear(g$simMethod, g$clustMethod, g$minClusterSize),
      error = function(e) e
    )
    if (inherits(evaluated, "error")) {
      rows[[i]] <- cbind(g, data.frame(
        ARI = NA_real_, n_clusters = NA_integer_,
        native_n_clusters = NA_integer_, native_n_assigned = NA_integer_,
        n_unassigned = NA_integer_, status = conditionMessage(evaluated)
      ))
    } else {
      rows[[i]] <- cbind(g, data.frame(
        ARI = evaluated$ari, n_clusters = evaluated$n_clusters,
        native_n_clusters = evaluated$native_n_clusters,
        native_n_assigned = evaluated$native_n_assigned,
        n_unassigned = evaluated$n_unassigned, status = "ok"
      ))
    }
  }
  results <- dplyr::bind_rows(rows)
  best_row <- results |>
    dplyr::filter(status == "ok", !is.na(ARI)) |>
    dplyr::arrange(dplyr::desc(ARI), simMethod, clustMethod, minClusterSize) |>
    dplyr::slice(1)
  best_eval <- evaluate_apear(
    best_row$simMethod, best_row$clustMethod, best_row$minClusterSize
  )
  default_row <- results |>
    dplyr::filter(
      simMethod == "jaccard", clustMethod == "markov", minClusterSize == 2L
    ) |>
    dplyr::slice(1)
  default_eval <- evaluate_apear("jaccard", "markov", 2L)

  list(
    results = results,
    best = best_row,
    best_labels = best_eval$labels,
    default = default_row,
    default_labels = default_eval$labels
  )
}

run_enrichplot_grid <- function(enriched_result, ground_truth_dt) {
  suppressPackageStartupMessages(library(enrichplot))

  n_terms <- nrow(enriched_result@result)
  n_cluster_max <- if (n_terms < 100L) n_terms - 1L else 100L
  n_cluster_values <- seq.int(2L, n_cluster_max)
  layouts <- c("kk", "fr")
  min_edges <- c(0.1, 0.2, 0.3, 0.4)
  all_ids <- ground_truth_dt$id
  gt <- ground_truth_dt$ground_truth_label

  edo <- enrichplot::pairwise_termsim(
    enriched_result,
    method = "JC",
    showCategory = n_terms
  )
  if (anyDuplicated(edo@result$Description)) {
    stop("enrichplot requires unique pathway descriptions for unambiguous mapping")
  }
  id_by_description <- setNames(edo@result$ID, edo@result$Description)

  rows <- list()
  assignments <- list()
  counter <- 0L
  for (layout_name in layouts) {
    for (edge_cutoff in min_edges) {
      set.seed(123)
      base_plot <- tryCatch(
        enrichplot::emapplot(
          x = edo,
          showCategory = n_terms,
          layout = layout_name,
          min_edge = edge_cutoff,
          group = FALSE
        ),
        error = function(e) e
      )

      if (inherits(base_plot, "error")) {
        for (n_cluster in n_cluster_values) {
          counter <- counter + 1L
          rows[[counter]] <- data.frame(
            layout = layout_name, min_edge = edge_cutoff,
            nCluster = n_cluster, clusterFunction = "kmeans",
            ARI = NA_real_, n_clusters = NA_integer_,
            native_n_assigned = NA_integer_, n_unassigned = NA_integer_,
            status = conditionMessage(base_plot)
          )
        }
        next
      }

      plot_data <- base_plot$data
      coords <- data.frame(x = plot_data$x, y = plot_data$y)
      pathway_ids <- unname(id_by_description[plot_data$name])
      ## emapplot() computes the network layout first and then calls K-means.
      ## Preserve the RNG state immediately after layout generation so every
      ## candidate exactly matches a standalone set.seed(123); emapplot(...)
      ## call without recomputing the same layout for each nCluster.
      rng_after_layout <- get(".Random.seed", envir = .GlobalEnv)
      for (n_cluster in n_cluster_values) {
        counter <- counter + 1L
        assign(".Random.seed", rng_after_layout, envir = .GlobalEnv)
        clustered <- tryCatch(
          stats::kmeans(coords, centers = n_cluster),
          error = function(e) e
        )
        key <- paste(layout_name, edge_cutoff, n_cluster, sep = "|")
        if (inherits(clustered, "error")) {
          rows[[counter]] <- data.frame(
            layout = layout_name, min_edge = edge_cutoff,
            nCluster = n_cluster, clusterFunction = "kmeans",
            ARI = NA_real_, n_clusters = NA_integer_,
            native_n_assigned = nrow(plot_data),
            n_unassigned = length(all_ids) - length(unique(stats::na.omit(pathway_ids))),
            status = conditionMessage(clustered)
          )
        } else {
          complete_labels <- relabel_complete_partition(
            all_ids = all_ids,
            assigned_ids = pathway_ids,
            assigned_labels = clustered$cluster
          )
          scored <- score_partition(complete_labels, gt)
          rows[[counter]] <- data.frame(
            layout = layout_name, min_edge = edge_cutoff,
            nCluster = n_cluster, clusterFunction = "kmeans",
            ARI = scored$ari, n_clusters = scored$n_clusters,
            native_n_assigned = nrow(plot_data),
            n_unassigned = length(all_ids) - length(unique(stats::na.omit(pathway_ids))),
            status = "ok"
          )
          assignments[[key]] <- complete_labels
        }
      }
    }
  }
  results <- dplyr::bind_rows(rows)
  best_row <- results |>
    dplyr::filter(status == "ok", !is.na(ARI)) |>
    dplyr::arrange(dplyr::desc(ARI), layout, min_edge, nCluster) |>
    dplyr::slice(1)
  best_key <- paste(best_row$layout, best_row$min_edge,
                    best_row$nCluster, sep = "|")

  default_n_cluster <- floor(sqrt(n_terms))
  default_row <- results |>
    dplyr::filter(
      layout == "kk", min_edge == 0.2, nCluster == default_n_cluster
    ) |>
    dplyr::slice(1)
  default_key <- paste("kk", 0.2, default_n_cluster, sep = "|")

  list(
    results = results,
    best = best_row,
    best_labels = assignments[[best_key]],
    default = default_row,
    default_labels = assignments[[default_key]],
    n_terms = n_terms,
    n_cluster_max = n_cluster_max
  )
}

jaccard_index <- function(a, b) {
  union_size <- length(union(a, b))
  if (union_size == 0L) return(0)
  length(intersect(a, b)) / union_size
}

build_jaccard_heatmap_data <- function(assignments, tool_order) {
  all_ids <- assignments$ID
  gt <- assignments$ground_truth_label
  gt_levels <- sort(unique(gt))

  output <- list()
  counter <- 0L
  for (tool in tool_order) {
    labels <- assignments[[tool]]
    tool_levels <- sort(unique(labels))
    mat <- matrix(0, nrow = length(tool_levels), ncol = length(gt_levels),
                  dimnames = list(as.character(tool_levels), as.character(gt_levels)))
    for (i in seq_along(tool_levels)) {
      ids_tool <- all_ids[labels == tool_levels[[i]]]
      for (j in seq_along(gt_levels)) {
        ids_gt <- all_ids[gt == gt_levels[[j]]]
        mat[i, j] <- jaccard_index(ids_tool, ids_gt)
      }
    }
    best_gt <- apply(mat, 1, function(x) gt_levels[which.max(x)])
    ord <- order(best_gt, tool_levels)
    for (display_row in seq_along(ord)) {
      i <- ord[[display_row]]
      counter <- counter + 1L
      output[[counter]] <- data.frame(
        tool = tool,
        tool_module = tool_levels[[i]],
        display_row = display_row,
        best_ground_truth = best_gt[[i]],
        ground_truth_label = gt_levels,
        jaccard = as.numeric(mat[i, ]),
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(output)
}

plot_jaccard_heatmap <- function(hm_data, output_pdf, dataset_title,
                                 ari_by_tool = NULL) {
  suppressPackageStartupMessages({
    library(cowplot)
    library(patchwork)
  })

  tool_order <- c("MAPA", "PAVER", "aPEAR", "enrichplot")
  if (!is.null(ari_by_tool)) {
    available_ari <- ari_by_tool[tool_order]
    ranked_tools <- names(sort(
      available_ari[is.finite(available_ari)],
      decreasing = TRUE
    ))
    tool_order <- c(ranked_tools, setdiff(tool_order, ranked_tools))
  }
  tool_colors <- c(
    MAPA = "#DD5129", PAVER = "#43B284",
    aPEAR = "#0F7BA2", enrichplot = "#FAB255"
  )
  gt_levels <- sort(unique(hm_data$ground_truth_label))
  plot_list <- list()
  row_counts <- integer()

  for (tool in tool_order) {
    df <- hm_data |>
      dplyr::filter(.data$tool == .env$tool) |>
      dplyr::mutate(
        gt_factor = factor(ground_truth_label, levels = gt_levels),
        row_factor = factor(display_row, levels = rev(sort(unique(display_row))))
      )
    stopifnot(
      nrow(df) == sum(hm_data$tool == tool),
      dplyr::n_distinct(df$tool) == 1L
    )
    n_rows <- dplyr::n_distinct(df$display_row)
    row_counts <- c(row_counts, n_rows)
    show_y <- n_rows <= 80
    y_labels <- setNames(
      paste0("M", df$best_ground_truth[match(sort(unique(df$display_row)), df$display_row)]),
      as.character(sort(unique(df$display_row)))
    )
    panel_title <- tool
    if (!is.null(ari_by_tool) &&
        tool %in% names(ari_by_tool) &&
        is.finite(ari_by_tool[[tool]])) {
      panel_title <- sprintf("%s (ARI = %.3f)", tool, ari_by_tool[[tool]])
    }

    plot_list[[tool]] <- ggplot(df, aes(gt_factor, row_factor, fill = jaccard)) +
      geom_tile(color = "white", linewidth = if (n_rows <= 80) 0.25 else 0) +
      scale_fill_gradient(low = "white", high = tool_colors[[tool]],
                          limits = c(0, 1), guide = "none") +
      scale_y_discrete(labels = if (show_y) y_labels else NULL) +
      labs(title = panel_title, x = NULL, y = NULL) +
      theme_bw(base_size = 9) +
      theme(
        plot.title = element_text(color = tool_colors[[tool]], face = "bold",
                                  hjust = 0.5, size = 11),
        axis.text.y = if (show_y) element_text(size = 5) else element_blank(),
        axis.ticks.y = if (show_y) element_line(linewidth = 0.2) else element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(2, 4, 2, 4)
      )
  }

  bottom <- tool_order[[length(tool_order)]]
  plot_list[[bottom]] <- plot_list[[bottom]] +
    labs(x = "Ground-truth module") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 size = if (length(gt_levels) > 30) 5 else 7),
      axis.ticks.x = element_line(linewidth = 0.2)
    )

  stacked <- patchwork::wrap_plots(
    plot_list[tool_order], ncol = 1, heights = pmax(row_counts, 4)
  ) + patchwork::plot_annotation(title = dataset_title)

  legend_plot <- ggplot(data.frame(x = 1, y = 1, z = 0.5), aes(x, y, fill = z)) +
    geom_point(shape = 22, size = 4) +
    scale_fill_gradient(
      low = "white", high = "grey20", limits = c(0, 1),
      name = "Jaccard\nindex", breaks = seq(0, 1, 0.25)
    ) +
    theme_void() + theme(legend.position = "right")
  legend <- cowplot::get_legend(legend_plot)
  combined <- cowplot::plot_grid(stacked, legend, nrow = 1, rel_widths = c(0.94, 0.06))

  total_rows <- sum(row_counts)
  plot_height <- min(30, max(9, 4 + total_rows * 0.055))
  ggsave(output_pdf, combined, width = 11, height = plot_height,
         limitsize = FALSE)
}

run_complete_benchmark <- function(config) {
  root <- find_manuscript_root()
  input_dir <- file.path(root, config$input_dir)
  output_dir <- file.path(root, config$output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  enriched_result <- load_single_object(
    file.path(input_dir, "enriched_result.rda"), "enriched_result"
  )
  ground_truth_dt <- load_single_object(
    file.path(input_dir, "ground_truth_dt.rda"), "ground_truth_dt"
  ) |>
    dplyr::arrange(match(id, enriched_result@result$ID))
  if (nrow(ground_truth_dt) != nrow(enriched_result@result) ||
      anyNA(match(enriched_result@result$ID, ground_truth_dt$id))) {
    stop("Ground truth and enrichment result do not contain identical pathway IDs")
  }
  ground_truth_dt <- ground_truth_dt[match(enriched_result@result$ID,
                                           ground_truth_dt$id), ]

  message("Running PAVER grid for ", config$dataset_name)
  paver <- run_paver_grid(input_dir, ground_truth_dt)
  message("Running aPEAR grid for ", config$dataset_name)
  apear <- run_apear_grid(enriched_result, ground_truth_dt)
  message("Running enrichplot grid for ", config$dataset_name)
  enrichplot_res <- run_enrichplot_grid(enriched_result, ground_truth_dt)

  readr::write_csv(paver$results, file.path(output_dir, "paver_sensitivity.csv"))
  readr::write_csv(apear$results, file.path(output_dir, "apear_sensitivity.csv"))
  readr::write_csv(enrichplot_res$results,
                   file.path(output_dir, "enrichplot_sensitivity.csv"))
  save(paver, file = file.path(output_dir, "paver_sensitivity.rda"))
  save(apear, file = file.path(output_dir, "apear_sensitivity.rda"))
  save(enrichplot_res, file = file.path(output_dir, "enrichplot_sensitivity.rda"))

  mapa <- config$get_mapa_results(root, ground_truth_dt)

  paver_best <- paver$best
  apear_best <- apear$best
  enrich_best <- enrichplot_res$best
  summary <- dplyr::bind_rows(
    mapa$summary,
    data.frame(
      method = "PAVER",
      default_ari = paver$default$ARI,
      default_n_clusters = paver$default$n_clusters,
      default_native_n_clusters = paver$default$native_n_clusters,
      default_n_unassigned = paver$default$n_unassigned,
      default_settings = paste0(
        "hclust_method=ward.D2; cutreeDynamic defaults; minClusterSize=20",
        "; unassigned_as_singletons"
      ),
      optimized_ari = paver_best$ARI,
      optimized_n_clusters = paver_best$n_clusters,
      optimized_native_n_clusters = paver_best$native_n_clusters,
      optimized_n_unassigned = paver_best$n_unassigned,
      optimized_settings = sprintf(
        "hclust_method=%s; minClusterSize=%d; maxCoreScatter=%.2f; minGap=%.4f",
        paver_best$hclust_method, paver_best$minClusterSize,
        paver_best$maxCoreScatter, paver_best$minGap
      ),
      default_status = paver$default$status
    ),
    data.frame(
      method = "aPEAR",
      default_ari = apear$default$ARI,
      default_n_clusters = apear$default$n_clusters,
      default_native_n_clusters = apear$default$native_n_clusters,
      default_n_unassigned = apear$default$n_unassigned,
      default_settings = "simMethod=jaccard; clustMethod=markov; minClusterSize=2; excluded_as_singletons",
      optimized_ari = apear_best$ARI,
      optimized_n_clusters = apear_best$n_clusters,
      optimized_native_n_clusters = apear_best$native_n_clusters,
      optimized_n_unassigned = apear_best$n_unassigned,
      optimized_settings = sprintf(
        "simMethod=%s; clustMethod=%s; minClusterSize=%d; excluded_as_singletons",
        apear_best$simMethod, apear_best$clustMethod,
        apear_best$minClusterSize
      ),
      default_status = apear$default$status
    ),
    data.frame(
      method = "enrichplot",
      default_ari = enrichplot_res$default$ARI,
      default_n_clusters = enrichplot_res$default$n_clusters,
      default_native_n_clusters = enrichplot_res$default$n_clusters,
      default_n_unassigned = enrichplot_res$default$n_unassigned,
      default_settings = sprintf(
        "layout=kk; nCluster=floor(sqrt(n))=%d; min_edge=0.2; clusterFunction=kmeans",
        floor(sqrt(enrichplot_res$n_terms))
      ),
      optimized_ari = enrich_best$ARI,
      optimized_n_clusters = enrich_best$n_clusters,
      optimized_native_n_clusters = enrich_best$n_clusters,
      optimized_n_unassigned = enrich_best$n_unassigned,
      optimized_settings = sprintf(
        "layout=%s; nCluster=%d; min_edge=%.1f; clusterFunction=kmeans",
        enrich_best$layout, enrich_best$nCluster, enrich_best$min_edge
      ),
      default_status = enrichplot_res$default$status
    )
  )
  readr::write_csv(summary, file.path(output_dir, "default_optimized_ari_settings.csv"))
  benchmark_summary <- summary
  save(benchmark_summary,
       file = file.path(output_dir, "default_optimized_ari_settings.rda"))

  assignments <- data.frame(
    ID = ground_truth_dt$id,
    ground_truth_label = ground_truth_dt$ground_truth_label,
    MAPA = mapa$optimized_labels,
    PAVER = paver$best_labels,
    aPEAR = apear$best_labels,
    enrichplot = enrichplot_res$best_labels,
    check.names = FALSE
  )
  stopifnot(!anyNA(assignments), nrow(assignments) == nrow(ground_truth_dt))
  readr::write_csv(assignments,
                   file.path(output_dir, "optimized_clustering_assignments.csv"))
  optimized_clustering_assignments <- assignments
  save(optimized_clustering_assignments,
       file = file.path(output_dir, "optimized_clustering_assignments.rda"))

  default_assignments <- data.frame(
    ID = ground_truth_dt$id,
    ground_truth_label = ground_truth_dt$ground_truth_label,
    MAPA = mapa$default_labels,
    PAVER = paver$default_labels,
    aPEAR = apear$default_labels,
    enrichplot = enrichplot_res$default_labels,
    check.names = FALSE
  )
  stopifnot(!anyNA(default_assignments),
            nrow(default_assignments) == nrow(ground_truth_dt))
  readr::write_csv(default_assignments,
                   file.path(output_dir, "default_clustering_assignments.csv"))
  default_clustering_assignments <- default_assignments
  save(default_clustering_assignments,
       file = file.path(output_dir, "default_clustering_assignments.rda"))

  hm_data <- build_jaccard_heatmap_data(
    assignments, c("MAPA", "PAVER", "aPEAR", "enrichplot")
  )
  readr::write_csv(hm_data,
                   file.path(output_dir, "optimized_jaccard_heatmap_data.csv"))
  plot_jaccard_heatmap(
    hm_data,
    file.path(output_dir, "optimized_heatmap_consistency.pdf"),
    paste0(config$dataset_name, ": optimized clustering consistency"),
    stats::setNames(summary$optimized_ari, summary$method)
  )

  writeLines(capture.output(sessionInfo()),
             file.path(output_dir, "sessionInfo.txt"))
  message("Completed benchmark: ", output_dir)
  print(summary)
  invisible(list(summary = summary, assignments = assignments))
}
