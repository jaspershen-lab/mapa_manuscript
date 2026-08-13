## Benchmarking on the informative GO set benchmark (k = 20)
## Step 2: Run enrichplot, aPEAR and PAVER, scan their tunable parameters and
## record the ARI against the ground truth modules.
## Adapted from 1_code/2_control_data/05_benchmarking/02_run_analysis.R
## (nCluster ranges widened for 313 terms / 57 ground-truth modules)

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(patchwork)
library(ggplot2)

setwd("3_data_analysis/11_informative_go_benchmarking/")
load("enriched_result.rda")
load("ground_truth_dt.rda")

n_terms <- nrow(enriched_result@result)

# enrichplot ====
library(enrichplot)
edo <- enrichplot::pairwise_termsim(enriched_result, showCategory = n_terms)

evaluate_clustering_ari <- function(edo,
                                    ground_truth_dt,
                                    nCluster,
                                    showCategory = n_terms,
                                    layout = "fr",
                                    min_edge = 0.2,
                                    color_edge = "grey",
                                    clusterFunction = stats::kmeans,
                                    verbose = TRUE) {

  ari_results <- data.frame(
    nCluster = integer(),
    ARI = numeric()
  )

  for (n_clust in nCluster) {

    if (verbose) {
      cat("Testing nCluster =", n_clust, "...\n")
    }

    tryCatch({
      set.seed(123)
      p <- enrichplot::emapplot(
        x = edo,
        showCategory = showCategory,
        layout = layout,
        min_edge = min_edge,
        color_edge = color_edge,
        node_label = "group",
        group = TRUE,
        clusterFunction = clusterFunction,
        nCluster = n_clust
      )

      clustering_data <- p$data
      cluster_results <- clustering_data[, c("name", "color2")]
      colnames(cluster_results) <- c("Description", "cluster_label")

      cluster_results <- cluster_results |>
        dplyr::left_join(edo@result[, c("ID", "Description")],
                         by = "Description")

      unique_clusters <- unique(cluster_results$cluster_label)
      cluster_mapping <- data.frame(
        cluster_label = unique_clusters,
        cluster_id = 1:length(unique_clusters)
      )

      cluster_results <- merge(cluster_results, cluster_mapping, by = "cluster_label")

      enrichplot_label <- cluster_results |>
        dplyr::select(ID, cluster_id) |>
        dplyr::rename(id = ID, enrichplot_label = cluster_id)

      combined_cluster_res <- ground_truth_dt |>
        dplyr::left_join(enrichplot_label, by = "id")

      combined_cluster_res_clean <- combined_cluster_res |>
        dplyr::filter(!is.na(enrichplot_label) & !is.na(ground_truth_label))

      if (nrow(combined_cluster_res_clean) > 0) {
        ari_value <- mclust::adjustedRandIndex(
          combined_cluster_res_clean$enrichplot_label,
          combined_cluster_res_clean$ground_truth_label
        )

        ari_results <- rbind(ari_results, data.frame(
          nCluster = n_clust,
          ARI = ari_value
        ))

        if (verbose) {
          cat("  ARI =", round(ari_value, 4),
              "| Terms clustered:", nrow(combined_cluster_res_clean),
              "| Clusters found:", length(unique(cluster_results$cluster_id)), "\n")
        }
      } else {
        if (verbose) {
          cat("  No matching terms found for ARI calculation\n")
        }
      }

    }, error = function(e) {
      if (verbose) {
        cat("  Error with nCluster =", n_clust, ":", e$message, "\n")
      }
    })
  }

  return(ari_results)
}

enrichplot_ari <- evaluate_clustering_ari(edo = edo,
                                          ground_truth_dt = ground_truth_dt,
                                          nCluster = 2:100)

print(enrichplot_ari[which.max(enrichplot_ari$ARI), ])
save(enrichplot_ari, file = "enrichplot_ari.rda")

# aPEAR ====
library(aPEAR)

calculate_apear_ari <- function(enriched_result,
                                ground_truth_dt,
                                simMethod = "jaccard",
                                innerCutoff = 0.1,
                                outerCutoff = 0.5,
                                clustNameMethod = "pagerank",
                                colorBy = "pvalue",
                                colorType = "pval",
                                nodeSize = "Count",
                                drawEllipses = TRUE,
                                minClusterSize = 2,
                                fontSize = 3) {

  set.seed(135)

  cluster_methods <- c("markov", "hier") # Bug in "spectral" (https://github.com/kerseviciute/aPEAR/issues/9)

  ari_results <- list()

  for (method in cluster_methods) {

    cat("Processing clustering method:", method, "\n")

    tryCatch({
      apear_res <- enrichmentNetwork(
        enrichment = enriched_result@result,
        simMethod = simMethod,
        innerCutoff = innerCutoff,
        outerCutoff = outerCutoff,
        clustMethod = method,
        clustNameMethod = clustNameMethod,
        colorBy = colorBy,
        colorType = colorType,
        nodeSize = nodeSize,
        plotOnly = FALSE,
        drawEllipses = drawEllipses,
        minClusterSize = minClusterSize,
        fontSize = fontSize
      )

      apear_cluster <- apear_res$clusters |>
        dplyr::rename(Description = ID)

      pathway_id_desc <- enriched_result@result |>
        dplyr::select(ID, Description) |>
        dplyr::left_join(apear_cluster, by = "Description")

      cluster_mapper <- data.frame(
        Cluster = unique(apear_cluster$Cluster),
        cluster_label = 1:length(unique(apear_cluster$Cluster))
      )

      all_cluster_label <- pathway_id_desc |>
        dplyr::left_join(cluster_mapper, by = "Cluster") |>
        dplyr::arrange(cluster_label)

      ## Assign sequential labels to singletons (NA values)
      na_indices <- which(is.na(all_cluster_label$cluster_label))
      if (length(na_indices) > 0) {
        max_label <- max(all_cluster_label$cluster_label, na.rm = TRUE)
        singleton_labels <- seq(max_label + 1, max_label + length(na_indices))
        all_cluster_label$cluster_label[na_indices] <- singleton_labels
      }

      combined_cluster_res <- all_cluster_label |>
        dplyr::select(ID, cluster_label) |>
        dplyr::rename(id = ID) |>
        dplyr::right_join(ground_truth_dt, by = "id")

      ari_value <- mclust::adjustedRandIndex(
        combined_cluster_res$cluster_label,
        combined_cluster_res$ground_truth_label
      )

      n_clusters <- length(unique(all_cluster_label$cluster_label))

      ari_results[[method]] <- list(
        ari = ari_value,
        n_clusters = n_clusters
      )

      cat("ARI for", method, ":", round(ari_value, 4), "| Clusters:", n_clusters, "\n")

    }, error = function(e) {
      cat("Error in method", method, ":", e$message, "\n")
      ari_results[[method]] <<- list(ari = NA, n_clusters = NA)
    })
  }

  ari_summary <- data.frame(
    Method = names(ari_results),
    ARI = sapply(ari_results, function(x) x$ari),
    N_Clusters = sapply(ari_results, function(x) x$n_clusters),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(desc(ARI))

  return(ari_summary)
}

apear_ari <- calculate_apear_ari(enriched_result = enriched_result,
                                 ground_truth_dt = ground_truth_dt)
print(apear_ari)
save(apear_ari, file = "apear_ari.rda")

# PAVER ====
library(PAVER)
load("PAVER_01_input.rda")
load("PAVER_02_embedding_matrix.rda")
load("PAVER_03_term2name.rda")

PAVER_result <- prepare_data(input, embedding_matrix, term2name)

find_optimal_paver_clustering <- function(PAVER_result,
                                          ground_truth_dt,
                                          step_size = 0.01,
                                          minClusterSize = 1,
                                          verbose = TRUE,
                                          ...) {

  set.seed(234)

  maxCoreScatter_values <- seq(0.01, 1, by = step_size)

  results <- data.frame(
    maxCoreScatter = numeric(),
    nClust = integer(),
    ari = numeric(),
    stringsAsFactors = FALSE
  )

  cat("Testing", length(maxCoreScatter_values), "maxCoreScatter values...\n")

  for (i in seq_along(maxCoreScatter_values)) {
    maxCoreScatter <- maxCoreScatter_values[i]
    minGap <- (1 - maxCoreScatter) * 3 / 4

    tryCatch({
      current_result <- generate_themes(
        PAVER_result,
        maxCoreScatter = maxCoreScatter,
        minGap = minGap,
        minClusterSize = minClusterSize,
        ...
      )

      clustered_input <- PAVER_export(current_result)

      unique_clusters <- unique(clustered_input$Cluster)
      nClust <- length(unique_clusters)

      if (nClust < 2) {
        cat("Skipping maxCoreScatter =", maxCoreScatter, "(insufficient clusters)\n")
        next
      }

      cluster_mapping <- data.frame(
        Cluster = unique_clusters,
        paver_cluster_label = seq_len(nClust)
      )

      clustered_input <- clustered_input |>
        dplyr::left_join(cluster_mapping, by = "Cluster")

      combined_cluster_res <- clustered_input |>
        dplyr::select(GOID, paver_cluster_label) |>
        dplyr::rename(id = GOID) |>
        dplyr::left_join(ground_truth_dt, by = "id")

      combined_cluster_res <- combined_cluster_res[!is.na(combined_cluster_res$ground_truth_label), ]

      if (nrow(combined_cluster_res) > 0) {
        ari_score <- mclust::adjustedRandIndex(
          combined_cluster_res$paver_cluster_label,
          combined_cluster_res$ground_truth_label
        )

        results <- rbind(results, data.frame(
          maxCoreScatter = maxCoreScatter,
          nClust = nClust,
          ari = ari_score
        ))

        if (i %% 10 == 0) {
          cat("Completed", i, "/", length(maxCoreScatter_values),
              "- Current: nClust =", nClust, ", ARI =", round(ari_score, 4), "\n")
        }
      }

    }, error = function(e) {
      if (verbose) {
        cat("Skipping maxCoreScatter =", maxCoreScatter, "(error:", e$message, ")\n")
      }
    })
  }

  optimal_results <- results |>
    dplyr::group_by(nClust) |>
    dplyr::slice_max(ari, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(nClust) |>
    dplyr::select(nClust, maxCoreScatter, ari)

  cat("\nOptimization complete! Found", nrow(optimal_results), "unique cluster configurations.\n")
  cat("Best overall ARI:", max(optimal_results$ari),
      "with", optimal_results$nClust[which.max(optimal_results$ari)], "clusters\n")

  return(optimal_results)
}

paver_clustering_res <- find_optimal_paver_clustering(
  PAVER_result = PAVER_result,
  ground_truth_dt = ground_truth_dt,
  step_size = 0.01,
  minClusterSize = 1
)

save(paver_clustering_res, file = "paver_ari.rda")
