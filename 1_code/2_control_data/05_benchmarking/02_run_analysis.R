library(r4projects)
setwd(get_project_wd())
rm(list = ls())
# source('1_code/100_tools.R')
library(patchwork)
library(ggplot2)

# setwd(get_project_wd())
# control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)
# ground_truth_dt <-
#   control_dt |>
#   dplyr::mutate(ground_truth_label = as.numeric(sub(pattern = "Functional_module_", replacement = "", x = expected_module)))
# ground_truth_dt <- ground_truth_dt |> dplyr::select(id, ground_truth_label)
# save(ground_truth_dt, file = "ground_truth_dt.rda")

setwd("3_data_analysis/02_control_data/05_benchmarking/")
load("enriched_result.rda")
load("ground_truth_dt.rda")

# enrichplot ====
# https://bioconductor.org/packages/devel/bioc/manuals/enrichplot/man/enrichplot.pdf
# library(clusterProfiler)
# BiocManager::install("enrichplot")
library(enrichplot)
edo <- enrichplot::pairwise_termsim(enriched_result)
evaluate_clustering_ari <- function(edo,
                                    ground_truth_dt,
                                    nCluster,
                                    showCategory = 44,
                                    layout = "fr",
                                    min_edge = 0.2,
                                    color_edge = "grey",
                                    clusterFunction = stats::kmeans,
                                    verbose = TRUE) {

  # Initialize results dataframe
  ari_results <- data.frame(
    nCluster = integer(),
    ARI = numeric(),
    n_terms_clustered = integer(),
    n_clusters_found = integer()
  )

  # Loop through different nCluster values
  for (n_clust in nCluster) {

    if (verbose) {
      cat("Testing nCluster =", n_clust, "...\n")
    }

    tryCatch({
      set.seed(123)
      # Generate emapplot with current nCluster
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

      # Extract clustering results
      clustering_data <- p$data
      cluster_results <- clustering_data[, c("name", "color2")]
      colnames(cluster_results) <- c("Description", "cluster_label")

      # Join with enriched result to get IDs
      cluster_results <- cluster_results |>
        dplyr::left_join(edo@result[, c("ID", "Description")],
                         by = "Description")

      # Create numerical cluster IDs
      unique_clusters <- unique(cluster_results$cluster_label)
      cluster_mapping <- data.frame(
        cluster_label = unique_clusters,
        cluster_id = 1:length(unique_clusters)
      )

      cluster_results <- merge(cluster_results, cluster_mapping, by = "cluster_label")

      # Prepare data for ARI calculation
      enrichplot_label <- cluster_results |>
        dplyr::select(ID, cluster_id) |>
        dplyr::rename(id = ID, enrichplot_label = cluster_id)

      combined_cluster_res <- ground_truth_dt |>
        dplyr::left_join(enrichplot_label, by = "id")

      # Remove rows with missing cluster assignments
      combined_cluster_res_clean <- combined_cluster_res |>
        dplyr::filter(!is.na(enrichplot_label) & !is.na(ground_truth_label))

      # Calculate ARI
      if (nrow(combined_cluster_res_clean) > 0) {
        ari_value <- mclust::adjustedRandIndex(
          combined_cluster_res_clean$enrichplot_label,
          combined_cluster_res_clean$ground_truth_label
        )

        # Store results
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
                                          nCluster = 2:43)

enrichplot_ari[which.max(enrichplot_ari$ARI),]
plot(enrichplot_ari$nCluster, enrichplot_ari$ARI)

p <- enrichplot::emapplot(
  x = edo,
  showCategory = 44,
  layout = "fr",
  min_edge = 0.2,
  color_edge = "grey",
  node_label = "group",
  group = TRUE,
  clusterFunction = stats::kmeans,
  nCluster = 17
)
p

ggsave(plot = p, filename = "enrichplot_best_clustering_network.pdf",
       width = 10, height = 8)
save(enrichplot_ari, file = "enrichplot_ari.rda")

# aPEAR ====
# library(devtools)
# install_github('ievaKer/aPEAR')

library(aPEAR)

# library(Spectrum)
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

  # Define clustering methods to test
  cluster_methods <- c("markov", "hier") # Bug in "spectral" (https://github.com/kerseviciute/aPEAR/issues/9)

  # Initialize results list
  ari_results <- list()
  cluster_results <- list()

  # Loop through each clustering method
  for (method in cluster_methods) {

    cat("Processing clustering method:", method, "\n")

    tryCatch({
      # Run enrichmentNetwork for current method
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

      # Extract cluster assignments
      apear_cluster <- apear_res$clusters |>
        dplyr::rename(Description = ID)

      # Create pathway ID to description mapping
      pathway_id_desc <- enriched_result@result |>
        dplyr::select(ID, Description) |>
        dplyr::left_join(apear_cluster, by = "Description")

      # Create cluster mapper (assign numeric labels to clusters)
      cluster_mapper <- data.frame(
        Cluster = unique(apear_cluster$Cluster),
        cluster_label = 1:length(unique(apear_cluster$Cluster))
      )

      # Join with cluster mapper
      all_cluster_label <- pathway_id_desc |>
        dplyr::left_join(cluster_mapper, by = "Cluster") |>
        dplyr::arrange(cluster_label)

      # Assign sequential labels to singletons (NA values)
      na_indices <- which(is.na(all_cluster_label$cluster_label))
      if (length(na_indices) > 0) {
        max_label <- max(all_cluster_label$cluster_label, na.rm = TRUE)
        singleton_labels <- seq(max_label + 1, max_label + length(na_indices))
        all_cluster_label$cluster_label[na_indices] <- singleton_labels
      }

      # Combine with ground truth
      combined_cluster_res <- all_cluster_label |>
        dplyr::select(ID, cluster_label) |>
        dplyr::rename(id = ID) |>
        dplyr::right_join(ground_truth_dt, by = "id")

      # Calculate ARI
      ari_value <- mclust::adjustedRandIndex(
        combined_cluster_res$cluster_label,
        combined_cluster_res$ground_truth_label
      )

      # Count number of clusters (including singletons)
      n_clusters <- length(unique(all_cluster_label$cluster_label))

      # Store results
      ari_results[[method]] <- list(
        ari = ari_value,
        n_clusters = n_clusters
      )
      cluster_results[[method]] <- list(
        clusters = all_cluster_label,
        combined = combined_cluster_res,
        plot = apear_res$plot
      )

      cat("ARI for", method, ":", round(ari_value, 4), "| Clusters:", n_clusters, "\n")

    }, error = function(e) {
      cat("Error in method", method, ":", e$message, "\n")
      ari_results[[method]] <- list(ari = NA, n_clusters = NA)
      cluster_results[[method]] <- NULL
    })
  }

  # Create summary results
  ari_summary <- data.frame(
    Method = names(ari_results),
    ARI = sapply(ari_results, function(x) x$ari),
    N_Clusters = sapply(ari_results, function(x) x$n_clusters),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(desc(ARI))

  # Return just the ARI summary dataframe
  return(ari_summary)
}

apear_ari <- calculate_apear_ari(enriched_result = enriched_result,
                                 ground_truth_dt = ground_truth_dt)
best_param <- apear_ari[which.max(apear_ari$ARI),]

apear_res_best <- enrichmentNetwork(
  enrichment = enriched_result@result,
  simMethod = "jaccard",
  innerCutoff = 0.1,
  outerCutoff = 0.5,
  clustMethod = best_param$Method,
  clustNameMethod = "pagerank",
  colorBy = "pvalue",
  colorType = "pval",
  nodeSize = "Count",
  plotOnly = FALSE,
  drawEllipses = TRUE,
  minClusterSize = 2,
  fontSize = 3
)

p <- apear_res_best$plot
apear_best_clustering <- apear_res_best$clusters

ggplot2::ggsave(plot = p,
                filename = "apear_best_clustering_network.pdf",
                width = 10, height = 8)
save(apear_ari, file = "apear_ari.rda")

# PAVER ====
# remotes::install_github("willgryan/PAVER", build_vignettes = TRUE, force = TRUE)
# browseVignettes("PAVER") # See the instruction for using PAVER
library(PAVER)
load("PAVER_01_input.rda")
load("PAVER_02_embedding_matrix.rda")
load("PAVER_03_term2name.rda")

# PAVER example data
# input <- gsea_example
# embeddings <- readRDS(url("https://github.com/willgryan/PAVER_embeddings/raw/main/2023-03-06/embeddings_2023-03-06.RDS"))
# term2name <- readRDS(url("https://github.com/willgryan/PAVER_embeddings/raw/main/2023-03-06/term2name_2023-03-06.RDS"))

PAVER_result <- prepare_data(input, embedding_matrix, term2name)

find_optimal_paver_clustering <- function(PAVER_result,
                                          ground_truth_dt,
                                          step_size = 0.01,
                                          minClusterSize = 1,
                                          verbose = TRUE,
                                          ...) {

  # Set seed for reproducibility
  set.seed(234)

  # Generate sequence of maxCoreScatter values
  maxCoreScatter_values <- seq(0.01, 1, by = step_size)  # Start from 0.01 to avoid edge cases

  # Initialize results storage
  results <- data.frame(
    maxCoreScatter = numeric(),
    nClust = integer(),
    ari = numeric(),
    stringsAsFactors = FALSE
  )

  # Progress tracking
  cat("Testing", length(maxCoreScatter_values), "maxCoreScatter values...\n")
  cat("Note: Some parameter combinations may fail due to clustering constraints - these will be skipped automatically.\n\n")

  for (i in seq_along(maxCoreScatter_values)) {
    maxCoreScatter <- maxCoreScatter_values[i]

    # Calculate minGap based on maxCoreScatter
    minGap <- (1 - maxCoreScatter) * 3 / 4

    tryCatch({
      if (verbose) {
        cat("Testing maxCoreScatter =", maxCoreScatter, ", minGap =", minGap, "\n")
      }

      # Generate themes with current parameters
      current_result <- generate_themes(
        PAVER_result,
        maxCoreScatter = maxCoreScatter,
        minGap = minGap,
        minClusterSize = minClusterSize,
        ...
      )

      if (verbose) cat("  - generate_themes completed\n")

      # Export clustered data
      clustered_input <- PAVER_export(current_result)

      if (verbose) cat("  - PAVER_export completed\n")

      # Check if clustering was successful (at least 2 clusters)
      unique_clusters <- unique(clustered_input$Cluster)
      nClust <- length(unique_clusters)

      if (nClust < 2) {
        cat("Skipping maxCoreScatter =", maxCoreScatter, "(insufficient clusters)\n")
        next
      }

      # Create cluster mapping
      cluster_mapping <- data.frame(
        Cluster = unique_clusters,
        paver_cluster_label = seq_len(nClust)
      )

      # Join with cluster mapping
      clustered_input <- clustered_input |>
        dplyr::left_join(cluster_mapping, by = "Cluster")

      # Prepare data for ARI calculation
      combined_cluster_res <- clustered_input |>
        dplyr::select(GOID, paver_cluster_label) |>
        dplyr::rename(id = GOID) |>
        dplyr::left_join(ground_truth_dt, by = "id")

      # Remove any missing ground truth labels
      combined_cluster_res <- combined_cluster_res[!is.na(combined_cluster_res$ground_truth_label), ]

      # Calculate ARI
      if (nrow(combined_cluster_res) > 0) {
        ari_score <- mclust::adjustedRandIndex(
          combined_cluster_res$paver_cluster_label,
          combined_cluster_res$ground_truth_label
        )

        # Store results
        results <- rbind(results, data.frame(
          maxCoreScatter = maxCoreScatter,
          nClust = nClust,
          ari = ari_score
        ))

        # Progress update
        if (i %% 10 == 0) {
          cat("Completed", i, "/", length(maxCoreScatter_values),
              "- Current: nClust =", nClust, ", ARI =", round(ari_score, 4), "\n")
        }
      }

    }, error = function(e) {
      # Silently skip errors and continue optimization
      if (verbose) {
        cat("Skipping maxCoreScatter =", maxCoreScatter, "(error:", e$message, ")\n")
      }
    })
  }

  # Filter to keep only the best ARI for each number of clusters
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

## Best params for PAVER
minClusterSize <- 1
maxCoreScatter <- max(paver_clustering_res$maxCoreScatter)
minGap <- (1 - maxCoreScatter) * 3 / 4
PAVER_result <- generate_themes(
  PAVER_result,
  maxCoreScatter = maxCoreScatter,
  minGap = minGap,
  minClusterSize = minClusterSize
)

## Result of PAVER
### 1. Result report
clustered_input <- PAVER_export(PAVER_result)

### 2. Theme plot
paver_theme_plot <- PAVER_theme_plot(PAVER_result)
paver_theme_plot

### 3. Interpretation plot
interpreted_plot <- PAVER_interpretation_plot(PAVER_result)
interpreted_plot

p <- interpreted_plot + paver_theme_plot
p

ggsave(filename = "paver_best_clustering_network.pdf",
       plot = p,
       width = 12, height = 6)


## MAPA ====
load("3_data_analysis/07_example_output/openai_semantic_sim_matrix.rda")
load("3_data_analysis/02_control_data/05_benchmarking/enriched_result.rda")

enriched_res_temp <- openai_semantic_sim_matrix$enriched_pathway

enriched_res_temp@enrichment_go_result@result <-
  enriched_result@result |>
  dplyr::filter(grepl("^GO:", ID)) |>
  dplyr::mutate(p_adjust = 0.01)

enriched_res_temp@enrichment_kegg_result@result <-
  enriched_result@result |>
  dplyr::filter(grepl("^hsa", ID)) |>
  dplyr::mutate(p_adjust = 0.01)

enriched_res_temp@enrichment_reactome_result@result <-
  enriched_result@result |>
  dplyr::filter(grepl("^R-HSA", ID)) |>
  dplyr::mutate(p_adjust = 0.01)

enriched_res_control_dt_1 <- enriched_res_temp
save(enriched_res_control_dt_1, file = "3_data_analysis/02_control_data/05_benchmarking/enriched_res_control_dt_1.rda")

bioembed_sim_res <-
  get_bioembedsim(
    object = enriched_res_control_dt_1,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg", "reactome"),
    count.cutoff.go = 0,
    count.cutoff.kegg = 0,
    count.cutoff.reactome = 0
  )

save(bioembed_sim_res, file = "3_data_analysis/02_control_data/05_benchmarking/mapa_bioembed_sim_res.rda")

bioembed_sim_res$enriched_pathway@enrichment_go_result <-
  bioembed_sim_res$enriched_pathway@enrichment_go_result |>
  dplyr::mutate(ONTOLOGY = "")

bioembed_sim_res$enriched_pathway@enrichment_kegg_result <-
  bioembed_sim_res$enriched_pathway@enrichment_kegg_result |>
  dplyr::mutate(category = "") |>
  dplyr::mutate(subcategory = "")

# fm_result <- get_functional_modules(
#   object = bioembed_sim_res,
#   sim.cutoff = 0.55,
#   cluster_method = "louvain"
# )

object$enriched_pathway@merged_module$functional_module_result <-
  object$enriched_pathway@merged_module$functional_module_result |>
  dplyr::rename(module_content_number = module_content_number.x)

object$enriched_pathway@merged_module$result_with_module <-
  object$enriched_pathway@merged_module$result_with_module |>
  dplyr::rename(module_content_number = module_content_number.x)

save(object, file = "3_data_analysis/02_control_data/05_benchmarking/mapa_fm_result.rda")

# plot <-
#   plot_similarity_network(
#     object = object$enriched_pathway,
#     level = "functional_module",
#     database = c("go", "kegg", "reactome")
#   )

library(Cairo)
CairoPDF("3_data_analysis/02_control_data/05_benchmarking/mapa_network.pdf", width = 10, height = 8)
plot
dev.off()
