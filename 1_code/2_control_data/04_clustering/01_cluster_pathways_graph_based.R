library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(igraph)
library(tidygraph)

# Edge Data
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_df.rda")
# Node Attribute Data
control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

network_summary_df <- bind_rows(network_summary_list)

analyze_community_detection <- function(embedding_sim_df,
                                        control_dt,
                                        cutoff_range = seq(0.1, 0.9, by = 0.1),
                                        verbose = TRUE) {

  all_results <- list()

  for (sim.cutoff in cutoff_range) {

    if (verbose) {
      cat("Processing cutoff:", sim.cutoff, "\n")
    }

    # Filter edges based on cutoff
    edge_data_with_cutoff <-
      embedding_sim_df |>
      rename(weight = sim) |>
      dplyr::filter(weight > sim.cutoff)

    # Check if we have any edges left
    if (nrow(edge_data_with_cutoff) == 0) {
      if (verbose) {
        cat("No edges remain at cutoff", sim.cutoff, "- skipping\n")
      }
      next
    }

    # Prepare node data
    node_data <-
      control_dt |>
      dplyr::rename(node = id) |>
      dplyr::select(node, expected_module, expected_count, database, name)

    # Create graph
    graph_data_with_cutoff <-
      tbl_graph(
        nodes = node_data,
        edges = edge_data_with_cutoff,
        directed = FALSE,
        node_key = "node"
      ) |>
      mutate(degree = centrality_degree())

    graph_data <- graph_data_with_cutoff

    # Initialize community results for this cutoff
    community_results <- list()

    ## Fast Greedy
    fast_greedy_res <- tryCatch(cluster_fast_greedy(graph_data), error = function(e) NA)
    community_results[["fast_greedy"]] <- fast_greedy_res

    ## Louvain
    louvain_res <- tryCatch(cluster_louvain(graph_data), error = function(e) NA)
    community_results[["louvain"]] <- louvain_res

    ## Walktrap
    walktrap_res <- tryCatch(cluster_walktrap(graph_data), error = function(e) NA)
    community_results[["walktrap"]] <- walktrap_res

    ## InfoMap
    infomap_res <- tryCatch(cluster_infomap(graph_data), error = function(e) NA)
    community_results[["infomap"]] <- infomap_res

    ## Label Prop
    label_prop_res <- tryCatch(cluster_label_prop(graph_data), error = function(e) NA)
    community_results[["label_prop"]] <- label_prop_res

    ## Leading Eigen
    leading_eigen_res <- tryCatch(cluster_leading_eigen(graph_data), error = function(e) NA)
    community_results[["leading_eigen"]] <- leading_eigen_res

    ## Leiden
    leiden_res <- tryCatch(cluster_leiden(graph_data), error = function(e) NA)
    community_results[["leiden"]] <- leiden_res

    ## Edge Betweenness
    distance_weights <- 1 - E(graph_data)$weight
    edge_betweenness_res <- tryCatch(cluster_edge_betweenness(graph_data, weights = distance_weights), error = function(e) NA)
    community_results[["edge_betweenness"]] <- edge_betweenness_res

    ## Optimal
    optimal_res <- tryCatch(cluster_optimal(graph_data), error = function(e) NA)
    community_results[["optimal"]] <- optimal_res

    ## Spinglass
    spinglass_res <- tryCatch(cluster_spinglass(graph_data), error = function(e) NA)
    community_results[["spinglass"]] <- spinglass_res

    ## Fluid Communities
    best_modularity <- -1
    fluid_community_best_result <- NULL
    max_k <- min(43, vcount(graph_data))  # Don't exceed number of nodes

    for (k in 2:max_k) {
      temp_result <- try(cluster_fluid_communities(graph_data, no.of.communities = k), silent = TRUE)
      if (!inherits(temp_result, "try-error")) {
        mod <- modularity(x = graph_data, membership = membership(temp_result))
        if (mod > best_modularity) {
          best_modularity <- mod
          fluid_community_best_result <- temp_result
        }
      }
    }

    if (is.null(fluid_community_best_result)) {
      fluid_community_best_result <- NA
    }
    community_results[["fluid_community"]] <- fluid_community_best_result

    # Prepare ground truth
    ground_truth <-
      node_data |>
      dplyr::mutate(ground_truth = as.numeric(sub(pattern = "Functional_module_", replacement = "", x = expected_module))) |>
      pull(ground_truth)
    names(ground_truth) <- node_data$name

    # Calculate summary statistics for all algorithms
    network_summary_list <- lapply(names(community_results), function(algo_name) {
      tryCatch(
        expr = {
          comm <- community_results[[algo_name]]
          if (is.na(comm)[1]) {
            return(data.frame(
              Cutoff = sim.cutoff,
              Algorithm = algo_name,
              Type = "Network-Based",
              Num_Clusters = NA,
              Modularity = NA,
              ARI = NA,
              RI = NA,
              NMI = NA,
              VI = NA,
              SplitJoin = NA
            ))
          }

          detected_membership <- membership(comm)
          data.frame(
            Cutoff = sim.cutoff,
            Algorithm = algo_name,
            Type = "Network-Based",
            Num_Clusters = length(unique(detected_membership)),
            Modularity = modularity(x = graph_data, membership = detected_membership),
            ARI = compare(detected_membership, ground_truth, method = "adjusted.rand"),
            RI = compare(detected_membership, ground_truth, method = "rand"),
            NMI = compare(detected_membership, ground_truth, method = "nmi"),
            VI = compare(detected_membership, ground_truth, method = "vi"),
            SplitJoin = compare(detected_membership, ground_truth, method = "split.join")
          )
        },
        error = function(e) {
          if (verbose) {
            message("Error with algorithm ", algo_name, " at cutoff ", sim.cutoff, ": ", e$message)
          }
          data.frame(
            Cutoff = sim.cutoff,
            Algorithm = algo_name,
            Type = "Network-Based",
            Num_Clusters = NA,
            Modularity = NA,
            ARI = NA,
            RI = NA,
            NMI = NA,
            VI = NA,
            SplitJoin = NA
          )
        }
      )
    })

    network_summary_df <- bind_rows(network_summary_list)
    all_results[[as.character(sim.cutoff)]] <- network_summary_df
  }

  # Combine all results
  final_results <- bind_rows(all_results)

  return(final_results)
}

results_0.01 <- analyze_community_detection(
  embedding_sim_df = embedding_sim_df,
  control_dt = control_dt,
  cutoff_range = seq(0.2, 0.9, by = 0.01),
  verbose = TRUE
)

save(results_0.01, file = "3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_graph_based/results_0.01.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_graph_based/results_0.01.rda")



