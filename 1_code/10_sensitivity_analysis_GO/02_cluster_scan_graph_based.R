## Sensitivity analysis on the informative GO set benchmark (k = 20)
## Step 2: Scan similarity cutoffs x graph-based community detection
## algorithms and evaluate against the ground truth modules.
## Adapted from 1_code/2_control_data/04_clustering/01_cluster_pathways_graph_based.R

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(igraph)
library(tidygraph)

# Edge data (embedding cosine similarity)
load("3_data_analysis/10_sensitivity_res_GO/embedding_sim_df.rda")

# Node data + ground truth
mapa_input <- readr::read_csv("2_data/informative_go_set/k_20/mapa_input.csv",
                              show_col_types = FALSE)
ground_truth_dt <- readr::read_csv("2_data/informative_go_set/k_20/ground_truth.csv",
                                   show_col_types = FALSE)
ground_truth_dt <-
  ground_truth_dt |>
  dplyr::mutate(ground_truth_label = as.numeric(factor(module_id)))

node_data <-
  mapa_input |>
  dplyr::rename(node = go_id) |>
  dplyr::left_join(ground_truth_dt[, c("go_id", "ground_truth_label")],
                   by = c("node" = "go_id")) |>
  dplyr::select(node, name, ground_truth_label)

stopifnot(!anyNA(node_data$ground_truth_label))

analyze_community_detection <- function(embedding_sim_df,
                                        node_data,
                                        cutoff_range = seq(0.2, 0.9, by = 0.01),
                                        max_k_fluid = 100,
                                        edge_betweenness_max_edges = 8000,
                                        verbose = TRUE) {

  all_results <- list()

  for (sim.cutoff in cutoff_range) {

    if (verbose) {
      cat("Processing cutoff:", sim.cutoff, "\n")
    }

    # Filter edges based on cutoff
    edge_data_with_cutoff <-
      embedding_sim_df |>
      dplyr::rename(weight = sim) |>
      dplyr::filter(weight > sim.cutoff)

    if (nrow(edge_data_with_cutoff) == 0) {
      if (verbose) {
        cat("No edges remain at cutoff", sim.cutoff, "- skipping\n")
      }
      next
    }

    # Create graph (all benchmark terms as nodes, including isolated ones)
    graph_data <-
      tbl_graph(
        nodes = node_data,
        edges = edge_data_with_cutoff,
        directed = FALSE,
        node_key = "node"
      ) |>
      mutate(degree = centrality_degree())

    community_results <- list()

    ## Fast Greedy
    community_results[["fast_greedy"]] <-
      tryCatch(cluster_fast_greedy(graph_data), error = function(e) NA)

    ## Louvain
    community_results[["louvain"]] <-
      tryCatch(cluster_louvain(graph_data), error = function(e) NA)

    ## Walktrap
    community_results[["walktrap"]] <-
      tryCatch(cluster_walktrap(graph_data), error = function(e) NA)

    ## InfoMap
    community_results[["infomap"]] <-
      tryCatch(cluster_infomap(graph_data), error = function(e) NA)

    ## Label Prop
    community_results[["label_prop"]] <-
      tryCatch(cluster_label_prop(graph_data), error = function(e) NA)

    ## Leading Eigen
    community_results[["leading_eigen"]] <-
      tryCatch(cluster_leading_eigen(graph_data), error = function(e) NA)

    ## Leiden
    community_results[["leiden"]] <-
      tryCatch(cluster_leiden(graph_data), error = function(e) NA)

    ## Edge Betweenness (Girvan-Newman); with 313 nodes this is only feasible
    ## on sparser graphs, so skip when the graph is too dense
    if (ecount(graph_data) <= edge_betweenness_max_edges) {
      distance_weights <- 1 - E(graph_data)$weight
      community_results[["edge_betweenness"]] <-
        tryCatch(cluster_edge_betweenness(graph_data, weights = distance_weights),
                 error = function(e) NA)
    } else {
      community_results[["edge_betweenness"]] <- NA
    }

    ## Optimal: skipped — exact modularity optimization is computationally
    ## infeasible for 313 nodes (was only run on the 44-pathway control data)

    ## Spinglass (only works on connected graphs)
    community_results[["spinglass"]] <-
      tryCatch(cluster_spinglass(graph_data), error = function(e) NA)

    ## Fluid Communities (select k by modularity, as in the original scan)
    best_modularity <- -1
    fluid_community_best_result <- NULL
    max_k <- min(max_k_fluid, vcount(graph_data))

    for (k in 2:max_k) {
      temp_result <- try(cluster_fluid_communities(graph_data, no.of.communities = k),
                         silent = TRUE)
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

    # Ground truth aligned with node order
    ground_truth <- node_data$ground_truth_label
    names(ground_truth) <- node_data$node

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
            message("Error with algorithm ", algo_name, " at cutoff ",
                    sim.cutoff, ": ", e$message)
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

    all_results[[as.character(sim.cutoff)]] <- bind_rows(network_summary_list)
  }

  bind_rows(all_results)
}

set.seed(23)
graph_based_results <- analyze_community_detection(
  embedding_sim_df = embedding_sim_df,
  node_data = node_data,
  cutoff_range = seq(0.2, 0.9, by = 0.01),
  verbose = TRUE
)

save(graph_based_results,
     file = "3_data_analysis/10_sensitivity_res_GO/graph_based_results.rda")

## Best result per algorithm
graph_based_results |>
  dplyr::filter(!is.na(ARI)) |>
  dplyr::group_by(Algorithm) |>
  dplyr::slice_max(ARI, n = 1, with_ties = FALSE) |>
  dplyr::arrange(dplyr::desc(ARI)) |>
  print(n = 20)
