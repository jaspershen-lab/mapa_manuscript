## Sensitivity analysis on the informative GO set benchmark (k = 20)
## Step 4: Combine graph-based and distance-based scan results and report
## the best ARI overall and per algorithm.

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

res_dir <- "3_data_analysis/10_sensitivity_res_GO"

load(file.path(res_dir, "graph_based_results.rda"))
load(file.path(res_dir, "kmeans_results_df.rda"))
load(file.path(res_dir, "hcluster_results.rda"))
load(file.path(res_dir, "binary_cut_result_df.rda"))
load(file.path(res_dir, "hdb_result_df.rda"))
load(file.path(res_dir, "ap_results_df.rda"))
load(file.path(res_dir, "ms_result_df.rda"))
load(file.path(res_dir, "gmm_result_df.rda"))

all_ari_long <- dplyr::bind_rows(
  ## Graph-based: parameter is the similarity cutoff
  graph_based_results |>
    dplyr::transmute(
      Class = "Graph_based",
      Algorithm = Algorithm,
      Parameter = paste0("sim.cutoff=", Cutoff),
      ARI = ARI,
      Num_Clusters = Num_Clusters
    ),
  ## Distance-based
  kmeans_results_df |>
    dplyr::transmute(
      Class = "Distance_based",
      Algorithm = "kmeans",
      Parameter = paste0("centers=", n_centers),
      ARI = ari_score,
      Num_Clusters = Num_Clusters
    ),
  hcluster_results |>
    dplyr::transmute(
      Class = "Distance_based",
      Algorithm = paste0("hclust_", Method),
      Parameter = paste0("k=", K),
      ARI = ARI,
      Num_Clusters = K
    ),
  binary_cut_result_df |>
    dplyr::transmute(
      Class = "Distance_based",
      Algorithm = "binary_cut",
      Parameter = paste0("cutoff=", cutoff),
      ARI = ari_score,
      Num_Clusters = Num_Clusters
    ),
  hdb_result_df |>
    dplyr::transmute(
      Class = "Distance_based",
      Algorithm = "hdbscan",
      Parameter = paste0("minPts=", minpts),
      ARI = ari_score_with_noise,
      Num_Clusters = NumClusters
    ),
  ap_results_df |>
    dplyr::transmute(
      Class = "Distance_based",
      Algorithm = "affinity_propagation",
      Parameter = paste0("p=", p),
      ARI = ari_score,
      Num_Clusters = Num_Clusters
    ),
  ms_result_df |>
    dplyr::transmute(
      Class = "Distance_based",
      Algorithm = "mean_shift",
      Parameter = paste0("bandwidth=", bandwidth),
      ARI = ari_score,
      Num_Clusters = Num_Clusters
    ),
  gmm_result_df |>
    dplyr::transmute(
      Class = "Distance_based",
      Algorithm = "gmm",
      Parameter = paste0("G=", G),
      ARI = ari_score,
      Num_Clusters = Num_Clusters
    )
)

save(all_ari_long, file = file.path(res_dir, "all_ari_long.rda"))
readr::write_csv(all_ari_long, file.path(res_dir, "all_ari_long.csv"))

## Best per algorithm
best_ari_per_algorithm <-
  all_ari_long |>
  dplyr::filter(!is.na(ARI)) |>
  dplyr::group_by(Class, Algorithm) |>
  dplyr::slice_max(ARI, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::arrange(dplyr::desc(ARI))

save(best_ari_per_algorithm, file = file.path(res_dir, "best_ari_per_algorithm.rda"))
readr::write_csv(best_ari_per_algorithm, file.path(res_dir, "best_ari_per_algorithm.csv"))

cat("\n===== Best ARI per algorithm (k = 20 informative GO benchmark) =====\n")
print(as.data.frame(best_ari_per_algorithm), row.names = FALSE)

cat("\n===== Best overall =====\n")
print(as.data.frame(best_ari_per_algorithm[1, ]), row.names = FALSE)
