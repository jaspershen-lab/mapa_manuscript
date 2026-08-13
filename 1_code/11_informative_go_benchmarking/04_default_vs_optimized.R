## Benchmarking on the informative GO set benchmark (k = 20)
## Step 4: Default vs optimized ARI for the four methods.
## Default configurations:
##   - MAPA:       get_functional_modules(sim.cutoff = 0.55, cluster_method = "louvain")
##   - enrichplot: emapplot() clustering defaults (nCluster = NULL ->
##                 floor(sqrt(n)) = 17, clusterFunction = stats::kmeans);
##                 showCategory raised to 313 so all terms are clustered
##                 (the pure default of 30 would only cover 30/313 terms)
##   - aPEAR:      enrichmentNetwork() defaults (jaccard, clustMethod = "markov",
##                 innerCutoff = 0.1, outerCutoff = 0.5, minClusterSize = 2)
##   - PAVER:      generate_themes() defaults (ward.D2 + dynamicTreeCut defaults)
## Optimized results are taken from the parameter scans
## (1_code/10_sensitivity_analysis_GO and 02_run_analysis.R of this folder).

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/11_informative_go_benchmarking/ground_truth_dt.rda")
load("3_data_analysis/11_informative_go_benchmarking/enriched_result.rda")

# ---- MAPA ----
## Default: actual package call with default parameters
load("3_data_analysis/10_sensitivity_res_GO/bioembed_sim_res.rda")

## Fill in the real gene sets (UniProt accessions) so that the module-result
## assembly inside merge_pathways_bioembedsim() can map gene IDs; the
## benchmark object was built with empty geneID fields
load("3_data_analysis/11_informative_go_benchmarking/pathway_genes.rda")
go_result <- bioembed_sim_res$enriched_pathway@enrichment_go_result@result
go_result$geneID <- sapply(pathway_genes[go_result$ID], paste, collapse = "/")
bioembed_sim_res$enriched_pathway@enrichment_go_result@result <- go_result

all_accessions <- unique(unlist(pathway_genes))
bioembed_sim_res$enriched_pathway@variable_info <-
  data.frame(
    ensembl = all_accessions,
    entrezid = all_accessions,
    uniprot = all_accessions,
    symbol = all_accessions,
    variable_id = paste0("gene_", seq_along(all_accessions))
  )

## Drop the template object's leftover KEGG/Reactome example results so only
## the 313 benchmark GO terms enter the module network
bioembed_sim_res$enriched_pathway@enrichment_kegg_result <- NULL
bioembed_sim_res$enriched_pathway@enrichment_reactome_result <- NULL

set.seed(23)
mapa_default_res <- get_functional_modules(
  object = bioembed_sim_res,
  sim.cutoff = 0.55,
  cluster_method = "louvain"
)

mapa_default_modules <-
  mapa_default_res@merged_module$result_with_module |>
  dplyr::select(node, module) |>
  dplyr::rename(id = node) |>
  dplyr::right_join(ground_truth_dt, by = "id")
stopifnot(!anyNA(mapa_default_modules$module))

mapa_default_ari <- mclust::adjustedRandIndex(mapa_default_modules$module,
                                              mapa_default_modules$ground_truth_label)
mapa_default_nclust <- length(unique(mapa_default_modules$module))
cat("MAPA default (louvain, sim.cutoff = 0.55): ARI =", round(mapa_default_ari, 4),
    "| Clusters:", mapa_default_nclust, "\n")

## Optimized: best of the clustering scan on MAPA's embedding similarity matrix
load("3_data_analysis/10_sensitivity_res_GO/all_ari_long.rda")
mapa_best <- all_ari_long |>
  dplyr::slice_max(order_by = ARI, n = 1, with_ties = FALSE)

# ---- enrichplot ----
library(enrichplot)
n_terms <- nrow(enriched_result@result)
edo <- enrichplot::pairwise_termsim(enriched_result, showCategory = n_terms)

set.seed(123)
p_default <- enrichplot::emapplot(
  x = edo,
  showCategory = n_terms,
  node_label = "group",
  group = TRUE
)

enrichplot_default_clusters <-
  p_default$data[, c("name", "color2")] |>
  dplyr::rename(Description = name, cluster_label = color2) |>
  dplyr::left_join(edo@result[, c("ID", "Description")], by = "Description") |>
  dplyr::select(ID, cluster_label) |>
  dplyr::rename(id = ID) |>
  dplyr::right_join(ground_truth_dt, by = "id")
stopifnot(!anyNA(enrichplot_default_clusters$cluster_label))

enrichplot_default_ari <- mclust::adjustedRandIndex(
  enrichplot_default_clusters$cluster_label,
  enrichplot_default_clusters$ground_truth_label
)
enrichplot_default_nclust <- length(unique(enrichplot_default_clusters$cluster_label))
cat("enrichplot default (nCluster = floor(sqrt(n)) =", floor(sqrt(n_terms)),
    "): ARI =", round(enrichplot_default_ari, 4),
    "| Clusters:", enrichplot_default_nclust, "\n")

load("3_data_analysis/11_informative_go_benchmarking/enrichplot_ari.rda")
enrichplot_best <- enrichplot_ari[which.max(enrichplot_ari$ARI), ]

# ---- aPEAR ----
## Default = clustMethod "markov" with the enrichmentNetwork() defaults,
## already computed in 02_run_analysis.R
load("3_data_analysis/11_informative_go_benchmarking/apear_ari.rda")
apear_default <- apear_ari[apear_ari$Method == "markov", ]
apear_best <- apear_ari[which.max(apear_ari$ARI), ]

# ---- PAVER ----
library(PAVER)
load("3_data_analysis/11_informative_go_benchmarking/PAVER_01_input.rda")
load("3_data_analysis/11_informative_go_benchmarking/PAVER_02_embedding_matrix.rda")
load("3_data_analysis/11_informative_go_benchmarking/PAVER_03_term2name.rda")

PAVER_result <- prepare_data(input, embedding_matrix, term2name)

set.seed(234)
PAVER_default_result <- generate_themes(PAVER_result)
paver_default_export <- PAVER_export(PAVER_default_result)

paver_default_clusters <-
  paver_default_export |>
  dplyr::select(GOID, Cluster) |>
  dplyr::rename(id = GOID) |>
  dplyr::right_join(ground_truth_dt, by = "id")
stopifnot(!anyNA(paver_default_clusters$Cluster))

paver_default_ari <- mclust::adjustedRandIndex(paver_default_clusters$Cluster,
                                               paver_default_clusters$ground_truth_label)
paver_default_nclust <- length(unique(paver_default_clusters$Cluster))
cat("PAVER default (generate_themes defaults): ARI =", round(paver_default_ari, 4),
    "| Clusters:", paver_default_nclust, "\n")

load("3_data_analysis/11_informative_go_benchmarking/paver_ari.rda")
paver_best <- paver_clustering_res[which.max(paver_clustering_res$ari), ]

# ---- Combine ----
default_vs_optimized_ari <- data.frame(
  methods = c("MAPA", "enrichplot", "aPEAR", "PAVER"),
  default_ari = c(mapa_default_ari,
                  enrichplot_default_ari,
                  apear_default$ARI,
                  paver_default_ari),
  default_params = c("sim.cutoff=0.55, louvain",
                     paste0("nCluster=floor(sqrt(n))=", floor(sqrt(n_terms)), ", kmeans"),
                     "jaccard, markov",
                     "ward.D2, cutreeDynamic defaults"),
  default_n_clusters = c(mapa_default_nclust,
                         enrichplot_default_nclust,
                         apear_default$N_Clusters,
                         paver_default_nclust),
  optimized_ari = c(mapa_best$ARI,
                    enrichplot_best$ARI,
                    apear_best$ARI,
                    paver_best$ari),
  optimized_params = c(paste0(mapa_best$Algorithm, ", ", mapa_best$Parameter),
                       paste0("nCluster=", enrichplot_best$nCluster, ", kmeans"),
                       paste0("jaccard, ", apear_best$Method),
                       paste0("maxCoreScatter=", paver_best$maxCoreScatter)),
  optimized_n_clusters = c(mapa_best$Num_Clusters,
                           enrichplot_best$nCluster,
                           apear_best$N_Clusters,
                           paver_best$nClust)
)

save(default_vs_optimized_ari,
     file = "3_data_analysis/11_informative_go_benchmarking/default_vs_optimized_ari.rda")
readr::write_csv(default_vs_optimized_ari,
                 "3_data_analysis/11_informative_go_benchmarking/default_vs_optimized_ari.csv")

cat("\n===== Default vs optimized ARI (k = 20 informative GO benchmark) =====\n")
print(default_vs_optimized_ari, row.names = FALSE)
