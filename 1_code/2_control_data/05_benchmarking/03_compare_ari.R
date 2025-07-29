library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# mapa best result ===
load("3_data_analysis/02_control_data/04_clustering/all_results.rda")
load("3_data_analysis/02_control_data/04_clustering/01_cluster_pathways_graph_based_result/results_0.01.rda")
mapa_best_result <-
  all_results |>
  slice_max(order_by = ari_score, n = 1, with_ties = TRUE)
mapa_best_ari <- unique(mapa_best_result$ari_score)

mapa_best_params <-
  results_0.01 |>
  filter(Algorithm == "louvain" & Num_Clusters == 14)
# mapa best params: sim.cutoff = 0.55, nCluster = 14, louvain

# enrichmap best result ====
load("3_data_analysis/02_control_data/05_benchmarking/enrichplot_ari.rda")
enrichplot_best_result <- enrichplot_ari[which.max(enrichplot_ari$ARI), ]
enrichplot_best_ari <- enrichplot_best_result$ARI
# enrichplot best params: nCluster = 17, kmeans

# aPEAR best result ====
load("3_data_analysis/02_control_data/05_benchmarking/apear_ari.rda")
apear_best_result <- apear_ari[which.max(apear_ari$ARI), ]
apear_best_ari <- apear_best_result$ARI
# aPEAR best params: nCluster = 25, hierarchical clustering

# PAVER best result ===
load("3_data_analysis/02_control_data/05_benchmarking/paver_ari.rda")
paver_best_result <- paver_clustering_res[which.max(paver_clustering_res$ari), ]
paver_best_ari <- paver_best_result$ari

all_best_ari <- data.frame(methods = c("MAPA", "enrichmap", "aPEAR", "PAVER"),
                           best_ari = c(mapa_best_ari,
                                        enrichplot_best_ari,
                                        apear_best_ari,
                                        paver_best_ari))
save(all_best_ari, file = "3_data_analysis/02_control_data/05_benchmarking/4_methods_all_best_ari.rda")

load("3_data_analysis/02_control_data/05_benchmarking/4_methods_all_best_ari.rda")
setwd("3_data_analysis/02_control_data/05_benchmarking/comparison_result/")

names(methods) <- c("MAPA", "aPEAR", "PAVER", "enrichplot")
methods
all_best_ari$methods[2] <- "enrichplot"

p <- ggplot(all_best_ari, aes(x = reorder(methods, best_ari, decreasing = TRUE), y = best_ari, fill = methods)) +
  geom_col() +
  scale_fill_manual(values = methods) +
  geom_text(aes(label = round(best_ari, 3)), vjust = -0.3) +
  labs(x = "Methods",
       y = "Adjusted Rand Index (ARI)") +
  theme_bw()
p

ggsave(plot = p, filename = "compare_best_ari.pdf",
       width = 8, height = 6)
