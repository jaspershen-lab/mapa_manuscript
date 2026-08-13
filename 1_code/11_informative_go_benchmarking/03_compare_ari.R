## Benchmarking on the informative GO set benchmark (k = 20)
## Step 3: Compare MAPA's optimized clustering result with the best results
## of enrichplot, aPEAR and PAVER.
## Adapted from 1_code/2_control_data/05_benchmarking/03_compare_ari.R

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# mapa best result ===
## Best ARI over the clustering scan on MAPA's biotext embedding similarity
## matrix (1_code/10_sensitivity_analysis_GO)
load("3_data_analysis/10_sensitivity_res_GO/all_ari_long.rda")
mapa_best_result <-
  all_ari_long |>
  dplyr::slice_max(order_by = ARI, n = 1, with_ties = TRUE)
print(mapa_best_result)
mapa_best_ari <- unique(mapa_best_result$ARI)

# enrichplot best result ====
load("3_data_analysis/11_informative_go_benchmarking/enrichplot_ari.rda")
enrichplot_best_result <- enrichplot_ari[which.max(enrichplot_ari$ARI), ]
print(enrichplot_best_result)
enrichplot_best_ari <- enrichplot_best_result$ARI

# aPEAR best result ====
load("3_data_analysis/11_informative_go_benchmarking/apear_ari.rda")
apear_best_result <- apear_ari[which.max(apear_ari$ARI), ]
print(apear_best_result)
apear_best_ari <- apear_best_result$ARI

# PAVER best result ===
load("3_data_analysis/11_informative_go_benchmarking/paver_ari.rda")
paver_best_result <- paver_clustering_res[which.max(paver_clustering_res$ari), ]
print(paver_best_result)
paver_best_ari <- paver_best_result$ari

all_best_ari <- data.frame(methods = c("MAPA", "enrichplot", "aPEAR", "PAVER"),
                           best_ari = c(mapa_best_ari,
                                        enrichplot_best_ari,
                                        apear_best_ari,
                                        paver_best_ari))
save(all_best_ari,
     file = "3_data_analysis/11_informative_go_benchmarking/4_methods_all_best_ari.rda")
readr::write_csv(all_best_ari,
                 "3_data_analysis/11_informative_go_benchmarking/4_methods_all_best_ari.csv")

cat("\n===== Best ARI of the 4 methods (k = 20 informative GO benchmark) =====\n")
print(all_best_ari)

# Bar plot ====
names(methods) <- c("MAPA", "aPEAR", "PAVER", "enrichplot")

p <- ggplot(all_best_ari, aes(x = reorder(methods, best_ari, decreasing = TRUE), y = best_ari, fill = methods)) +
  geom_col() +
  scale_fill_manual(values = methods) +
  geom_text(aes(label = round(best_ari, 3)), vjust = -0.3) +
  labs(x = "Methods",
       y = "Adjusted Rand Index (ARI)") +
  theme_bw()
p

## Barplot kept as an alternative view; the main figure (compare_best_ari.pdf)
## is the default-vs-optimized arrow plot from 05_plot_default_vs_optimized.R
ggsave(plot = p,
       filename = "3_data_analysis/11_informative_go_benchmarking/compare_best_ari_barplot.pdf",
       width = 8, height = 6)
