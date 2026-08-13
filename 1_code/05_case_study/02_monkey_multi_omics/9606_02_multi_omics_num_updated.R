library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
library(org.Hs.eg.db)
source("1_code/100_tools.R")

setwd("2_data/case_study/02_monkey_multi_omics/")

# Keep all results generated with the updated MAPA version separate from the
# original taxid_9606 results. Enrichment inputs are still loaded from the
# original directory below.
output_root <- "taxid_9606-multi_omics_num_updated"
output_dirs <- file.path(
  output_root,
  rep(c("thymus", "spleen", "ovary", "stomach", "kidney", "aortic_arch", "thyroid"), each = 2),
  rep(c("up", "down"), 7)
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Diagnostic plots from the original script are collected here instead of
# creating Rplots.pdf outside the updated-results directory. A stale module ID
# should not prevent the remaining tissue/direction calculations from running.
grDevices::pdf(file.path(output_root, "diagnostic_plots.pdf"))
plot_module_info_mapa <- mapa::plot_module_info
plot_module_info <- function(...) {
  tryCatch(
    {
      diagnostic_plot <- plot_module_info_mapa(...)
      print(diagnostic_plot)
      invisible(NULL)
    },
    error = function(e) warning("Skipping diagnostic module plot: ", conditionMessage(e))
  )
}

# MAPA 3.0.0 currently maps gene symbols through the live Ensembl service when
# constructing Reactome enzyme-metabolite edges. Use the equivalent local
# org.Hs.eg.db mapping so a temporary Ensembl outage cannot abort the run.
map_symbol_to_uniprot_local <- function(gene_symbols) {
  valid_symbols <- intersect(
    unique(gene_symbols),
    AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
  )
  if (length(valid_symbols) == 0) {
    return(tibble::tibble(symbol = character(), uniprotkb = character()))
  }

  AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = valid_symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "UNIPROT")
  ) |>
    dplyr::transmute(symbol = SYMBOL, uniprotkb = UNIPROT) |>
    dplyr::filter(!is.na(uniprotkb), uniprotkb != "") |>
    dplyr::distinct()
}

mapa_namespace <- asNamespace("mapa")
unlockBinding(".map_symbol_to_uniprot", mapa_namespace)
assign(".map_symbol_to_uniprot", map_symbol_to_uniprot_local, envir = mapa_namespace)
lockBinding(".map_symbol_to_uniprot", mapa_namespace)

# Thymus ====
load("taxid_9606/thymus/met_enriched_pathways.rda")
load("taxid_9606/thymus/up/up_rna_enriched_pathways.rda")
load("taxid_9606/thymus/up/up_prot_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/thymus/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/thymus/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/thymus/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.45,
                                         cluster_method = "louvain")
save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/thymus/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_46",
                 show_rwr_edge = TRUE)

# Thymus (down) ====
load("taxid_9606/thymus/down/down_rna_enriched_pathways.rda")
load("taxid_9606/thymus/down/down_prot_enriched_pathways.rda")
load("taxid_9606/thymus/met_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/thymus/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/thymus/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/thymus/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")
save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/thymus/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_5",
                 node_size = 5,
                 show_rwr_edge = TRUE)

# Spleen ====
load("taxid_9606/spleen/met_enriched_pathways.rda")
load("taxid_9606/spleen/up/up_rna_enriched_pathways.rda")
load("taxid_9606/spleen/up/up_prot_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/spleen/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/spleen/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/spleen/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/spleen/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_65")

# Spleen (down) ====
load("taxid_9606/spleen/met_enriched_pathways.rda")
load("taxid_9606/spleen/down/down_rna_enriched_pathways.rda")
load("taxid_9606/spleen/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/spleen/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/spleen/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/spleen/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/spleen/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_16")

# Ovary ====
load("taxid_9606/ovary/met_enriched_pathways.rda")
load("taxid_9606/ovary/up/up_rna_enriched_pathways.rda")
load("taxid_9606/ovary/up/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/ovary/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/ovary/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/ovary/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.45,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/ovary/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_61")

# Ovary (down) ====
load("taxid_9606/ovary/met_enriched_pathways.rda")
load("taxid_9606/ovary/down/down_rna_enriched_pathways.rda")
load("taxid_9606/ovary/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/ovary/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/ovary/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/ovary/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/ovary/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_21")

# Stomach ====
setwd(r4projects::get_project_wd())
# Keep output configuration and diagnostic plot device across sections.
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/stomach/met_enriched_pathways.rda")
load("taxid_9606/stomach/up/up_rna_enriched_pathways.rda")
load("taxid_9606/stomach/up/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/stomach/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/stomach/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/stomach/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/stomach/up/multi_omics_fm.rda")

# plot_module_info(multi_omics_fm,
#                  module_id = "Functional_module_1")

# Stomach (down) ====
setwd(r4projects::get_project_wd())
# Keep output configuration and diagnostic plot device across sections.
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/stomach/met_enriched_pathways.rda")
load("taxid_9606/stomach/down/down_rna_enriched_pathways.rda")
load("taxid_9606/stomach/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/stomach/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/stomach/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/stomach/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/stomach/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_6")

# Kidney ====
setwd(r4projects::get_project_wd())
# Keep output configuration and diagnostic plot device across sections.
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/kidney/met_enriched_pathways.rda")
load("taxid_9606/kidney/up/up_rna_enriched_pathways.rda")
load("taxid_9606/kidney/up/up_prot_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/kidney/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/kidney/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.5, delta1 = 0.3)
save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/kidney/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.5,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/kidney/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_1")

# Kidney (down) ====
setwd(r4projects::get_project_wd())
# Keep output configuration and diagnostic plot device across sections.
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/kidney/met_enriched_pathways.rda")
load("taxid_9606/kidney/down/down_rna_enriched_pathways.rda")
load("taxid_9606/kidney/down/down_prot_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/kidney/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/kidney/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.5, lambda = 0.5, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/kidney/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.45,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/kidney/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_1")

# Aortic arch ====
setwd(r4projects::get_project_wd())
# Keep output configuration and diagnostic plot device across sections.
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/aortic_arch/met_enriched_pathways.rda")
load("taxid_9606/aortic_arch/up/up_rna_enriched_pathways.rda")
load("taxid_9606/aortic_arch/up/up_prot_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/aortic_arch/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/aortic_arch/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/aortic_arch/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/aortic_arch/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_68")

# Aortic arch (down) ====
setwd(r4projects::get_project_wd())
# Keep output configuration and diagnostic plot device across sections.
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/aortic_arch/met_enriched_pathways.rda")
load("taxid_9606/aortic_arch/down/down_rna_enriched_pathways.rda")
load("taxid_9606/aortic_arch/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/aortic_arch/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/aortic_arch/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/aortic_arch/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/aortic_arch/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_8")

# Thyroid ====
setwd(r4projects::get_project_wd())
# Keep output configuration and diagnostic plot device across sections.
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/thyroid/met_enriched_pathways.rda")
load("taxid_9606/thyroid/up/up_rna_enriched_pathways.rda")
load("taxid_9606/thyroid/up/up_prot_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/thyroid/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/thyroid/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/thyroid/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/thyroid/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 show_rwr_edge = TRUE,
                 module_id = "Functional_module_24")

# Thyroid (down) ====
setwd(r4projects::get_project_wd())
# Keep output configuration and diagnostic plot device across sections.
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/thyroid/met_enriched_pathways.rda")
load("taxid_9606/thyroid/down/down_rna_enriched_pathways.rda")
load("taxid_9606/thyroid/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606-multi_omics_num_updated/thyroid/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           embedding_source = "local",
                           db_version = "v1")
save(mnet_obj, file = "taxid_9606-multi_omics_num_updated/thyroid/down/mnet_obj.rda")

load("taxid_9606-multi_omics_num_updated/thyroid/down/mnet_obj.rda")
mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5, min_path_sim = 0)
save(mnet_sim, file = "taxid_9606-multi_omics_num_updated/thyroid/down/mnet_sim.rda")
mnet_sim@x |> hist()

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606-multi_omics_num_updated/thyroid/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_1")

grDevices::dev.off()
