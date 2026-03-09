library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
library(org.Hs.eg.db)
source("1_code/100_tools.R")

setwd("2_data/case_study/02_monkey_multi_omics/")

# Thymus ====
load("taxid_9606/thymus/met_enriched_pathways.rda")
load("taxid_9606/thymus/up/up_rna_enriched_pathways.rda")
load("taxid_9606/thymus/up/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/thymus/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/thymus/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606/thymus/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")
save(multi_omics_fm, file = "taxid_9606/thymus/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_16")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_68")


# Thymus (down) ====
load("taxid_9606/thymus/down/down_rna_enriched_pathways.rda")
load("taxid_9606/thymus/down/down_prot_enriched_pathways.rda")
load("taxid_9606/thymus/met_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/thymus/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/thymus/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606/thymus/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")
save(multi_omics_fm, file = "taxid_9606/thymus/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_4")

# Spleen ====
load("taxid_9606/spleen/met_enriched_pathways.rda")
load("taxid_9606/spleen/up/up_rna_enriched_pathways.rda")
load("taxid_9606/spleen/up/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/spleen/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/spleen/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606/spleen/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/spleen/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_51")

# Spleen (down) ====
load("taxid_9606/spleen/met_enriched_pathways.rda")
load("taxid_9606/spleen/down/down_rna_enriched_pathways.rda")
load("taxid_9606/spleen/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/spleen/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/spleen/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606/spleen/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/spleen/down/multi_omics_fm.rda")

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
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/ovary/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/ovary/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606/ovary/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.45,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/ovary/up/multi_omics_fm.rda")

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
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/ovary/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/ovary/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606/ovary/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/ovary/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_21")

# Stomach ====
load("stomach/met_enriched_pathways.rda")
load("stomach/up_rna_enriched_pathways.rda")
load("stomach/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9544,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.9,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "stomach/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "ovary/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "ovary/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "ovary/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_22")

# Kidney ====
load("kidney/met_enriched_pathways.rda")
load("kidney/up_rna_enriched_pathways.rda")
load("kidney/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9544,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.9,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "kidney/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "kidney/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.3)
save(mnet_sim, file = "kidney/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "kidney/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_22")

# Aortic arch ====
load("aortic_arch/met_enriched_pathways.rda")
load("aortic_arch/up_rna_enriched_pathways.rda")
load("aortic_arch/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9544,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.9,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "aortic_arch/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "aortic_arch/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.3)
save(mnet_sim, file = "aortic_arch/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "aortic_arch/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_22")

# Aortic arch (down) ====
load("aortic_arch/met_enriched_pathways.rda")
load("aortic_arch_down/down_rna_enriched_pathways.rda")
load("aortic_arch_down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9544,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.9,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "aortic_arch_down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "aortic_arch_down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.3)
save(mnet_sim, file = "aortic_arch_down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "aortic_arch_down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_7")

# Thyroid ====
load("thyroid/met_enriched_pathways.rda")
load("thyroid/up_rna_enriched_pathways.rda")
load("thyroid/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9544,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.9,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "thyroid/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "thyroid/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.3)
save(mnet_sim, file = "thyroid/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "thyroid/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_11")


# Thyroid (down) ====
load("thyroid/met_enriched_pathways.rda")
load("thyroid_down/down_rna_enriched_pathways.rda")
load("thyroid_down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9544,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.9,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "thyroid_down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "thyroid_down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.3)
save(mnet_sim, file = "aortic_arch_down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "aortic_arch_down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_7")

