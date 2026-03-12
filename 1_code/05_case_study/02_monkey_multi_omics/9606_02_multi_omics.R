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
                 module_id = "Functional_module_65",
                 node_size = 5)

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
                 module_id = "Functional_module_65")

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
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/stomach/met_enriched_pathways.rda")
load("taxid_9606/stomach/up/up_rna_enriched_pathways.rda")
load("taxid_9606/stomach/up/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/stomach/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/stomach/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606/stomach/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/stomach/up/multi_omics_fm.rda")

# plot_module_info(multi_omics_fm,
#                  module_id = "Functional_module_1")

# Stomach (down) ====
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/stomach/met_enriched_pathways.rda")
load("taxid_9606/stomach/down/down_rna_enriched_pathways.rda")
load("taxid_9606/stomach/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/stomach/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/stomach/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)

save(mnet_sim, file = "taxid_9606/stomach/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/stomach/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_6")

# Kidney ====
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/kidney/met_enriched_pathways.rda")
load("taxid_9606/kidney/up/up_rna_enriched_pathways.rda")
load("taxid_9606/kidney/up/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/kidney/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/kidney/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606/kidney/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/kidney/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_1")

# Kidney (down) ====
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/kidney/met_enriched_pathways.rda")
load("taxid_9606/kidney/down/down_rna_enriched_pathways.rda")
load("taxid_9606/kidney/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/kidney/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/kidney/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606/kidney/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/kidney/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_1")

# Aortic arch ====
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/aortic_arch/met_enriched_pathways.rda")
load("taxid_9606/aortic_arch/up/up_rna_enriched_pathways.rda")
load("taxid_9606/aortic_arch/up/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/aortic_arch/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/aortic_arch/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606/aortic_arch/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/aortic_arch/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_68")

# Aortic arch (down) ====
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/aortic_arch/met_enriched_pathways.rda")
load("taxid_9606/aortic_arch/down/down_rna_enriched_pathways.rda")
load("taxid_9606/aortic_arch/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/aortic_arch/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/aortic_arch/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606/aortic_arch/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/aortic_arch/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_8")

# Thyroid ====
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/thyroid/met_enriched_pathways.rda")
load("taxid_9606/thyroid/up/up_rna_enriched_pathways.rda")
load("taxid_9606/thyroid/up/up_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = up_rna_enriched_pathways,
                                       proteome_enrich = up_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/thyroid/up/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/thyroid/up/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606/thyroid/up/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/thyroid/up/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_24")


# Thyroid (down) ====
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/")

load("taxid_9606/thyroid/met_enriched_pathways.rda")
load("taxid_9606/thyroid/down/down_rna_enriched_pathways.rda")
load("taxid_9606/thyroid/down/down_prot_enriched_pathways.rda")

network_tables <- build_network_tables(transcriptome_enrich = down_rna_enriched_pathways,
                                       proteome_enrich = down_prot_enriched_pathways,
                                       metabolome_enrich = met_enriched_pathways,
                                       taxon_id = 9606,
                                       reactome_dir = "reactome_db",
                                       input_directory = "string_db",
                                       string_score_cutoff = 0.8,
                                       tf_confidence_levels = c("A"))
save(network_tables, file = "taxid_9606/thyroid/down/network_tables.rda")

mnet_obj <- build_MNetwork(network_tables = network_tables,
                           api_provider = "openai",
                           text_embedding_model = "text-embedding-3-small",
                           api_key = api_key)
save(mnet_obj, file = "taxid_9606/thyroid/down/mnet_obj.rda")

mnet_sim <- get_multi_omics_sim(mnet_obj = mnet_obj,
                                r = 0.7, lambda = 0.6, delta1 = 0.5)
save(mnet_sim, file = "taxid_9606/thyroid/down/mnet_sim.rda")

multi_omics_fm <- get_functional_modules(object = mnet_obj,
                                         sim_matrix = mnet_sim,
                                         sim_cutoff = 0.55,
                                         cluster_method = "louvain")

save(multi_omics_fm, file = "taxid_9606/thyroid/down/multi_omics_fm.rda")

plot_module_info(multi_omics_fm,
                 module_id = "Functional_module_1")

