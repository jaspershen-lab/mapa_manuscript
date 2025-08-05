library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# gene info
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/all_variable_info.rda")
# functional module info
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/all_525_result_with_module_info.rda")
# meta fuctional module info
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/shared_tissues_fm/meta_module_llm_interpreted_result.rda")

mapped_genes <- all_info |> filter(tissue == "Intestine (Jejunum)" & fm_module == "Functional_module_171") |> pull(mapped_molecules)
mapped_genes <- strsplit(mapped_genes, split = "/")[[1]]
genes_info <- all_variable_info |> filter(Tissue == "Intestine (Jejunum)" & ensembl %in% mapped_genes)

mapped_genes_2 <- all_info |> filter(tissue == "Intestine (Duodenum)" & fm_module == "Functional_module_171") |> pull(mapped_molecules)
mapped_genes_2 <- strsplit(mapped_genes_2, split = "/")[[1]]
genes_info_2 <- all_variable_info |> filter(Tissue == "Intestine (Duodenum)" & ensembl %in% mapped_genes_2)

all_mapped_genes <- genes_info |> full_join(genes_info_2, by = "symbol")

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/embedding_res_object/embedsim_res_Intestine__Duodenum_.rda")

functional_module_res <- get_functional_modules(
  object = embedsim_res,
  sim.cutoff = 0.55,
  cluster_method = "louvain"
)
