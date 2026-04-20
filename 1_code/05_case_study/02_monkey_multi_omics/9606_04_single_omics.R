library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
library(org.Hs.eg.db)
source("1_code/100_tools.R")

setwd("2_data/case_study/02_monkey_multi_omics/")

# spleen up ====
load("taxid_9606/spleen/met_enriched_pathways.rda")
load("taxid_9606/spleen/up/up_rna_enriched_pathways.rda")
load("taxid_9606/spleen/up/up_prot_enriched_pathways.rda")

## mRNA
T_pth_sim <-
  get_bioembedsim(
    object = up_rna_enriched_pathways,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05
  )

save(T_pth_sim, file = "taxid_9606/spleen/up/T_pth_sim.rda")

T_fm <- get_functional_modules(
  object = T_pth_sim,
  sim.cutoff = 0.7,
  cluster_method = "louvain"
)

save(T_fm, file = "taxid_9606/spleen/up/T_fm.rda")
load("taxid_9606/spleen/up/T_fm.rda")

T_fm@merged_module$functional_module_result <-
  T_fm@merged_module$functional_module_result |>
  dplyr::filter(module %in% c("Functional_module_2"))

llm_t_res <- llm_interpret_module(
  object = T_fm,
  phenotype = "aging",
  llm_model = "gpt-4o-mini-2024-07-18",
  embedding_model = "text-embedding-3-small",
  api_key = api_key,
  embedding_output_dir = "embedding_output/",
  module_content_number_cutoff = 2,
  orgdb = org.Hs.eg.db
)

save(llm_t_res, file = "taxid_9606/spleen/up/llm_t_res.rda")

plot_similarity_network(
  object = llm_t_res,
  level = "functional_module",
  module_id = c("Functional_module_2", "Functional_module_6", "Functional_module_10",
                "Functional_module_15", "Functional_module_16", "Functional_module_17"),
  llm_text = TRUE
)

## protein
P_pth_sim <-
  get_bioembedsim(
    object = up_prot_enriched_pathways,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05
  )

save(P_pth_sim, file = "taxid_9606/spleen/up/P_pth_sim.rda")

P_fm <- get_functional_modules(
  object = P_pth_sim,
  sim.cutoff = 0.7,
  cluster_method = "louvain"
)

save(P_fm, file = "taxid_9606/spleen/up/P_fm.rda")

load("taxid_9606/spleen/up/P_fm.rda")

P_fm@merged_module$functional_module_result <-
  P_fm@merged_module$functional_module_result |>
  dplyr::filter(module %in% c("Functional_module_6", "Functional_module_10", "Functional_module_21"))

llm_p_res <- llm_interpret_module(
  object = P_fm,
  phenotype = "aging",
  llm_model = "gpt-4o-mini-2024-07-18",
  embedding_model = "text-embedding-3-small",
  api_key = api_key,
  embedding_output_dir = "embedding_output/",
  module_content_number_cutoff = 0,
  orgdb = org.Hs.eg.db
)

save(llm_p_res, file = "taxid_9606/spleen/up/llm_p_res.rda")

load("taxid_9606/spleen/up/llm_p_res.rda")

set.seed(5)

plot <-
plot_similarity_network(
  object = llm_p_res,
  level = "functional_module",
  module_id = c("Functional_module_3", "Functional_module_6", "Functional_module_7", "Functional_module_12"),
  llm_text = TRUE
)

plot

library(Cairo)
CairoPDF("taxid_9606/spleen/up/spleen_up_prot_fms.pdf", width = 8, height = 7)
plot
dev.off()

load("taxid_9606/spleen/up/llm_t_res.rda")
object <- llm_t_res
level = "functional_module"
module_id = c("Functional_module_2", "Functional_module_3")
llm_text = TRUE

load("taxid_9606/spleen/up/llm_p_res.rda")
load("taxid_9606/spleen/met_enriched_pathways.rda")

## metabolite
met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)
met_enriched_pathways
m_pth_sim <-
  get_bioembedsim(
    object = met_enriched_pathways,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("metkegg"),
    p.adjust.cutoff.metkegg = 0.05
  )

# No enriched pathways found

# Thymus down ====
load("taxid_9606/thymus/down/down_rna_enriched_pathways.rda")
load("taxid_9606/thymus/down/down_prot_enriched_pathways.rda")
load("taxid_9606/thymus/met_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

down_rna_enriched_pathways

# T_pth_sim <-
#   get_bioembedsim(
#     object = down_rna_enriched_pathways,
#     api_provider = "openai",
#     text_embedding_model = "text-embedding-3-small",
#     api_key = api_key,
#     database = c("go", "kegg"),
#     p.adjust.cutoff.go = 0.05,
#     p.adjust.cutoff.kegg = 0.05
#   )

# No enriched pathways

# T_fm <- get_functional_modules(
#   object = T_pth_sim,
#   sim.cutoff = 0.55,
#   cluster_method = "louvain"
# )

## protein
P_pth_sim <-
  get_bioembedsim(
    object = down_prot_enriched_pathways,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05,
    count.cutoff.go = 5,
    count.cutoff.kegg = 5
  )
save(P_pth_sim, file = "taxid_9606/thymus/down/P_pth_sim.rda")

P_fm <- get_functional_modules(
  object = P_pth_sim,
  sim.cutoff = 0.55,
  cluster_method = "louvain"
)

save(P_fm, file = "taxid_9606/thymus/down/P_fm.rda")

## metabolite
# met_enriched_pathways
# m_pth_sim <-
#   get_bioembedsim(
#     object = met_enriched_pathways,
#     api_provider = "openai",
#     text_embedding_model = "text-embedding-3-small",
#     api_key = api_key,
#     database = c("metkegg"),
#     p.adjust.cutoff.metkegg = 0.05
#   )

# Thyroid up ====
load("taxid_9606/thyroid/met_enriched_pathways.rda")
load("taxid_9606/thyroid/up/up_rna_enriched_pathways.rda")
load("taxid_9606/thyroid/up/up_prot_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

## mRNA
T_pth_sim <-
  get_bioembedsim(
    object = up_rna_enriched_pathways,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05
  )
save(T_pth_sim, file = "taxid_9606/thyroid/up/T_pth_sim.rda")

T_fm <- get_functional_modules(
  object = T_pth_sim,
  sim.cutoff = 0.55,
  cluster_method = "louvain"
)
save(T_fm, file = "taxid_9606/thyroid/up/T_fm.rda")

## protein
P_pth_sim <-
  get_bioembedsim(
    object = up_prot_enriched_pathways,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05
  )
save(P_pth_sim, file = "taxid_9606/thyroid/up/P_pth_sim.rda")

P_fm <- get_functional_modules(
  object = P_pth_sim,
  sim.cutoff = 0.55,
  cluster_method = "louvain"
)
save(P_fm, file = "taxid_9606/thyroid/up/P_fm.rda")

plot_module_info(
  object = P_fm,
  level = "functional_module",
  module_id = c("Functional_module_14")
)

met_enriched_pathways

# Aortic arch up ====
load("taxid_9606/aortic_arch/met_enriched_pathways.rda")
load("taxid_9606/aortic_arch/up/up_rna_enriched_pathways.rda")
load("taxid_9606/aortic_arch/up/up_prot_enriched_pathways.rda")

met_enriched_pathways@variable_info <- met_enriched_pathways@variable_info |> dplyr::rename(diff_metric = beta)

## mRNA
T_pth_sim <-
  get_bioembedsim(
    object = up_rna_enriched_pathways,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05
  )
save(T_pth_sim, file = "taxid_9606/aortic_arch/up/T_pth_sim.rda")
load("taxid_9606/aortic_arch/up/T_pth_sim.rda")
T_fm <- get_functional_modules(
  object = T_pth_sim,
  sim.cutoff = 0.7,
  cluster_method = "louvain"
)

plot_module_info(
  object = T_fm,
  level = "functional_module",
  module_id = "Functional_module_60"
)

save(T_fm, file = "taxid_9606/aortic_arch/up/T_fm.rda")

load("taxid_9606/aortic_arch/up/T_fm.rda")
T_fm@merged_module$functional_module_result <- T_fm@merged_module$functional_module_result |>
  dplyr::filter(module == c("Functional_module_60"))

llm_t_res <- llm_interpret_module(
  object = T_fm,
  phenotype = "aging",
  llm_model = "gpt-4o-mini-2024-07-18",
  embedding_model = "text-embedding-3-small",
  api_key = api_key,
  embedding_output_dir = "embedding_output/",
  module_content_number_cutoff = 2,
  orgdb = org.Hs.eg.db
)

save(llm_t_res, file = "taxid_9606/aortic_arch/up/llm_t_res.rda")

plot_module_info(
  object = llm_t_res,
  level = "functional_module",
  module_id = c("Functional_module_8")
)

## protein
P_pth_sim <-
  get_bioembedsim(
    object = up_prot_enriched_pathways,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05
  )
save(P_pth_sim, file = "taxid_9606/aortic_arch/up/P_pth_sim.rda")

load("taxid_9606/aortic_arch/up/P_pth_sim.rda")
P_fm <- get_functional_modules(
  object = P_pth_sim,
  sim.cutoff = 0.7,
  cluster_method = "louvain"
)
save(P_fm, file = "taxid_9606/aortic_arch/up/P_fm.rda")

P_fm@merged_module$functional_module_result <- P_fm@merged_module$functional_module_result |>
  dplyr::filter(module == "Functional_module_68")

llm_p_res <- llm_interpret_module(
  object = P_fm,
  phenotype = "aging",
  llm_model = "gpt-4o-mini-2024-07-18",
  embedding_model = "text-embedding-3-small",
  api_key = api_key,
  embedding_output_dir = "embedding_output/",
  module_content_number_cutoff = 2,
  orgdb = org.Hs.eg.db
)

save(llm_p_res, file = "taxid_9606/aortic_arch/up/llm_p_res.rda")
load("taxid_9606/aortic_arch/up/llm_p_res.rda")

plot_module_info(
  object = llm_p_res,
  level = "functional_module",
  module_id = c("Functional_module_68")
)


