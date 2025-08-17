library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

met_data <- read_csv("2_data/case_study/01_monkey/m_dt_u_d_kegg_manual.csv")

dir.create("3_data_analysis/05_case_study/01_monkey/plasma_metabolomics/")
setwd("3_data_analysis/05_case_study/01_monkey/plasma_metabolomics/")

# only 9 metabolites in cluster U
met_data_u <- met_data |> filter(is.na(kegg_manual) & Cluster == "Cluster U")
# only 32 metabolites in cluster D -> No enriched pathways
met_data_d <- met_data |> filter(is.na(kegg_manual) & Cluster == "Cluster D")

met_ud <- met_data |> filter(Cluster %in% c("Cluster U", "Cluster D") & !is.na(keggid))

## Put them together -> dysregulated metabolites
enrich_pathway_res <- enrich_pathway(
  variable_info = met_ud,
  query_type = "metabolite",
  met_organism = "mcf",
  database = "metkegg"
)

enrich_pathway_res

save(enrich_pathway_res, file = "met_dys_enrich_pathway_res.rda")

embedsim_res <-
  get_bioembedsim(
    object = enrich_pathway_res,
    api_provider = "openai",
    api_key = api_key,
    text_embedding_model = "text-embedding-3-small",
    database = c("metkegg"),
    count.cutoff.metkegg = 0
  )

save(embedsim_res, file = "met_dys_embedsim_res.rda")

load("met_dys_embedsim_res.rda")

functional_module_res <- get_functional_modules(
  object = embedsim_res,
  sim.cutoff = 0.55,
  cluster_method = "louvain"
)

save(functional_module_res, file = "functional_module_res_u_and_d_met.rda")

functional_module <- functional_module_res@merged_module$functional_module_result
export(functional_module, file = "met_u_and_d_functional_module_result.csv")

plot <- plot_similarity_network(
  object = functional_module_res,
  level = "functional_module",
  database = c("metkegg"),
  degree_cutoff = 0,
  llm_text = FALSE,
  text_all = TRUE
)

plot
