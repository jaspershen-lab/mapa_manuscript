library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("2_data/case_study/01_monkey/02_proteomics/p_dt_cleaned_down.rda")

setwd("3_data_analysis/05_case_study/01_monkey/")

## Step1: data input and preprocessing ====
### Create an AnnotationHub object to access the database
# ah <- AnnotationHub::AnnotationHub()
# ### Query the object for Macaca fascicularis
# mf_query_result <- AnnotationHub::query(ah, c("OrgDb", "9541"))
# mf_query_result$ah_id
# ### Get the orgDb object
# mf.orgdb <- ah[["AH119900"]] #| Taxonomy ID: 9541

p_dt_cleaned_down <- p_dt_cleaned_down |> dplyr::rename(symbol = ID)
variable_info <- convert_id(data = p_dt_cleaned_down,
                            query_type = "gene",
                            from_id_type = "symbol",
                            ah_id = "AH119900",
                            return_orgdb = TRUE)
variable_info_data <- variable_info$data
save(variable_info_data, file = "02_proteomics_down_variable_info.rda")

## Step2: enrichment analysis ====
enrich_pathway_res <- enrich_pathway(
  variable_info = variable_info$data,
  query_type = "gene",
  database = c("go", "kegg"),
  save_to_local = FALSE,
  go.orgdb = variable_info$orgdb,
  go.keytype = "SYMBOL",
  go.ont = "ALL",
  kegg.organism = "mcf", # Taxid:9541
  kegg.keytype = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

save(enrich_pathway_res,
     file = "02_proteomics_down_enrich_pathway_res.rda")

# only 10 pathways with p.adjust < 0.05 and count > 5

## Step3: Similarity calculation and clustering ====
# load("01_bulk_rna_seq_up_enrich_pathway_res_cluster_u.rda")
### 3.1 Traditional similarity method ====
# merged_pathways <- merge_pathways(
#   object = enrich_pathway_res,
#   database = c("go", "kegg"),
#   sim.cutoff.go = 0.5,
#   sim.cutoff.kegg = 0.5,
#   measure.method.go = "Sim_Lin_1998",
#   go.orgdb = variable_info$orgdb,
#   measure.method.kegg = "jaccard"
# )
#
# save(merged_pathways, file = "01_bulk_rna_seq_up_merged_pathways.rda")
#
# load("01_bulk_rna_seq_up_merged_pathways.rda")
# optimal_clustering_test <- determine_optimal_clusters(
#   object = merged_pathways,
#   cutoff_range = c(0.1, 0.9),
#   cutoff_increment = 0.05,
#   methods = c("h_ward.D", "h_ward.D2", "h_complete",
#               "binary_cut", "louvain", "walktrap", "edge_betweenness",
#               "fast_greedy")
# )
#
# save(optimal_clustering_test, file = "01_bulk_rna_seq_up_optimal_clustering_test.rda")
# optimal_clustering_test$evaluation_plot
# optimal_clustering_test$best_combination # use louvain, cutoff: 0.45
#
# functional_module_res <- get_functional_modules(
#   object = merged_pathways,
#   sim.cutoff = 0.45,
#   cluster_method = "louvain"
# )
#
# save(functional_module_res, file = "01_bulk_rna_seq_up_functional_module_res.rda")

### 3.2 Biotext embedding method ====
embedsim_res <-
  get_bioembedsim(
    object = enrich_pathway_res,
    api_provider = "siliconflow",
    api_key = api_key,
    text_embedding_model = "Qwen/Qwen3-Embedding-8B",
    database = c("go", "kegg"),
    p.adjust.cutoff.go = 0.05,
    p.adjust.cutoff.kegg = 0.05,
    count.cutoff.go = 5
  )

save(embedsim_res, file = "02_proteomics_down_embedsim_res.rda")

optimal_clustering_test <- determine_optimal_clusters(
  object = embedsim_res,
  cutoff_range = c(0.1, 0.9),
  cutoff_increment = 0.02,
  methods = c("h_ward.D", "h_ward.D2", "h_complete",
              "binary_cut", "louvain", "walktrap",
              "fast_greedy")
)

save(optimal_clustering_test, file = "02_proteomics_down_optimal_clustering_test_biotext.rda")
optimal_clustering_test$evaluation_plot
optimal_clustering_test$best_combination # use louvain, cutoff: 0.88

# optimal: hierarchical_ward.D2 0.85 | louvain 0.8
# functional_module_res <- get_functional_modules(
#   object = embedsim_res,
#   sim.cutoff = 0.88,
#   cluster_method = "louvain"
# )
load("02_proteomics_down_embedsim_res.rda")

functional_module_res <- get_functional_modules(
  object = embedsim_res,
  sim.cutoff = 0.55,
  cluster_method = "h_ward.D2"
)

clustering_quality <- assess_clustering_quality(functional_module_res)
clustering_quality$size_plot

clustering_quality$evaluation_plot

save(functional_module_res, file = "02_proteomics_down_functional_module_res_biotext.rda")

## Step4: LLM interpretation ====
setwd("3_data_analysis/05_case_study/01_monkey/")
load("02_proteomics_down_functional_module_res_biotext.rda")
ah <- AnnotationHub::AnnotationHub()
mf.orgdb <- ah[["AH119900"]] #| Taxonomy ID: 9541

llm_interpreted_fm_res <-
  llm_interpret_module(
    object = functional_module_res,
    api_provider = "openai",
    llm_model = "gpt-4o-mini-2024-07-18",
    embedding_model = "text-embedding-3-small",
    api_key = api_key,
    embedding_output_dir = "embedding_output/",
    module_content_number_cutoff = 1,
    orgdb = mf.orgdb
  )

setwd(get_project_wd())
save(llm_interpreted_fm_res, file = "2_data/case_study/01_monkey/05_llm_interpreted_result/02_proteomics_down_llm_interpreted_fm_res.rda")

## Step5: Result visualization ====
load("2_data/case_study/01_monkey/05_llm_interpreted_result/02_proteomics_down_llm_interpreted_fm_res.rda")
setwd("3_data_analysis/05_case_study/01_monkey/")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  level = "functional_module",
  database = c("go", "kegg"),
  degree_cutoff = 1,
  llm_text = TRUE
)

plot

library(Cairo)
CairoPDF("result_visualization/02_proteomics_down_functional_module_network.pdf", width = 12, height = 8)
plot
dev.off()
