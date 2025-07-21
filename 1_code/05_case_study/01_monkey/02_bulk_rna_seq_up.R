library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# load("2_data/case_study/01_monkey/01_bulk_rna_seq/dt_u.rda")

setwd("3_data_analysis/05_case_study/01_monkey/")

# Step1: data input and preprocessing ====
### Create an AnnotationHub object to access the database
ah <- AnnotationHub::AnnotationHub()
### Query the object for Macaca fascicularis
mf_query_result <- AnnotationHub::query(ah, c("OrgDb", "9541"))
mf_query_result$ah_id
### Get the orgDb object
mf.orgdb <- ah[["AH119900"]] #| Taxonomy ID: 9541

variable_info <- convert_id(data = dt_u,
                            query_type = "gene",
                            from_id_type = "ensembl",
                            ah_id = "AH119900",
                            return_orgdb = TRUE)
variable_info$data <- variable_info$data |> dplyr::mutate(symbol = Symbol) |> dplyr::select(-Symbol)
variable_info$data$symbol[variable_info$data$symbol == "NA"] <- NA
save(variable_info, file = "01_bulk_rna_seq_up_variable_info.rda")

## Step2: enrichment analysis ====
## cannot save the variable info since the orgdb can not use after save in .rda file
enrich_pathway_res <- enrich_pathway(
  variable_info = variable_info$data,
  query_type = "gene",
  database = c("go", "kegg"),
  save_to_local = FALSE,
  go.orgdb = variable_info$orgdb,
  go.keytype = "ENSEMBL",
  go.ont = "ALL",
  kegg.organism = "mcf", # Taxid:9541
  kegg.keytype = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

save(enrich_pathway_res,
     file = "01_bulk_rna_seq_up_enrich_pathway_res_cluster_u.rda")

## Step3: Similarity calculation and clustering ====
# load("01_bulk_rna_seq_up_enrich_pathway_res_cluster_u.rda")
### 3.1 Traditional similarity method ====
merged_pathways <- merge_pathways(
  object = enrich_pathway_res,
  database = c("go", "kegg"),
  sim.cutoff.go = 0.5,
  sim.cutoff.kegg = 0.5,
  measure.method.go = "Sim_Lin_1998",
  go.orgdb = variable_info$orgdb,
  measure.method.kegg = "jaccard"
)

save(merged_pathways, file = "01_bulk_rna_seq_up_merged_pathways.rda")

load("01_bulk_rna_seq_up_merged_pathways.rda")
optimal_clustering_test <- determine_optimal_clusters(
  object = merged_pathways,
  cutoff_range = c(0.1, 0.9),
  cutoff_increment = 0.05,
  methods = c("h_ward.D", "h_ward.D2", "h_complete",
              "binary_cut", "louvain", "walktrap", "edge_betweenness",
              "fast_greedy")
)

save(optimal_clustering_test, file = "01_bulk_rna_seq_up_optimal_clustering_test.rda")
optimal_clustering_test$evaluation_plot
optimal_clustering_test$best_combination # use louvain, cutoff: 0.45

functional_module_res <- get_functional_modules(
  object = merged_pathways,
  sim.cutoff = 0.45,
  cluster_method = "louvain"
)

save(functional_module_res, file = "01_bulk_rna_seq_up_functional_module_res.rda")

### 3.2 Biotext embedding method ====
# load("01_bulk_rna_seq_up_enrich_pathway_res_cluster_u.rda")
embedsim_res <-
  get_bioembedsim(
    object = enrich_pathway_res,
    api_provider = "siliconflow",
    api_key = api_key,
    text_embedding_model = "Qwen/Qwen3-Embedding-8B",
    database = c("go", "kegg")
  )

save(embedsim_res, file = "01_bulk_rna_seq_up_embedsim_res.rda")

optimal_clustering_test <- determine_optimal_clusters(
  object = embedsim_res,
  cutoff_range = c(0.1, 0.9),
  cutoff_increment = 0.05,
  methods = c("h_ward.D", "h_ward.D2", "h_complete",
              "binary_cut", "louvain", "walktrap", "edge_betweenness",
              "fast_greedy")
)

save(optimal_clustering_test, file = "01_bulk_rna_seq_up_optimal_clustering_test_biotext.rda")
optimal_clustering_test$evaluation_plot
optimal_clustering_test$best_combination # use louvain, cutoff: 0.45

#
load("01_bulk_rna_seq_up_embedsim_res.rda")

functional_module_res <- get_functional_modules(
  object = embedsim_res,
  sim.cutoff = 0.55,
  cluster_method = "h_ward.D2"
)

functional_module_res <- get_functional_modules(
  object = embedsim_res,
  sim.cutoff = 0.8,
  cluster_method = "louvain"
)

clustering_quality <- assess_clustering_quality(functional_module_res)
clustering_quality$size_plot
clustering_quality$evaluation_plot

save(functional_module_res, file = "01_bulk_rna_seq_up_functional_module_res_biotext.rda")

## Step4: LLM interpretation ====
ah <- AnnotationHub::AnnotationHub()
# ### Query the object for Macaca fascicularis
# mf_query_result <- AnnotationHub::query(ah, c("OrgDb", "9541"))
# mf_query_result$ah_id
# ### Get the orgDb object
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

save(llm_interpreted_fm_res, file = "2_data/case_study/01_monkey/05_llm_interpreted_result/01_bulk_rna_seq_up_llm_interpreted_fm_res.rda")

## Step5: Result visualization ====
load("2_data/case_study/01_monkey/05_llm_interpreted_result/01_bulk_rna_seq_up_llm_interpreted_fm_res.rda")
setwd("3_data_analysis/05_case_study/01_monkey/")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  level = "functional_module",
  database = c("go", "kegg"),
  degree_cutoff = 2,
  llm_text = TRUE
)

plot

library(Cairo)
CairoPDF("result_visualization/01_bulk_rna_seq_up_functional_module_network.pdf", width = 12, height = 8)
plot
dev.off()
