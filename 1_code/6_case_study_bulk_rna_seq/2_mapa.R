library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Cluster U ====
## Step1: data input and preprocessing ====
load("2_data/case_study/bulk_rna_seq/dt_u.rda")

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
     file = "3_data_analysis/06_case_study_bulk_rna_seq/enrich_pathway_res_cluster_u.rda")

## Step3: Similarity calculation and clustering ====
load("3_data_analysis/06_case_study_bulk_rna_seq/enrich_pathway_res_cluster_u.rda")
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

save(merged_pathways, file = "3_data_analysis/06_case_study_bulk_rna_seq/merged_pathways_res_cluster_u.rda")

merged_modules <- merge_modules(
  object = merged_pathways,
  sim.cutoff = 0.3,
  measure_method = "jaccard",
  cluster_method = "girvan newman"
)

# merged_modules <- merge_modules(
#   object = merged_pathways,
#   sim.cutoff = 0.65,
#   measure_method = "jaccard",
#   cluster_method = "binary cut"
# )

# merged_modules <- merge_modules(
#   object = merged_pathways,
#   sim.cutoff = 0.95,
#   measure_method = "jaccard",
#   cluster_method = "hierarchical",
#   hclust.method = "complete"
# )

save(merged_modules, file = "3_data_analysis/06_case_study_bulk_rna_seq/merged_modules_res_cluster_u.rda")

### 3.2 Biotext embedding method ====
embedsim_res <-
  get_bioembedsim(
    object = enrich_pathway_res,
    api_provider = "openai",
    api_key = api_key,
    text_embedding_model = "text-embedding-3-small",
    database = c("go", "kegg")
  )

save(embedsim_res, file = "3_data_analysis/06_case_study_bulk_rna_seq/embedsim_res_cluster_u.rda")

load("3_data_analysis/06_case_study_bulk_rna_seq/embedsim_res_cluster_u.rda")

functional_module_res <-
  merge_pathways_bioembedsim(
    object = embedsim_res,
    sim.cutoff = 0.6,
    cluster_method = "binary cut"
  )

save(functional_module_res, file = "3_data_analysis/06_case_study_bulk_rna_seq/functional_module_res_cluster_u.rda")

## Step4: LLM interpretation ====
load("3_data_analysis/06_case_study_bulk_rna_seq/functional_module_res_cluster_u.rda")

llm_interpreted_res <-
  llm_interpret_module(
    object = functional_module_res,
    module_content_number_cutoff = 1,
    llm_model = "gpt-4o-mini-2024-07-18",
    embedding_model = "text-embedding-3-small",
    api_key = api_key,
    embedding_output_dir = "3_data_analysis/06_case_study_bulk_rna_seq/embedding_output/",
    phenotype = "Aging",
    orgdb = variable_info$orgdb,
    output_prompt = TRUE
  )

save(llm_interpreted_res, file = "3_data_analysis/06_case_study_bulk_rna_seq/llm_interpreted_res_biotext_cluster_u.rda")

llm_interpreted_overlap_res <-
  llm_interpret_module(
    object = merged_modules,
    module_content_number_cutoff = 1,
    llm_model = "gpt-4o-mini-2024-07-18",
    embedding_model = "text-embedding-3-small",
    api_key = api_key,
    embedding_output_dir = "3_data_analysis/06_case_study_bulk_rna_seq/embedding_output/",
    phenotype = "Aging",
    orgdb = variable_info$orgdb,
    output_prompt = TRUE
  )

save(llm_interpreted_overlap_res, file = "3_data_analysis/06_case_study_bulk_rna_seq/llm_interpreted_res_overlap_cluster_u.rda")

## Step5: Result visualization ====

