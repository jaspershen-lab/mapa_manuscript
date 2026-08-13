## Sensitivity analysis on the informative GO set benchmark (k = 20)
## Step 1: Build an enriched pathway object from
## 2_data/informative_go_set/k_20/mapa_input.csv (313 informative GO BP terms,
## treated as enriched pathways) and calculate the biotext embedding
## similarity matrix with mapa::get_bioembedsim().

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# OpenAI API key from ~/.Renviron
api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "") {
  stop("OPENAI_API_KEY not found in ~/.Renviron")
}

# 1. Read the k = 20 benchmark GO terms ====
mapa_input <- readr::read_csv("2_data/informative_go_set/k_20/mapa_input.csv",
                              show_col_types = FALSE)

# 2. Build the enriched pathway object ====
## Use the example output object as a template and replace the GO enrichment
## result with the benchmark terms (same strategy as
## 1_code/2_control_data/05_benchmarking/02_run_analysis.R)
load("3_data_analysis/07_example_output/openai_semantic_sim_matrix.rda")
enriched_res_go_k20 <- openai_semantic_sim_matrix$enriched_pathway

go_result <-
  data.frame(
    ONTOLOGY = "BP",
    ID = mapa_input$go_id,
    Description = mapa_input$name,
    GeneRatio = paste0(mapa_input$gene_set_size, "/17592"),
    BgRatio = paste0(mapa_input$gene_set_size, "/17592"),
    RichFactor = 1,
    FoldEnrichment = 1,
    zScore = 1,
    pvalue = 0.001,
    p_adjust = 0.01,
    qvalue = 0.01,
    geneID = "",
    Count = mapa_input$gene_set_size
  )
rownames(go_result) <- go_result$ID

enriched_res_go_k20@enrichment_go_result@result <- go_result

save(enriched_res_go_k20,
     file = "3_data_analysis/10_sensitivity_res_GO/enriched_res_go_k20.rda")

# 3. Biotext embedding + cosine similarity via mapa ====
bioembed_sim_res <-
  get_bioembedsim(
    object = enriched_res_go_k20,
    embedding_source = "realtime",
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = "go",
    p.adjust.cutoff.go = 0.05,
    count.cutoff.go = 0
  )

save(bioembed_sim_res,
     file = "3_data_analysis/10_sensitivity_res_GO/bioembed_sim_res.rda")

# 4. Extract similarity matrix / edge list for the clustering scans ====
embedding_sim_matrix <- bioembed_sim_res$sim_matrix

## Sanity check: all benchmark terms should be present
missing_ids <- setdiff(mapa_input$go_id, rownames(embedding_sim_matrix))
if (length(missing_ids) > 0) {
  warning("Terms missing from similarity matrix: ",
          paste(missing_ids, collapse = ", "))
}
message("Similarity matrix: ", nrow(embedding_sim_matrix), " x ",
        ncol(embedding_sim_matrix))

embedding_sim_df <-
  as.data.frame.table(embedding_sim_matrix, responseName = "sim") |>
  dplyr::filter(Var1 != Var2) |>
  dplyr::rename(from = Var1, to = Var2) |>
  dplyr::mutate(across(c(from, to), as.character)) |>
  dplyr::filter(from < to)

save(embedding_sim_matrix,
     file = "3_data_analysis/10_sensitivity_res_GO/embedding_sim_matrix.rda")
save(embedding_sim_df,
     file = "3_data_analysis/10_sensitivity_res_GO/embedding_sim_df.rda")
