## Benchmarking on the informative GO set benchmark (k = 20)
## Step 1: Prepare inputs for enrichplot, aPEAR and PAVER.
## Adapted from 1_code/2_control_data/05_benchmarking/01_data_preparation.R
## (313 informative GO BP terms; gene sets from module_gene_sets.tsv,
## deduplicated UniProtKB accessions)

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source("1_code/2_control_data/05_benchmarking/utils.R")

api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "") {
  stop("OPENAI_API_KEY not found in ~/.Renviron")
}

mapa_input <- readr::read_csv("2_data/informative_go_set/k_20/mapa_input.csv",
                              show_col_types = FALSE)
gene_sets <- readr::read_tsv("2_data/informative_go_set/k_20/module_gene_sets.tsv",
                             show_col_types = FALSE)

## Ground truth (57 modules, 313 member terms)
ground_truth_dt <- readr::read_csv("2_data/informative_go_set/k_20/ground_truth.csv",
                                   show_col_types = FALSE)
ground_truth_dt <-
  ground_truth_dt |>
  dplyr::mutate(ground_truth_label = as.numeric(factor(module_id))) |>
  dplyr::rename(id = go_id) |>
  dplyr::select(id, ground_truth_label)
save(ground_truth_dt,
     file = "3_data_analysis/11_informative_go_benchmarking/ground_truth_dt.rda")

setwd("3_data_analysis/11_informative_go_benchmarking/")

# 1. Pathway gene sets ====
pathway_genes <- strsplit(gene_sets$uniprot_accessions, ";")
names(pathway_genes) <- gene_sets$go_id
save(pathway_genes, file = "pathway_genes.rda")

term2gene <- data.frame(
  term = rep(names(pathway_genes), lengths(pathway_genes)),
  gene = unlist(pathway_genes)
)

# 2. Use enricher to generate an enrichment object ====
library(clusterProfiler)
enriched_result <- enricher(gene = unique(term2gene$gene),
                            TERM2GENE = term2gene,
                            pvalueCutoff = 1,  # keep all since already enriched
                            qvalueCutoff = 1,
                            minGSSize = 0,
                            maxGSSize = 2000)

enriched_result@result <- enriched_result@result |>
  dplyr::mutate(
    Description = mapa_input$name[match(ID, mapa_input$go_id)]
  )

stopifnot(nrow(enriched_result@result) == nrow(mapa_input))
save(enriched_result, file = "enriched_result.rda")

# 3. For PAVER ====
### Input dataframe and term2name
input <- enriched_result@result |> dplyr::select(ID, p.adjust)
term2name <- enriched_result@result |>
  dplyr::select(ID, Description) |>
  dplyr::rename(TERM = Description)

### Pathway text embeddings (same text format as the control-data preparation)
go_pathway_info <- get_go_term_info(go_ids = mapa_input$go_id)
save(go_pathway_info, file = "go_pathway_info.rda")

go_text <- lapply(go_pathway_info,
                  function(x) {
                    paste(
                      x$id,
                      x$sub_ontology,
                      x$term_name,
                      x$term_definition,
                      sep = "\n"
                    )
                  })

embedding_list <- lapply(go_text, function(x) {
  mapa::get_embedding(chunk = x,
                      api_provider = "openai",
                      api_key = api_key,
                      model_name = "text-embedding-3-small")
})

go_text_order <- lapply(go_pathway_info, function(x) {x$id}) |> unlist()
names(embedding_list) <- go_text_order

n_dims <- length(embedding_list[[1]])
embedding_matrix <- matrix(NA,
                           nrow = length(embedding_list),
                           ncol = n_dims,
                           dimnames = list(names(embedding_list),
                                           paste0("dim_", 1:n_dims)))
for (i in 1:length(embedding_list)) {
  pathway_id <- names(embedding_list)[i]
  embedding_matrix[i, ] <- embedding_list[[pathway_id]]
}

stopifnot(!anyNA(embedding_matrix))

save(input, file = "PAVER_01_input.rda")
save(embedding_matrix, file = "PAVER_02_embedding_matrix.rda")
save(term2name, file = "PAVER_03_term2name.rda")
