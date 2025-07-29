library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source("1_code/2_control_data/05_benchmarking/utils.R")

curated_dt <- readxl::read_xlsx("2_data/control_data.xlsx")

setwd("3_data_analysis/02_control_data/05_benchmarking/")

# aPEAR, clusterprofiler, PAVER, geneFeast
# √ aPEAR needs clusterprofiler enrichment output/ a correctly formatted data frame
## https://gitlab.com/vugene/aPEAR/-/blob/main/vignettes/aPEAR-vignette.Rmd?ref_type=heads
# √ enrichplot needs clusterprofiler enrichment output
## https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#enrichment-map
# √ PAVER: pathway analysis (pathway ID + metrics like p value) results and precomputed pathway embeddings
## browseVignettes("PAVER")
# √ geneFeast:Functional enrichment analysis (FEA) results file(s) Gol file (gene list csv file) + setup YAML file + pathway data frame (type + id + description + geneID)
## https://avigailtaylor.github.io/GeneFEAST/user_guide.html

## Get GO term information
library(org.Hs.eg.db)
go_ids <- curated_dt |> dplyr::filter(database == "GO") |> dplyr::pull(id)
go2egs <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go_ids, columns = "ENTREZID", keytype = "GOALL"))
go2egs_filtered <- go2egs |> dplyr::filter(!(EVIDENCEALL %in% c("IEA", "NAS")))

go_pathway_genes <- go2egs_filtered |>
  group_by(GOALL) |>
  summarise(genes = list(unique(ENTREZID)), .groups = 'drop') %>%
  deframe()

## less gene for GO term (GeneFeast)
library(org.Hs.eg.db)
go_ids <- curated_dt |> dplyr::filter(database == "GO") |> dplyr::pull(id)
go2egs <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go_ids, columns = "ENTREZID", keytype = "GO"))
go2egs_filtered <- go2egs |> dplyr::filter(!(EVIDENCE %in% c("IEA", "NAS")))

go_pathway_genes <- go2egs_filtered |>
  dplyr::group_by(GO) |>
  dplyr::summarise(genes = list(unique(ENTREZID)), .groups = 'drop') |>
  tibble::deframe()

go_genes_df <- data.frame()
for (i in go_ids) {
  if (i %in% names(go_pathway_genes)) {
    genes <- paste(go_pathway_genes[[i]], collapse = "/")
    go_genes_df <- rbind(go_genes_df, data.frame(ID = i,
                                                 GeneID = genes))
  } else {
    go_genes_df <- rbind(go_genes_df, data.frame(ID = i,
                                                 GeneID = NA))
  }
}

genefeast_enrichment_res <- read.csv("genefeast_01_enrichment_result.csv")
genefeast_enrichment_res <- genefeast_enrichment_res |>
  dplyr::mutate(GeneID = dplyr::if_else(
    Type == "GO",
    go_genes_df$GeneID[match(ID, go_genes_df$ID)],
    GeneID
  ))
write.csv(genefeast_enrichment_res, file = "genefeast_01_enrichment_result_short_gene.csv")
## end

go_pathway_info <- get_go_term_info(go_ids = go_ids)

go_pathway_info_term_name <- do.call(rbind, lapply(go_pathway_info, function(x) {
  data.frame(
    ID = x$id,
    Description = x$term_name,
    stringsAsFactors = FALSE
  )
}))

save(go_pathway_info, file = "go_pathway_info.rda")

## Get KEGG pathway information
kegg_ids <- curated_dt |> dplyr::filter(database == "KEGG") |> dplyr::pull(id)

kegg_pathway_info <- get_kegg_pathway_id_info(kegg_ids = kegg_ids)

kegg_pathway_genes <- lapply(kegg_pathway_info, function(pathway) {
  unlist(strsplit(pathway$annotated_genes, "/"))
})
names(kegg_pathway_genes) <- sapply(kegg_pathway_info, function(x) x$id)

kegg_pathway_info_term_name <- do.call(rbind, lapply(kegg_pathway_info, function(x) {
  data.frame(
    ID = x$id,
    Description = x$term_name,
    stringsAsFactors = FALSE
  )
}))

save(kegg_pathway_info, file = "kegg_pathway_info.rda")

## Get Reactome pathway information
reactome_ids <- curated_dt |> dplyr::filter(database == "Reactome") |> dplyr::pull(id)
library(reactome.db)
suppressMessages(reactome2egs <- AnnotationDbi::select(reactome.db, keys = reactome_ids, columns = "ENTREZID", keytype = "PATHID"))

reactome_pathway_genes <- reactome2egs |>
  group_by(PATHID) |>
  summarise(genes = list(unique(ENTREZID)), .groups = 'drop') %>%
  deframe()

reactome_pathway_info <- mapa::get_reactome_pathway_info(reactome_ids = reactome_ids)
reactome_pathway_info_term_name <- do.call(rbind, lapply(reactome_pathway_info, function(x) {
  data.frame(
    ID = x$id,
    Description = x$term_name,
    stringsAsFactors = FALSE
  )
}))

save(reactome_pathway_info, file = "reactome_pathway_info.rda")

## Collect all information
pathway_genes <- c(go_pathway_genes, kegg_pathway_genes, reactome_pathway_genes)
save(pathway_genes, file = "pathway_genes.rda")

term2gene <- data.frame(
  term = rep(names(pathway_genes), lengths(pathway_genes)),
  gene = unlist(pathway_genes)
)

## Use enricher to generate an enrichment object
library(clusterProfiler)
enriched_result <- enricher(gene = unique(term2gene$gene),
                            TERM2GENE = term2gene,
                            pvalueCutoff = 1,  # keep all since already enriched
                            qvalueCutoff = 1,
                            minGSSize = 0,
                            maxGSSize = 2000)

enriched_result@result <- enriched_result@result |>
  dplyr::mutate(
    Description = case_when(
      grepl("^GO", ID) ~ go_pathway_info_term_name$Description[match(ID, go_pathway_info_term_name$ID)],
      grepl("^hsa", ID) ~ kegg_pathway_info_term_name$Description[match(ID, kegg_pathway_info_term_name$ID)],
      grepl("^R-", ID) ~ reactome_pathway_info_term_name$Description[match(ID, reactome_pathway_info_term_name$ID)],
      TRUE ~ Description
    )
  )

save(enriched_result, file = "enriched_result.rda")

## For PAVER ====
### Create a input dataframe using pathway ids
input <- enriched_result@result |> dplyr::select(ID, p.adjust)
term2name <- enriched_result@result |> dplyr::select(ID, Description) |> dplyr::rename(TERM = Description)
text_list <- list()
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

kegg_text <- lapply(kegg_pathway_info,
                    function(x) {
                      paste(
                        x$id,
                        x$term_name,
                        x$term_definition,
                        x$class,
                        sep = "\n"
                      )
                    })

reactome_text <- lapply(reactome_pathway_info,
                        function(x) {
                          paste(
                            x$id,
                            x$term_name,
                            x$term_definition,
                            sep = "\n"
                          )
                        })

all_text_info <- c(go_text, kegg_text, reactome_text)

embedding_list <- lapply(all_text_info, function(x) {
  mapa::get_embedding(chunk = x,
                      api_provider = "openai",
                      api_key = api_key,
                      model_name = "text-embedding-3-small")
})

go_text_order <- lapply(go_pathway_info, function(x) {x$id}) |> unlist()
kegg_text_order <- lapply(kegg_pathway_info, function(x) {x$id}) |> unlist()
reactome_text_order <- lapply(reactome_pathway_info, function(x) {x$id}) |> unlist()
all_text_order <- c(go_text_order, kegg_text_order, reactome_text_order)

names(embedding_list) <- all_text_order

# Create matrix with functional module names as row names
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

save(input, file = "PAVER_01_input.rda")
save(embedding_matrix, file = "PAVER_02_embedding_matrix.rda")
save(term2name, file = "PAVER_03_term2name.rda")

# For GeneFEAST ====
## 1. FEA result
load("enriched_result.rda")
enrichment_result <- enriched_result@result |>
  dplyr::select(ID, Description, geneID) |>
  dplyr::rename(GeneID = geneID) |>
  dplyr::mutate(
    Type = case_when(
      grepl("^GO", ID) ~ "GO",
      grepl("^hsa", ID) ~ "KEGG",
      grepl("^R-", ID) ~ "Reactome",
      TRUE ~ Description
    )
  )

export(enrichment_result, file = "genefeast_01_enrichment_result.csv")

## 2. Gene of interest (GoI) file
goI <- data.frame(GeneID = unique(unlist(pathway_genes)),
                  Value = 1)

export(goI, file = "genefeast_02_gol.csv")

## 3. A YAML setup file
