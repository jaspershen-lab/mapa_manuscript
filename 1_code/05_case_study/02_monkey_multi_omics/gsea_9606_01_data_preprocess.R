library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
library(org.Hs.eg.db)
source("1_code/100_tools.R")
library(readxl)

# dir.create("2_data/case_study/02_monkey_multi_omics/gsea_taxid_9606")

{
  variable_info
  query_type = "gene"
  order_by = "fc"
  database = c("go", "kegg", "reactome")
  save_to_local = FALSE
  path = "result"
  # GO parameters
  go.orgdb = org.Hs.eg.db
  go.keytype = "SYMBOL"
  go.ont = "ALL"
  # KEGG parameters
  kegg.organism = "hsa"
  kegg.keytype = "kegg"
  reactome.organism = "human"
  # Common parameters
  exponent = 1
  eps = 1e-10
  verbose = TRUE
  seed = FALSE
  pvalueCutoff = 0.05
  pAdjustMethod = "BH"
  qvalueCutoff = 0.2
  minGSSize = 10
  maxGSSize = 500
  readable = FALSE
}

# Age associated: p < 0.05, |beta| > 0.008
## Thymus ====
# dir.create("2_data/case_study/02_monkey_multi_omics/gsea_taxid_9606/thymus/up", recursive = TRUE)
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/thymus/thymus_age_related_mrna.xlsx")
summary(rna_raw_dt)
rna <- rna_raw_dt |>
  dplyr::filter(Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

rna_variable_info <- convert_id(
  data = rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

sum(is.na(rna_variable_info$entrezid))
sum(duplicated(rna_variable_info$entrezid))

rna_variable_info <- rna_variable_info |> dplyr::rename(p_value_adjust = Pvalue)

gsea_rna_enriched_pathways <-
  do_gsea(
    variable_info = rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    order_by = "beta",
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "taxid_9606/thymus/up_rna_enriched_pathways.rda")
