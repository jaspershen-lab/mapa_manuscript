library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
library(org.Hs.eg.db)
source("1_code/100_tools.R")
library(readxl)

# Age associated: p < 0.05, |beta| > 0.008
## Thymus ====
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/thymus/up", recursive = TRUE)
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/thymus/thymus_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "taxid_9606/thymus/up_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/thymus/thymus_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways
save(up_prot_enriched_pathways, file = "taxid_9606/thymus/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/thymus/thymus_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

sum(annotated_met_db$kegg_id == "")
sum(is.na(annotated_met_db$chebi_id))

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::select(-name) |>
  dplyr::rename(keggid = kegg_id, cpd_name = ID) |>
  dplyr::select(keggid, cpd_name, everything()) |>
  dplyr::distinct(keggid, .keep_all = TRUE) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
met_enriched_pathways

save(met_enriched_pathways, file = "taxid_9606/thymus/met_enriched_pathways.rda")

## Thymus (down) ====
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/thymus/down")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/thymus/thymus_age_related_mrna.xlsx")
summary(rna_raw_dt)
down_rna <- rna_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_rna_variable_info <- convert_id(
  data = down_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "taxid_9606/thymus/down/down_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/thymus/thymus_age_related_protein.xlsx")
summary(prot_raw_dt)
down_prot <- prot_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_prot_variable_info <- convert_id(
  data = down_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways
save(down_prot_enriched_pathways, file = "taxid_9606/thymus/down/down_prot_enriched_pathways.rda")

## Spleen ====
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/spleen/up", recursive = TRUE)
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/spleen/spleen_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

rio::export(up_rna, file = "taxid_9606/spleen/up/demo_up_T_list.xlsx")

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "taxid_9606/spleen/up/up_rna_enriched_pathways.rda")

# ### add reactome
# up_rna_enriched_pathways <-
#   enrich_pathway(
#     variable_info = up_rna_variable_info,
#     query_type = "gene",
#     database = c("go", "kegg", "reactome"),
#     go.orgdb = org.Hs.eg.db,
#     go.keytype = "SYMBOL",
#     go.ont = "ALL",
#     kegg.organism = "hsa",
#     kegg.keytype = "kegg",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH"
#   )
# up_rna_enriched_pathways


# protein
prot_raw_dt <- read_xlsx("raw_data/spleen/spleen_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

rio::export(up_prot, file = "taxid_9606/spleen/up/demo_up_P_list.xlsx")

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "taxid_9606/spleen/up/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/spleen/spleen_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::select(-name) |>
  dplyr::rename(keggid = kegg_id, cpd_name = ID) |>
  dplyr::select(keggid, cpd_name, everything()) |>
  dplyr::distinct(keggid, .keep_all = TRUE) |>
  dplyr::filter(keggid != "") -> met_dt

rio::export(met_dt, file = "taxid_9606/spleen/demo_M_list.xlsx")

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "taxid_9606/spleen/met_enriched_pathways.rda")

## Spleen (down) ====
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/spleen/down")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/spleen/spleen_age_related_mrna.xlsx")
summary(rna_raw_dt)
down_rna <- rna_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_rna_variable_info <- convert_id(
  data = down_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "taxid_9606/spleen/down/down_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/spleen/spleen_age_related_protein.xlsx")
summary(prot_raw_dt)
down_prot <- prot_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_prot_variable_info <- convert_id(
  data = down_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways

save(down_prot_enriched_pathways, file = "taxid_9606/spleen/down/down_prot_enriched_pathways.rda")

## Ovary ====
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/ovary/up", recursive = TRUE)
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/ovary/ovary_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "taxid_9606/ovary/up/up_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/ovary/ovary_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "taxid_9606/ovary/up/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/ovary/ovary_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::select(-name) |>
  dplyr::rename(keggid = kegg_id, cpd_name = ID) |>
  dplyr::select(keggid, cpd_name, everything()) |>
  dplyr::distinct(keggid, .keep_all = TRUE) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "taxid_9606/ovary/met_enriched_pathways.rda")

## Ovary (down) ====
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/ovary/down")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/ovary/ovary_age_related_mrna.xlsx")
summary(rna_raw_dt)
down_rna <- rna_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_rna_variable_info <- convert_id(
  data = down_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "taxid_9606/ovary/down/down_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/ovary/ovary_age_related_protein.xlsx")
summary(prot_raw_dt)
down_prot <- prot_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_prot_variable_info <- convert_id(
  data = down_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways

save(down_prot_enriched_pathways, file = "taxid_9606/ovary/down/down_prot_enriched_pathways.rda")

## Stomach ====
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/stomach/up", recursive = TRUE)
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/stomach/stomach_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "taxid_9606/stomach/up/up_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/stomach/stomach_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "taxid_9606/stomach/up/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/stomach/stomach_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::select(-name) |>
  dplyr::rename(keggid = kegg_id, cpd_name = ID) |>
  dplyr::select(keggid, cpd_name, everything()) |>
  dplyr::distinct(keggid, .keep_all = TRUE) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "taxid_9606/stomach/met_enriched_pathways.rda")

## Stomach (down) ====
setwd(r4projects::get_project_wd())
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/stomach/down")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/stomach/stomach_age_related_mrna.xlsx")
summary(rna_raw_dt)
down_rna <- rna_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_rna_variable_info <- convert_id(
  data = down_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "taxid_9606/stomach/down/down_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/stomach/stomach_age_related_protein.xlsx")
summary(prot_raw_dt)
down_prot <- prot_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_prot_variable_info <- convert_id(
  data = down_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways

save(down_prot_enriched_pathways, file = "taxid_9606/stomach/down/down_prot_enriched_pathways.rda")

## Kidney ====
setwd(r4projects::get_project_wd())
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/kidney/up", recursive = TRUE)
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rm(list = ls())
rna_raw_dt <- read_xlsx("raw_data/kidney/kidney_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "taxid_9606/kidney/up/up_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/kidney/kidney_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "taxid_9606/kidney/up/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/kidney/kidney_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::select(-name) |>
  dplyr::rename(keggid = kegg_id, cpd_name = ID) |>
  dplyr::select(keggid, cpd_name, everything()) |>
  dplyr::distinct(keggid, .keep_all = TRUE) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "taxid_9606/kidney/met_enriched_pathways.rda")

## Kidney (down) ====
setwd(r4projects::get_project_wd())
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/kidney/down")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/kidney/kidney_age_related_mrna.xlsx")
summary(rna_raw_dt)
down_rna <- rna_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_rna_variable_info <- convert_id(
  data = down_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "taxid_9606/kidney/down/down_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/kidney/kidney_age_related_protein.xlsx")
summary(prot_raw_dt)
down_prot <- prot_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_prot_variable_info <- convert_id(
  data = down_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways

save(down_prot_enriched_pathways, file = "taxid_9606/kidney/down/down_prot_enriched_pathways.rda")

## Aortic arch ====
setwd(r4projects::get_project_wd())
rm(list = ls())
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/aortic_arch/up", recursive = TRUE)
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/aortic_arch/aortic_arch_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "taxid_9606/aortic_arch/up/up_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/aortic_arch/aortic_arch_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "taxid_9606/aortic_arch/up/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/aortic_arch/aortic_arch_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::select(-name) |>
  dplyr::rename(keggid = kegg_id, cpd_name = ID) |>
  dplyr::select(keggid, cpd_name, everything()) |>
  dplyr::distinct(keggid, .keep_all = TRUE) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "taxid_9606/aortic_arch/met_enriched_pathways.rda")

## Aortic arch (down) ====
setwd(r4projects::get_project_wd())
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/aortic_arch/down")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/aortic_arch/aortic_arch_age_related_mrna.xlsx")
summary(rna_raw_dt)
down_rna <- rna_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_rna_variable_info <- convert_id(
  data = down_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "taxid_9606/aortic_arch/down/down_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/aortic_arch/aortic_arch_age_related_protein.xlsx")
summary(prot_raw_dt)
down_prot <- prot_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_prot_variable_info <- convert_id(
  data = down_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways

save(down_prot_enriched_pathways, file = "taxid_9606/aortic_arch/down/down_prot_enriched_pathways.rda")

## Thyroid ====
setwd(r4projects::get_project_wd())
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/thyroid/up", recursive = TRUE)
setwd("2_data/case_study/02_monkey_multi_omics")
rm(list = ls())
# RNA
rna_raw_dt <- read_xlsx("raw_data/thyroid/thyroid_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "taxid_9606/thyroid/up/up_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/thyroid/thyroid_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "taxid_9606/thyroid/up/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/thyroid/thyroid_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::select(-name) |>
  dplyr::rename(keggid = kegg_id, cpd_name = ID) |>
  dplyr::select(keggid, cpd_name, everything()) |>
  dplyr::distinct(keggid, .keep_all = TRUE) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "taxid_9606/thyroid/met_enriched_pathways.rda")

## Thyroid (down) ====
setwd(r4projects::get_project_wd())
# dir.create("2_data/case_study/02_monkey_multi_omics/taxid_9606/thyroid/down")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/thyroid/thyroid_age_related_mrna.xlsx")
summary(rna_raw_dt)
down_rna <- rna_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_rna_variable_info <- convert_id(
  data = down_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "taxid_9606/thyroid/down/down_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/thyroid/thyroid_age_related_protein.xlsx")
summary(prot_raw_dt)
down_prot <- prot_raw_dt |>
  dplyr::filter(beta < -0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

down_prot_variable_info <- convert_id(
  data = down_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Hs.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Hs.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "hsa",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways

save(down_prot_enriched_pathways, file = "taxid_9606/thyroid/down/down_prot_enriched_pathways.rda")

