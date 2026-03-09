library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(readxl)
library(org.Mmu.eg.db)

# Age associated: p < 0.05, |beta| > 0.008
## Thymus ====
dir.create("2_data/case_study/02_monkey_multi_omics/thymus")
# RNA
rna_raw_dt <- read_xlsx("2_data/case_study/02_monkey_multi_omics/raw_data/thymus/thymus_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  filter(beta > 0.008 & Pvalue < 0.05) |>
  rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Mmu.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

save(up_rna_enriched_pathways, file = "2_data/case_study/02_monkey_multi_omics/thymus/up_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("2_data/case_study/02_monkey_multi_omics/raw_data/thymus/thymus_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  filter(beta > 0.008 & Pvalue < 0.05) |>
  rename(symbol = ID)

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Mmu.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

save(up_prot_enriched_pathways, file = "2_data/case_study/02_monkey_multi_omics/thymus/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("2_data/case_study/02_monkey_multi_omics/raw_data/thymus/thymus_age_related_metabolites.xlsx")
load("2_data/case_study/02_monkey_female/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

sum(annotated_met_db$kegg_id == "")
sum(is.na(annotated_met_db$chebi_id))

met_raw_dt |>
  filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  left_join(annotated_met_db, by = "ID") |>
  rename(keggid = kegg_id) |>
  select(keggid, everything()) |>
  filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "mcc",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

save(met_enriched_pathways, file = "2_data/case_study/02_monkey_multi_omics/thymus/met_enriched_pathways.rda")


## Spleen ====
# dir.create("2_data/case_study/02_monkey_multi_omics/spleen")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/spleen/spleen_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Mmu.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "spleen/up_rna_enriched_pathways.rda")

# protein
prot_raw_dt <- read_xlsx("raw_data/spleen/spleen_age_related_protein.xlsx")
summary(prot_raw_dt)
up_prot <- prot_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_prot_variable_info <- convert_id(
  data = up_prot,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Mmu.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "spleen/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/spleen/spleen_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::rename(keggid = kegg_id) |>
  dplyr::select(keggid, everything()) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "mcc",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "spleen/met_enriched_pathways.rda")


## Ovary ====
# dir.create("2_data/case_study/02_monkey_multi_omics/ovary")
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
  organism = org.Mmu.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "ovary/up_rna_enriched_pathways.rda")

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
  organism = org.Mmu.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "ovary/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/ovary/ovary_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::rename(keggid = kegg_id) |>
  dplyr::select(keggid, everything()) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "mcc",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "ovary/met_enriched_pathways.rda")

## Stomach ====
# dir.create("2_data/case_study/02_monkey_multi_omics/stomach")
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
  organism = org.Mmu.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "stomach/up_rna_enriched_pathways.rda")

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
  organism = org.Mmu.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "stomach/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/stomach/stomach_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::rename(keggid = kegg_id) |>
  dplyr::select(keggid, everything()) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "mcc",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "stomach/met_enriched_pathways.rda")

## Kidney ====
# dir.create("2_data/case_study/02_monkey_multi_omics/kidney")
setwd("2_data/case_study/02_monkey_multi_omics")
# RNA
rna_raw_dt <- read_xlsx("raw_data/kidney/kidney_age_related_mrna.xlsx")
summary(rna_raw_dt)
up_rna <- rna_raw_dt |>
  dplyr::filter(beta > 0.008 & Pvalue < 0.05) |>
  dplyr::rename(symbol = ID)

up_rna_variable_info <- convert_id(
  data = up_rna,
  query_type = "gene",
  from_id_type = "symbol",
  organism = org.Mmu.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "kidney/up_rna_enriched_pathways.rda")

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
  organism = org.Mmu.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "kidney/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/kidney/kidney_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::rename(keggid = kegg_id) |>
  dplyr::select(keggid, everything()) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "mcc",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "kidney/met_enriched_pathways.rda")

# Aortic arch ====
# dir.create("2_data/case_study/02_monkey_multi_omics/aortic_arch")
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
  organism = org.Mmu.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "aortic_arch/up_rna_enriched_pathways.rda")

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
  organism = org.Mmu.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "aortic_arch/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/aortic_arch/aortic_arch_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::rename(keggid = kegg_id) |>
  dplyr::select(keggid, everything()) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "mcc",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "aortic_arch/met_enriched_pathways.rda")

# Aortic arch (down) ====
# dir.create("2_data/case_study/02_monkey_multi_omics/aortic_arch_down")
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
  organism = org.Mmu.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "aortic_arch_down/down_rna_enriched_pathways.rda")

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
  organism = org.Mmu.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways

save(down_prot_enriched_pathways, file = "aortic_arch_down/down_prot_enriched_pathways.rda")

# Thyroid ====
# dir.create("2_data/case_study/02_monkey_multi_omics/thyroid")
setwd("2_data/case_study/02_monkey_multi_omics")
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
  organism = org.Mmu.eg.db
)

up_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = up_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_rna_enriched_pathways

save(up_rna_enriched_pathways, file = "thyroid/up_rna_enriched_pathways.rda")

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
  organism = org.Mmu.eg.db
)

up_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = up_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

up_prot_enriched_pathways

save(up_prot_enriched_pathways, file = "thyroid/up_prot_enriched_pathways.rda")

# metabolite
met_raw_dt <- read_xlsx("raw_data/thyroid/thyroid_age_related_metabolites.xlsx")
load("raw_data/metabolite_annotation.Rdata")
annotated_met_db <- met.header.all.hmdb |>
  rownames_to_column(var = "ID")

met_raw_dt |>
  dplyr::filter(abs(beta) > 0.008 & Pvalue < 0.05) |>
  dplyr::left_join(annotated_met_db, by = "ID") |>
  dplyr::rename(keggid = kegg_id) |>
  dplyr::select(keggid, everything()) |>
  dplyr::filter(keggid != "") -> met_dt

met_enriched_pathways <-
  enrich_pathway(
    variable_info = met_dt,
    query_type = "metabolite",
    database = c("metkegg"),
    met_organism = "mcc",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

met_enriched_pathways

save(met_enriched_pathways, file = "thyroid/met_enriched_pathways.rda")

# Thyroid (down) ====
# dir.create("2_data/case_study/02_monkey_multi_omics/thyroid_down")
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
  organism = org.Mmu.eg.db
)

down_rna_enriched_pathways <-
  enrich_pathway(
    variable_info = down_rna_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_rna_enriched_pathways

save(down_rna_enriched_pathways, file = "thyroid_down/down_rna_enriched_pathways.rda")

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
  organism = org.Mmu.eg.db
)

down_prot_enriched_pathways <-
  enrich_pathway(
    variable_info = down_prot_variable_info,
    query_type = "gene",
    database = c("go", "kegg"),
    go.orgdb = org.Mmu.eg.db,
    go.keytype = "SYMBOL",
    go.ont = "ALL",
    kegg.organism = "mcc",
    kegg.keytype = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )

down_prot_enriched_pathways

save(down_prot_enriched_pathways, file = "thyroid_down/down_prot_enriched_pathways.rda")

