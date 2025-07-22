library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

# 01_bulk_rna_seq ====
dt <- import("2_data/case_study/01_monkey/mmc2.xlsx",
             skip = 1,
             header = TRUE,
             sheet = 2)

setwd("2_data/case_study/01_monkey/")

dt <- as.data.table(dt)
tissue_specific_info <- dt[Cluster %in% c("Cluster U", "Cluster D"), .N, by = .(Tissue, Cluster)]

dt <- dt |> rename(ensembl = GeneID)

# cluster U
dt_u <- dt |> filter(Cluster == "Cluster U")
dt_u <- dt_u |> distinct(ensembl, .keep_all = TRUE)

# cluster D
dt_d <- dt |> filter(Cluster == "Cluster D")
dt_d <- dt_d |> distinct(ensembl, .keep_all = TRUE)

save(dt_u, file = "01_bulk_rna_seq/dt_u.rda")
save(dt_d, file = "01_bulk_rna_seq/dt_d.rda")

# 02_proteomics_data ====
p_dt <- import("2_data/case_study/01_monkey/mmc2.xlsx",
               skip = 1,
               header = TRUE,
               sheet = 6)

setwd("2_data/case_study/01_monkey/")
# convert gene variant id to common gene id
p_dt_cleaned <- p_dt |>
  mutate(ID = gsub("\\.\\d+$", "", ID))

p_dt_cleaned_up <- p_dt_cleaned |> dplyr::filter(Cluster == "Cluster U")
p_dt_cleaned_down <- p_dt_cleaned |> dplyr::filter(Cluster == "Cluster D")

save(p_dt_cleaned_up, file = "02_proteomics/p_dt_cleaned_up.rda")
save(p_dt_cleaned_down, file = "02_proteomics/p_dt_cleaned_down.rda")

# uniprot_pattern <- "^[A-NR-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2}$|^[OPQ][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9]$"
#
# uniprot_ids <- p_dt_cleaned$ID[grepl(uniprot_pattern, p_dt_cleaned$ID)]
# gene_names <-
#   p_dt_cleaned |>
#   dplyr::filter(!grepl(uniprot_pattern, ID)) |>
#   pull(ID)
# library(biomaRt)
# ensembl_mart <- useMart(biomart = "ensembl")
# mfa_mart <- useDataset(dataset = "mfascicularis_gene_ensembl", mart = ensembl_mart)
#
# converted_id_gn_2_ensesmbl <-
#   getBM(
#     attributes = c("ensembl_gene_id", "external_gene_name"),
#     filters    = "external_gene_name",
#     values = gene_names,
#     mart = mfa_mart
#   )
#
# converted_id_uniprot_2_ensembl <-
#   getBM(
#     attributes = c("ensembl_gene_id", "uniprotsptrembl"),
#     filters    = "uniprotsptrembl",
#     values = uniprot_ids,
#     mart = mfa_mart
#   )
#
#
# converted_id_gn_2_ensesmbl <- converted_id_gn_2_ensesmbl |>
#   rename(ensembl = ensembl_gene_id, ID = external_gene_name)
#
# converted_id_uniprot_2_ensembl <- converted_id_uniprot_2_ensembl |>
#   rename(ensembl = ensembl_gene_id, ID = uniprotsptrembl)
#
# converted_id <- rbind(converted_id_uniprot_2_ensembl, converted_id_gn_2_ensesmbl) |>
#   distinct(ID, .keep_all = TRUE)
#
# p_dt_cleaned_up <- p_dt_cleaned |>
#   dplyr::filter(Cluster == "Cluster U") |>
#   distinct(ID, .keep_all = TRUE)
#
# p_dt_cleaned_up_converted <-
#   p_dt_cleaned_up |>
#   left_join(converted_id, by = "ID")
#
# p_dt_cleaned_down <- p_dt_cleaned |>
#   dplyr::filter(Cluster == "Cluster D") |>
#   distinct(ID, .keep_all = TRUE)
#
# p_dt_cleaned_down_converted <-
#   p_dt_cleaned_down |>
#   left_join(converted_id, by = "ID")
#
# sum(is.na(p_dt_cleaned_up_converted$ensembl))
# sum(is.na(p_dt_cleaned_down_converted$ensembl))
#
# unconverted_id_down <- p_dt_cleaned_down_converted$ID[is.na(p_dt_cleaned_down_converted$ensembl)]
# # ORF name starting with EGM_ cannot match with ensembl id
# unconverted_id_down_except_orf <- unconverted_id_down[!grepl("^EGM_", unconverted_id_down)]
# paste0(unconverted_id_down_except_orf, collapse = ",")
# # manually searched https://www.bgee.org/search/genes for the ensembl id using gene name
#
# unconverted_id_up <- p_dt_cleaned_up_converted$ID[is.na(p_dt_cleaned_up_converted$ensembl)]
# unconverted_id_up_except_orf <- unconverted_id_up[!grepl("^EGM_", unconverted_id_up)]
# paste0(unconverted_id_up_except_orf, collapse = ",")
#
# p_dt_cleaned_down_converted <- p_dt_cleaned_down_converted |> dplyr::filter(!is.na(ensembl))
# save(p_dt_cleaned_down_converted, file = "02_proteomics/protein_dt_cleaned_down_with_ensembl.rda")
#
# p_dt_cleaned_up_converted <- p_dt_cleaned_up_converted |> dplyr::filter(!is.na(ensembl))
# save(p_dt_cleaned_up_converted, file = "02_proteomics/protein_dt_cleaned_up_with_ensembl.rda")

# 03_metabolomic_data ====
m_dt <- import("2_data/case_study/01_monkey/mmc2.xlsx",
               skip = 1,
               header = TRUE,
               sheet = 5)

m_dt_match_with_keggid <- m_dt |> dplyr::mutate(keggid = NA, kegg_cpd_name = NA)

for (i in 1:nrow(m_dt)) {
  query_cpd_name <- gsub("\\**", "", m_dt$cleaned_cpd_name[i])
  query_res <- KEGGREST::keggFind("compound", query_cpd_name)
  if (length(query_res) == 0) {
    m_dt_match_with_keggid$keggid[i] <- NA
    m_dt_match_with_keggid$kegg_cpd_name[i] <- NA
  } else {
    m_dt_match_with_keggid$keggid[i] <- sub("cpd:", "", names(query_res[1]))
    m_dt_match_with_keggid$kegg_cpd_name[i] <- query_res[1]
  }
}

m_dt_up <- m_dt_match_with_keggid |> dplyr::filter(Cluster == "Cluster U")
sum(is.na(m_dt_up$keggid))
m_dt_up <- m_dt_up |> dplyr::filter(!is.na(keggid)) |> distinct(keggid, .keep_all = TRUE)

m_dt_down <- m_dt_match_with_keggid |> dplyr::filter(Cluster == "Cluster D")
sum(is.na(m_dt_down$keggid))
m_dt_down <- m_dt_down |> dplyr::filter(!is.na(keggid)) |> distinct(keggid, .keep_all = TRUE)

m_dt_u_d <- m_dt_match_with_keggid |> dplyr::filter(Cluster %in% c("Cluster U", "Cluster D"))

export(m_dt_u_d, file = "2_data/case_study/01_monkey/03_metabolomics/m_dt_u_d.xlsx")
# 04_sn_rna_seq_brain ====
brain_rna_seq_dt <- import(
  "2_data/case_study/01_monkey/04_sn_rnaseq/mmc4.xlsx",
  skip = 1,
  header = TRUE,
  sheet = 1
)

setwd("2_data/case_study/01_monkey/")

brain_up_dt <- brain_rna_seq_dt |>
  dplyr::filter((Group == "Young" & avg_log2FC < 0) | (Group == "Old" & avg_log2FC > 0))
range(brain_up_dt$p_val_adj)
range(brain_rna_seq_dt$avg_log2FC)

brain_up_dt_ensembl <-
  brain_up_dt |>
  dplyr::filter(grepl("^ENSMFAG", Gene))

ah <- AnnotationHub::AnnotationHub()
mf.orgdb <- ah[["AH119900"]]

converted_id <- clusterProfiler::bitr(
  geneID = brain_up_dt_ensembl$Gene,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = mf.orgdb
)

converted_id <- converted_id |> distinct(ENSEMBL, .keep_all = TRUE)
converted_id <- converted_id |> dplyr::rename(Gene = ENSEMBL)

brain_up_dt_ensembl <- brain_up_dt_ensembl |> left_join(converted_id, by = "Gene")
brain_up_dt_ensembl <- brain_up_dt_ensembl |>
  dplyr::filter(!is.na(SYMBOL)) |>
  group_by(Gene) |>
  arrange(p_val_adj) |>
  distinct(Gene, .keep_all = TRUE) |>
  ungroup()

brain_up_dt_ensembl <- brain_up_dt_ensembl |> mutate(Gene = SYMBOL) |> dplyr::select(-SYMBOL)

brain_up_dt_converted <- brain_up_dt |>
  dplyr::filter(!grepl("^ENSMFAG", Gene)) |>
  rbind(brain_up_dt_ensembl)

save(brain_up_dt_converted, file = "04_sn_rnaseq/brain_up_dt_converted.rda")

## down
brain_down_dt <- brain_rna_seq_dt |>
  dplyr::filter((Group == "Young" & avg_log2FC > 0) | (Group == "Old" & avg_log2FC < 0))
range(brain_down_dt$p_val_adj)
range(brain_down_dt$avg_log2FC)

brain_down_dt_ensembl <-
  brain_down_dt |>
  dplyr::filter(grepl("^ENSMFAG", Gene))

ah <- AnnotationHub::AnnotationHub()
mf.orgdb <- ah[["AH119900"]]

converted_id <- clusterProfiler::bitr(
  geneID = brain_down_dt_ensembl$Gene,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = mf.orgdb
)

converted_id <- converted_id |> distinct(ENSEMBL, .keep_all = TRUE)
converted_id <- converted_id |> dplyr::rename(Gene = ENSEMBL)

brain_down_dt_ensembl <- brain_down_dt_ensembl |> left_join(converted_id, by = "Gene")
brain_down_dt_ensembl <- brain_down_dt_ensembl |>
  dplyr::filter(!is.na(SYMBOL)) |>
  group_by(Gene) |>
  arrange(p_val_adj) |>
  distinct(Gene, .keep_all = TRUE) |>
  ungroup()

brain_down_dt_ensembl <- brain_down_dt_ensembl |> mutate(Gene = SYMBOL) |> dplyr::select(-SYMBOL)

brain_down_dt_converted <- brain_down_dt |>
  dplyr::filter(!grepl("^ENSMFAG", Gene)) |>
  rbind(brain_down_dt_ensembl)

sum(grepl("^ENSMFAG", brain_down_dt_converted$Gene))

save(brain_down_dt_converted, file = "04_sn_rnaseq/brain_down_dt_converted.rda")

# 04_sn_rna_seq_liver ====
liver_rna_seq_dt <- import(
  "2_data/case_study/01_monkey/04_sn_rnaseq/mmc4.xlsx",
  skip = 1,
  header = TRUE,
  sheet = 2
)

setwd("2_data/case_study/01_monkey/")

liver_up_dt <- liver_rna_seq_dt |>
  dplyr::filter((Group == "Young" & avg_log2FC < 0) | (Group == "Old" & avg_log2FC > 0))
range(liver_up_dt$p_val_adj)
range(liver_rna_seq_dt$avg_log2FC)

liver_up_dt_ensembl <-
  liver_up_dt |>
  dplyr::filter(grepl("^ENSMFAG", Gene))

converted_id <- clusterProfiler::bitr(
  geneID = liver_up_dt_ensembl$Gene,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = mf.orgdb
)

converted_id <- converted_id |> distinct(ENSEMBL, .keep_all = TRUE)
converted_id <- converted_id |> dplyr::rename(Gene = ENSEMBL)

liver_up_dt_ensembl <- liver_up_dt_ensembl |> left_join(converted_id, by = "Gene")
liver_up_dt_ensembl <- liver_up_dt_ensembl |>
  dplyr::filter(!is.na(SYMBOL)) |>
  group_by(Gene) |>
  arrange(p_val_adj) |>
  distinct(Gene, .keep_all = TRUE) |>
  ungroup()

liver_up_dt_ensembl <- liver_up_dt_ensembl |> mutate(Gene = SYMBOL) |> dplyr::select(-SYMBOL)

liver_up_dt_converted <- liver_up_dt |>
  dplyr::filter(!grepl("^ENSMFAG", Gene)) |>
  rbind(liver_up_dt_ensembl)

save(liver_up_dt_converted, file = "04_sn_rnaseq/liver_up_dt_converted.rda")

## down
liver_down_dt <- liver_rna_seq_dt |>
  dplyr::filter((Group == "Young" & avg_log2FC > 0) | (Group == "Old" & avg_log2FC < 0))
range(liver_down_dt$p_val_adj)
range(liver_down_dt$avg_log2FC)

liver_down_dt_ensembl <-
  liver_down_dt |>
  dplyr::filter(grepl("^ENSMFAG", Gene))

converted_id <- clusterProfiler::bitr(
  geneID = liver_down_dt_ensembl$Gene,
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = mf.orgdb
)

converted_id <- converted_id |> distinct(ENSEMBL, .keep_all = TRUE)
converted_id <- converted_id |> dplyr::rename(Gene = ENSEMBL)

liver_down_dt_ensembl <- liver_down_dt_ensembl |> left_join(converted_id, by = "Gene")
liver_down_dt_ensembl <- liver_down_dt_ensembl |>
  dplyr::filter(!is.na(SYMBOL)) |>
  group_by(Gene) |>
  arrange(p_val_adj) |>
  distinct(Gene, .keep_all = TRUE) |>
  ungroup()

liver_down_dt_ensembl <- liver_down_dt_ensembl |> mutate(Gene = SYMBOL) |> dplyr::select(-SYMBOL)

liver_down_dt_converted <- liver_down_dt |>
  dplyr::filter(!grepl("^ENSMFAG", Gene)) |>
  rbind(liver_down_dt_ensembl)

sum(grepl("^ENSMFAG", liver_down_dt_converted$Gene))

save(liver_down_dt_converted, file = "04_sn_rnaseq/liver_down_dt_converted.rda")
