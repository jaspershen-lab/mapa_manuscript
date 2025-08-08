library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# sn_rna_seq_brain ====
brain_rna_seq_dt <- import(
  "2_data/case_study/01_monkey/mmc4.xlsx",
  skip = 1,
  header = TRUE,
  sheet = 1
)

dir.create("3_data_analysis/05_case_study/01_monkey/brain_sn_rna_seq")
setwd("3_data_analysis/05_case_study/01_monkey/brain_sn_rna_seq")

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

save(brain_up_dt_converted, file = "brain_up_dt_converted.rda")

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

save(brain_down_dt_converted, file = "brain_down_dt_converted.rda")

# sn_rna_seq_liver ====
liver_rna_seq_dt <- import(
  "2_data/case_study/01_monkey/04_sn_rnaseq/mmc4.xlsx",
  skip = 1,
  header = TRUE,
  sheet = 2
)

dir.create("3_data_analysis/05_case_study/01_monkey/liver_sn_rna_seq")
setwd("3_data_analysis/05_case_study/01_monkey/liver_sn_rna_seq")

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

save(liver_up_dt_converted, file = "liver_up_dt_converted.rda")

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

save(liver_down_dt_converted, file = "liver_down_dt_converted.rda")
