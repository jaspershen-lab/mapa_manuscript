library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

# 1. Get annotated gene information for pathways ====
## 22 GO terms
go_info <-
  control_dt |>
  dplyr::filter(database == "GO") |>
  dplyr::pull(id) |>
  get_go_info()

## 14 KEGG pathways
kegg_info <-
  control_dt |>
  dplyr::filter(database == "KEGG") |>
  dplyr::pull(id) |>
  get_kegg_pathway_info()

## 8 Ractome pathways
reactome_info <-
  control_dt |>
  dplyr::filter(database == "Reactome") |>
  dplyr::pull(id) |>
  get_reactome_pathway_info()

all_text_info <- c(go_info, kegg_info, reactome_info)
all_combined_info <- combine_info(info = all_text_info)

# 2. Get semantic similarity using embedding model ====
embedding_matrix <- get_embedding_matrix(
  text = all_combined_info,
  api_provider = "openai",
  text_embedding_model = "text-embedding-3-small",
  api_key = api_key
)
save(embedding_matrix, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/embedding_matrix.rda")

embedding_sim_matrix <- calculate_cosine_sim(m = embedding_matrix)
save(embedding_sim_matrix, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/embedding_sim_matrix.rda")

embedding_sim_df <-
  as.data.frame.table(embedding_sim_matrix, responseName = "sim") |>
  dplyr::filter(Var1 != Var2) |>                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) |>
  dplyr::mutate(across(c(from, to), as.character)) |>
  dplyr::filter(from < to)
save(embedding_sim_df, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/embedding_sim_df.rda")
