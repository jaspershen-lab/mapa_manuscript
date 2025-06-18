library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

load("2_data/curated_pathway_data_2/pathways_table.rda")
setDT(pathways_table)
load("2_data/curated_pathway_data_2/hmdb_kegg_matched.rda")
load("2_data/curated_pathway_data_2/hmdb_reactome_matched.rda")
load("2_data/curated_pathway_data_2/kegg_reactome_matched.rda")

get_all_met_sim_results <- function(matched_df, db1, db2, api_key) {
  matched_df_with_sim <- matched_df |> mutate(biotext_sim = 0,
                                              jaccard_sim = 0,
                                              overlap_sim = 0,
                                              kappa_sim = 0,
                                              dice_sim = 0)
  for (i in 1:nrow(matched_df)) {
    pathway_pair <- matched_df[i,]
    db1_id <- pathway_pair[[db1]]
    db2_id <- pathway_pair[[db2]]

    # Calculate biotext embedding similarity
    db1_info <- switch(db1,
                       "hmdb_pathway_id" = get_smpdb_info(db1_id),
                       "kegg_pathway_id" = get_kegg_pathway_info(db1_id),
                       "reactome_pathway_id" = get_reactome_pathway_info(db1_id))
    db2_info <- switch(db2,
                       "hmdb_pathway_id" = get_smpdb_info(db2_id),
                       "kegg_pathway_id" = get_kegg_pathway_info(db2_id),
                       "reactome_pathway_id" = get_reactome_pathway_info(db2_id))

    all_combined_info <- combine_info(c(db1_info, db2_info))

    embedding_matrix <- get_embedding_matrix(
      text = all_combined_info,
      api_provider = "openai",
      text_embedding_model = "text-embedding-3-small",
      api_key = api_key
    )

    embedding_sim_matrix <- calculate_cosine_sim(m = embedding_matrix)
    embedding_sim_df <-
      as.data.frame.table(embedding_sim_matrix, responseName = "sim") |>
      dplyr::filter(Var1 != Var2) |>                 # Remove self-edges
      dplyr::rename(from = Var1, to = Var2) |>
      dplyr::mutate(across(c(from, to), as.character)) |>
      dplyr::filter(from < to)

    matched_df_with_sim$biotext_sim[i] <- embedding_sim_df$sim

    # Calculate similarity based on metabolite overlap
    met_list <- list()
    met_list[[db1_id]] <- pathways_table[pathway_id == db1_id, HMDB_ID] |> strsplit("\\{\\}") |> unlist()
    met_list[[db2_id]] <- pathways_table[pathway_id == db2_id, HMDB_ID] |> strsplit("\\{\\}") |> unlist()

    jaccard_sim_matrix <- term_similarity_internal(gl = met_list,
                                                   measure.method = "jaccard")
    jaccard_sim_df <-
      as.data.frame.table(jaccard_sim_matrix, responseName = "sim") %>%
      dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
      dplyr::rename(from = Var1, to = Var2) %>%
      dplyr::mutate(across(c(from, to), as.character)) %>%
      dplyr::filter(from < to)
    matched_df_with_sim$jaccard_sim[i] <- jaccard_sim_df$sim

    op_sim_matrix <- term_similarity_internal(gl = met_list,
                                              measure.method = "overlap")
    op_sim_df <-
      as.data.frame.table(op_sim_matrix, responseName = "sim") %>%
      dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
      dplyr::rename(from = Var1, to = Var2) %>%
      dplyr::mutate(across(c(from, to), as.character)) %>%
      dplyr::filter(from < to)
    matched_df_with_sim$overlap_sim[i] <- op_sim_df$sim

    kappa_sim_matrix <- term_similarity_internal(gl = met_list,
                                                 measure.method = "kappa")
    kappa_sim_df <-
      as.data.frame.table(kappa_sim_matrix, responseName = "sim") %>%
      dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
      dplyr::rename(from = Var1, to = Var2) %>%
      dplyr::mutate(across(c(from, to), as.character)) %>%
      dplyr::filter(from < to)
    matched_df_with_sim$kappa_sim[i] <- kappa_sim_df$sim

    dice_sim_matrix <- term_similarity_internal(gl = met_list,
                                                measure.method = "dice")
    dice_sim_df <-
      as.data.frame.table(dice_sim_matrix, responseName = "sim") %>%
      dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
      dplyr::rename(from = Var1, to = Var2) %>%
      dplyr::mutate(across(c(from, to), as.character)) %>%
      dplyr::filter(from < to)
    matched_df_with_sim$dice_sim[i] <- dice_sim_df$sim
  }

  matched_df_with_sim
}

hmdb_kegg_matched_with_sim <- get_all_met_sim_results(matched_df = hmdb_kegg_matched,
                                                      db1 = "hmdb_pathway_id",
                                                      db2 = "kegg_pathway_id",
                                                      api_key = api_key)

save(hmdb_kegg_matched_with_sim, file = "3_data_analysis/03_curated_pathway_dataset2/hmdb_kegg_matched_with_sim.rda")

hmdb_reactome_matched_with_sim <- get_all_met_sim_results(matched_df = hmdb_reactome_matched,
                                                          db1 = "hmdb_pathway_id",
                                                          db2 = "reactome_pathway_id",
                                                          api_key = api_key)

save(hmdb_reactome_matched_with_sim, file = "3_data_analysis/03_curated_pathway_dataset2/hmdb_reactome_matched_with_sim.rda")

kegg_reactome_matched_with_sim <- get_all_met_sim_results(matched_df = kegg_reactome_matched,
                                                          db1 = "kegg_pathway_id",
                                                          db2 = "reactome_pathway_id",
                                                          api_key = api_key)
save(kegg_reactome_matched_with_sim, file = "3_data_analysis/03_curated_pathway_dataset2/kegg_reactome_matched_with_sim.rda")
