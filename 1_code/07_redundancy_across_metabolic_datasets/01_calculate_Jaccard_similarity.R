library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)

# load three metabolic pathway datasets
load("2_data/hmdb_pathway.rda")
load("2_data/kegg_hsa_pathway.rda")
load("2_data/reactome_human_pathway.rda")

smpdb_database <- hmdb_pathway
kegg_database <- kegg_hsa_pathway
reactome_database <- reactome_human_pathway_final

# remove pathways with no HMDB IDs
smpdb_database_filtered <- smpdb_database[!is.na(smpdb_database$HMDB_ID), ]
kegg_database_filtered <- kegg_database[!is.na(kegg_database$HMDB_ID), ]
reactome_database_filtered <- reactome_database[!is.na(reactome_database$HMDB_ID), ]

# function to split hmdb id
split_hmdb <- function(hmdb_str) {
  hmdb_str %>%
    str_split("\\{\\}") %>%
    unlist() %>%
    unique()
}

# get the list of HMDB IDs for each pathway
smpdb_list <- smpdb_database_filtered %>%
  mutate(HMDB_vec = map(HMDB_ID, split_hmdb)) %>%
  transmute(
    pathway_id_smpdb = pathway_id,
    HMDB_ID_smpdb   = HMDB_ID,
    HMDB_vec_smpdb  = HMDB_vec
  )


kegg_list <- kegg_database_filtered %>%
  mutate(HMDB_vec = map(HMDB_ID, split_hmdb)) %>%
  transmute(
    pathway_id_kegg = pathway_id,
    HMDB_ID_kegg   = HMDB_ID,
    HMDB_vec_kegg  = HMDB_vec
  )

reactome_list <- reactome_database_filtered %>%
  mutate(HMDB_vec = map(HMDB_ID, split_hmdb)) %>%
  transmute(
    pathway_id_reactome = pathway_id,
    HMDB_ID_reactome   = HMDB_ID,
    HMDB_vec_reactome  = HMDB_vec
  )


# function to calculate the Jaccard similarity
jaccard_similarity <- function(vec1, vec2) {
  inter <- length(intersect(vec1, vec2))
  union <- length(union(vec1, vec2))
  if (union == 0) return(NA_real_)
  inter / union
}

# result for pathways
smpdb_kegg_result <- expand_grid(smpdb_list, kegg_list) %>%
  mutate(similarity = map2_dbl(HMDB_vec_smpdb, HMDB_vec_kegg, jaccard_similarity)) %>%
  select(
    pathway_id_smpdb, HMDB_ID_smpdb,
    pathway_id_kegg,  HMDB_ID_kegg,
    similarity
  )

reactome_kegg_result <- expand_grid(reactome_list, kegg_list) %>%
  mutate(similarity = map2_dbl(HMDB_vec_reactome, HMDB_vec_kegg, jaccard_similarity)) %>%
  select(
    pathway_id_reactome, HMDB_ID_reactome,
    pathway_id_kegg,  HMDB_ID_kegg,
    similarity
  )

smpdb_reactome_result <- expand_grid(smpdb_list, reactome_list) %>%
  mutate(similarity = map2_dbl(HMDB_vec_smpdb, HMDB_vec_reactome, jaccard_similarity)) %>%
  select(
    pathway_id_smpdb, HMDB_ID_smpdb,
    pathway_id_reactome,  HMDB_ID_reactome,
    similarity
  )

# save the results
save(smpdb_kegg_result, file = "2_data/smpdb_kegg_similarity.rda")
save(reactome_kegg_result, file = "2_data/reactome_kegg_similarity.rda")
save(smpdb_reactome_result, file = "2_data/smpdb_reactome_similarity.rda")


