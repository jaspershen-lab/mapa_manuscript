library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

library(org.Hs.eg.db)
library(AnnotationDbi)
library(reactome.db)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

# 1. Get similarity based on gene overlap ====
## Get gene id for each pathway
# go_id <-
#   control_dt |>
#   dplyr::filter(database == "GO") |>
#   dplyr::pull(id)
#
# go2egs <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go_id, columns = "ENTREZID", keytype = "GOALL"))
#
# go_gene_lists <- list()
#
# go_pathways <-
#   control_dt %>%
#   dplyr::filter(database == "GO")
#
# for (i in 1:nrow(go_pathways)) {
#   pathway_id <- go_pathways$id[i]
#
#   # Add this gene set to the list, using the pathway_id as the name
#   go_gene_lists[[pathway_id]] <-
#     go2egs %>%
#     dplyr::filter(GOALL == pathway_id) %>%
#     pull(ENTREZID) %>%
#     unique()
# }
#
# ## reactome
# reactome_id <-
#   control_dt %>%
#   dplyr::filter(database == "Reactome") %>%
#   dplyr::pull(id)
#
# reactome2egs <- AnnotationDbi::select(reactome.db, keys = reactome_id, columns = "ENTREZID", keytype = "PATHID")
#
# reactome_gene_lists <- list()
#
# reactome_pathways <-
#   control_dt %>%
#   dplyr::filter(database == "Reactome")
#
# for (i in 1:nrow(reactome_pathways)) {
#   pathway_id <- reactome_pathways$id[i]
#
#   # Add this gene set to the list, using the pathway_id as the name
#   reactome_gene_lists[[pathway_id]] <-
#     reactome2egs %>%
#     dplyr::filter(PATHID == pathway_id) %>%
#     pull(ENTREZID) %>%
#     unique()
# }
#
# ## KEGG
# kegg_id <-
#   control_dt %>%
#   dplyr::filter(database == "KEGG") %>%
#   dplyr::pull(id)
#
# chunk_size <- 10
# chunks <- split(kegg_id, ceiling(seq_along(kegg_id) / chunk_size))
# kegg_info <- list()
#
# for (i in 1:length(chunks)) {
#   sub_kegg_info <-
#     KEGGREST::keggGet(dbentries = chunks[[i]]) %>%
#     purrr::map(
#       function(x) {
#         # Initialize kegg2genename as empty
#         kegg2genename <- NA
#
#         # Only get annotated Entrez ID if include_gene_name is TRUE
#         if ("GENE" %in% names(x)) {
#           kegg2genename <-
#             x$GENE[seq(1, length(x$GENE), 2)]
#         }
#
#         all_info <- list(
#           "id" = unname(x$ENTRY),
#           "genes" = kegg2genename
#         )
#
#         return(all_info)
#       }
#     )
#   kegg_info <- c(kegg_info, sub_kegg_info)
# }
#
# kegg_gene_lists <- list()
#
# for (i in 1:length(kegg_info)) {
#   pathway_id <- kegg_info[[i]]$id
#
#   kegg_gene_lists[[pathway_id]] <- kegg_info[[i]]$genes
# }
#
# all_gene_list <- c(go_gene_lists, kegg_gene_lists, reactome_gene_lists)
#
# save(all_gene_list, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/all_gene_list.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/all_gene_list.rda")

# 2. Calculate similarity ====
# jaccard_sim_matrix <- term_similarity_internal(gl = all_gene_list,
#                                                measure.method = "jaccard")
# jaccard_sim_df <-
#   as.data.frame.table(jaccard_sim_matrix, responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(jaccard_sim_df, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/jaccard_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/jaccard_sim_df.rda")

# op_sim_matrix <- term_similarity_internal(gl = all_gene_list,
#                                           measure.method = "overlap")
# op_sim_df <-
#   as.data.frame.table(op_sim_matrix, responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(op_sim_df, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/op_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/op_sim_df.rda")

# kappa_sim_matrix <- term_similarity_internal(gl = all_gene_list,
#                                              measure.method = "kappa")
#
# kappa_sim_df <-
#   as.data.frame.table(kappa_sim_matrix, responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(kappa_sim_df, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/kappa_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/kappa_sim_df.rda")

# dice_sim_matrix <- term_similarity_internal(gl = all_gene_list,
#                                             measure.method = "dice")
# dice_sim_df <-
#   as.data.frame.table(dice_sim_matrix, responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(dice_sim_df, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/dice_sim_df.rda")
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/dice_sim_df.rda")
