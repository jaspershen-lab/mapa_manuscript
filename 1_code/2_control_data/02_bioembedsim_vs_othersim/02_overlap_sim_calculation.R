library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

library(org.Hs.eg.db)
library(AnnotationDbi)
library(reactome.db)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

setwd("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/overlap_similarity")

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
# save(all_gene_list, file = "all_gene_list.rda")
load("all_gene_list.rda")

# 2. Calculate similarity ====
jaccard_sim_matrix <- term_similarity_internal(gl = all_gene_list,
                                               measure.method = "jaccard")
jaccard_sim_df <-
  as.data.frame.table(jaccard_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

save(jaccard_sim_matrix, file = "jaccard_sim_matrix.rda")
save(jaccard_sim_df, file = "jaccard_sim_df.rda")
load("jaccard_sim_df.rda")

op_sim_matrix <- term_similarity_internal(gl = all_gene_list,
                                          measure.method = "overlap")
op_sim_df <-
  as.data.frame.table(op_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)
save(op_sim_matrix, file = "op_sim_matrix.rda")
save(op_sim_df, file = "op_sim_df.rda")

kappa_sim_matrix <- term_similarity_internal(gl = all_gene_list,
                                             measure.method = "kappa")

kappa_sim_df <-
  as.data.frame.table(kappa_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

save(kappa_sim_matrix, file = "kappa_sim_matrix.rda")
save(kappa_sim_df, file = "kappa_sim_df.rda")
load("kappa_sim_df.rda")

dice_sim_matrix <- term_similarity_internal(gl = all_gene_list,
                                            measure.method = "dice")
dice_sim_df <-
  as.data.frame.table(dice_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

save(dice_sim_matrix, file = "dice_sim_matrix.rda")
save(dice_sim_df, file = "dice_sim_df.rda")
load("dice_sim_df.rda")

###

###Jaccard index
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  pull(id)

temp_data <-
  jaccard_sim_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::select(pathway, all_of(remain_pathways)) %>%
  dplyr::filter(pathway %in% remain_pathways) %>%
  tibble::column_to_rownames(var = "pathway")

colnames(temp_data) <- NULL

pca_res <-
  prcomp(temp_data, center = TRUE, scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x)

pca_scores <-
  pca_scores %>%
  dplyr::select(PC1, PC2) %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  dplyr::mutate(
    expected_module = factor(expected_module, levels = stringr::str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

plot <-
  pca_scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = expected_module, shape = database), size = 6) +
  theme_bw() +
  stat_ellipse(
    aes(color = expected_module),
    type = "t",
    geom = "polygon",
    fill = NA
  ) +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right")
plot
ggsave(plot,
       filename = "jaccard_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(temp_data, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(unique(control_dt$expected_module), numeric = TRUE)),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database), size = 6, alpha = 0.8) +
  stat_ellipse(
    aes(color = expected_module),
    type = "t",
    geom = "polygon",
    fill = NA
  ) +
  theme_bw() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right")

# Show the plot
umap_plot
ggsave(umap_plot,
       filename = "jaccard_umap_plot.pdf",
       width = 9,
       height = 7)






###Overlap coefficient
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  pull(id)

temp_data <-
  op_sim_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::select(pathway, all_of(remain_pathways)) %>%
  dplyr::filter(pathway %in% remain_pathways) %>%
  tibble::column_to_rownames(var = "pathway")

colnames(temp_data) <- NULL

pca_res <-
  prcomp(temp_data, center = TRUE, scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x)

pca_scores <-
  pca_scores %>%
  dplyr::select(PC1, PC2) %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  dplyr::mutate(
    expected_module = factor(expected_module, levels = stringr::str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

plot <-
  pca_scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = expected_module, shape = database), size = 6) +
  theme_bw() +
  stat_ellipse(
    aes(color = expected_module),
    type = "t",
    geom = "polygon",
    fill = NA
  ) +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right")
plot
ggsave(plot,
       filename = "op_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(temp_data, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(unique(control_dt$expected_module), numeric = TRUE)),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database), size = 6, alpha = 0.8) +
  stat_ellipse(
    aes(color = expected_module),
    type = "t",
    geom = "polygon",
    fill = NA
  ) +
  theme_bw() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right")

# Show the plot
umap_plot
ggsave(umap_plot,
       filename = "op_umap_plot.pdf",
       width = 9,
       height = 7)






###Kappa
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  pull(id)

temp_data <-
  kappa_sim_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::select(pathway, all_of(remain_pathways)) %>%
  dplyr::filter(pathway %in% remain_pathways) %>%
  tibble::column_to_rownames(var = "pathway")

colnames(temp_data) <- NULL

pca_res <-
  prcomp(temp_data, center = TRUE, scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x)

pca_scores <-
  pca_scores %>%
  dplyr::select(PC1, PC2) %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  dplyr::mutate(
    expected_module = factor(expected_module, levels = stringr::str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

plot <-
  pca_scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = expected_module, shape = database), size = 6) +
  theme_bw() +
  stat_ellipse(
    aes(color = expected_module),
    type = "t",
    geom = "polygon",
    fill = NA
  ) +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right")
plot
ggsave(plot,
       filename = "kappa_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(temp_data, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(unique(control_dt$expected_module), numeric = TRUE)),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database), size = 6, alpha = 0.8) +
  stat_ellipse(
    aes(color = expected_module),
    type = "t",
    geom = "polygon",
    fill = NA
  ) +
  theme_bw() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right")

# Show the plot
umap_plot
ggsave(umap_plot,
       filename = "kappa_umap_plot.pdf",
       width = 9,
       height = 7)






###Dice
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  pull(id)

temp_data <-
  dice_sim_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::select(pathway, all_of(remain_pathways)) %>%
  dplyr::filter(pathway %in% remain_pathways) %>%
  tibble::column_to_rownames(var = "pathway")

colnames(temp_data) <- NULL

pca_res <-
  prcomp(temp_data, center = TRUE, scale. = TRUE)

pca_scores <- as.data.frame(pca_res$x)

pca_scores <-
  pca_scores %>%
  dplyr::select(PC1, PC2) %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  dplyr::mutate(
    expected_module = factor(expected_module, levels = stringr::str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

plot <-
  pca_scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = expected_module, shape = database), size = 6) +
  theme_bw() +
  stat_ellipse(
    aes(color = expected_module),
    type = "t",
    geom = "polygon",
    fill = NA
  ) +
  labs(x = "PC1", y = "PC2") +
  theme(legend.position = "right")
plot
ggsave(plot,
       filename = "dice_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(temp_data, n_neighbors = 15, min_dist = 0.1, metric = "cosine")

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(unique(control_dt$expected_module), numeric = TRUE)),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database), size = 6, alpha = 0.8) +
  stat_ellipse(
    aes(color = expected_module),
    type = "t",
    geom = "polygon",
    fill = NA
  ) +
  theme_bw() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right")

# Show the plot
umap_plot
ggsave(umap_plot,
       filename = "dice_umap_plot.pdf",
       width = 9,
       height = 7)

