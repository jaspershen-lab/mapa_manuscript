library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
# source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

library(GOSemSim)
library(org.Hs.eg.db)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

setwd("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/GO_terms/")

# ## BP
# go_bp_id <- control_dt %>% dplyr::filter(database == "GO" &
#                                            go_ontology == "BP") %>% pull(id)
# bp_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "BP")
#
# ##Wang
# bp_sim_matrix_wang <-
#   termSim(
#     t1 = go_bp_id,
#     t2 = go_bp_id,
#     semData = bp_semgodata,
#     method = "Wang"
#   )
#
# bp_sim_df_wang <-
#   bp_sim_matrix_wang %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(bp_sim_df_wang, file = "bp_sim_df_wang.rda")
# save(bp_sim_matrix_wang, file = "bp_sim_matrix_wang.rda")

load("bp_sim_df_wang.rda")
load("bp_sim_matrix_wang.rda")

# ##Resnik
# bp_sim_matrix_resnik <-
#   termSim(
#     t1 = go_bp_id,
#     t2 = go_bp_id,
#     semData = bp_semgodata,
#     method = "Resnik"
#   )
#
# bp_sim_df_resnik <-
#   bp_sim_matrix_resnik %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(bp_sim_df_resnik, file = "bp_sim_df_resnik.rda")
# save(bp_sim_matrix_resnik, file = "bp_sim_matrix_resnik.rda")
load("bp_sim_df_resnik.rda")
load("bp_sim_matrix_resnik.rda")

# ##Rel
# bp_sim_matrix_rel <-
#   termSim(
#     t1 = go_bp_id,
#     t2 = go_bp_id,
#     semData = bp_semgodata,
#     method = "Rel"
#   )
#
# bp_sim_df_rel <-
#   bp_sim_matrix_rel %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(bp_sim_df_rel, file = "bp_sim_df_rel.rda")
# save(bp_sim_matrix_rel, file = "bp_sim_matrix_rel.rda")

load("bp_sim_df_rel.rda")
load("bp_sim_matrix_rel.rda")

# ##Jiang
# bp_sim_matrix_jiang <-
#   termSim(
#     t1 = go_bp_id,
#     t2 = go_bp_id,
#     semData = bp_semgodata,
#     method = "Jiang"
#   )
#
# bp_sim_df_jiang <-
#   bp_sim_matrix_jiang %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(bp_sim_df_jiang, file = "bp_sim_df_jiang.rda")
# save(bp_sim_matrix_jiang, file = "bp_sim_matrix_jiang.rda")

load("bp_sim_df_jiang.rda")
load("bp_sim_matrix_jiang.rda")

# ##Lin
# bp_sim_matrix_lin <-
#   termSim(
#     t1 = go_bp_id,
#     t2 = go_bp_id,
#     semData = bp_semgodata,
#     method = "Lin"
#   )
#
# bp_sim_df_lin <-
#   bp_sim_matrix_lin %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(bp_sim_df_lin, file = "bp_sim_df_lin.rda")
# save(bp_sim_matrix_lin, file = "bp_sim_matrix_lin.rda")
load("bp_sim_df_lin.rda")
load("bp_sim_matrix_lin.rda")

# ##TCSS
# bp_semgodata_tcss <- godata(annoDb = org.Hs.eg.db, ont = "BP", processTCSS = TRUE)
# bp_sim_df_tcss <-
#   termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata_tcss, method = "TCSS") %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)

# ## MF
# go_mf_id <- control_dt %>% dplyr::filter(database == "GO" &
#                                            go_ontology == "MF") %>% pull(id)
# mf_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "MF")
#
# ##wang
# ##Wang
# mf_sim_matrix_wang <-
#   termSim(
#     t1 = go_mf_id,
#     t2 = go_mf_id,
#     semData = mf_semgodata,
#     method = "Wang"
#   )
#
# mf_sim_df_wang <-
#   mf_sim_matrix_wang %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(mf_sim_df_wang, file = "mf_sim_df_wang.rda")
# save(mf_sim_matrix_wang, file = "mf_sim_matrix_wang.rda")
load("mf_sim_df_wang.rda")
load("mf_sim_matrix_wang.rda")

# ##Resnik
# mf_sim_matrix_resnik <-
#   termSim(
#     t1 = go_mf_id,
#     t2 = go_mf_id,
#     semData = mf_semgodata,
#     method = "Resnik"
#   )
#
# mf_sim_df_resnik <-
#   mf_sim_matrix_resnik %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(mf_sim_df_resnik, file = "mf_sim_df_resnik.rda")
# save(mf_sim_matrix_resnik, file = "mf_sim_matrix_resnik.rda")
load("mf_sim_df_resnik.rda")
load("mf_sim_matrix_resnik.rda")

# ##Rel
# mf_sim_matrix_rel <-
#   termSim(
#     t1 = go_mf_id,
#     t2 = go_mf_id,
#     semData = mf_semgodata,
#     method = "Rel"
#   )
#
# mf_sim_df_rel <-
#   mf_sim_matrix_rel %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(mf_sim_df_rel, file = "mf_sim_df_rel.rda")
# save(mf_sim_matrix_rel, file = "mf_sim_matrix_rel.rda")
load("mf_sim_df_rel.rda")
load("mf_sim_matrix_rel.rda")

# ##Jiang
# mf_sim_matrix_jiang <-
#   termSim(
#     t1 = go_mf_id,
#     t2 = go_mf_id,
#     semData = mf_semgodata,
#     method = "Jiang"
#   )
#
# mf_sim_df_jiang <-
#   mf_sim_matrix_jiang %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(mf_sim_df_jiang, file = "mf_sim_df_jiang.rda")
# save(mf_sim_matrix_jiang, file = "mf_sim_matrix_jiang.rda")
load("mf_sim_df_jiang.rda")
load("mf_sim_matrix_jiang.rda")

# ##Lin
# mf_sim_matrix_lin <-
#   termSim(
#     t1 = go_mf_id,
#     t2 = go_mf_id,
#     semData = mf_semgodata,
#     method = "Lin"
#   )
#
# mf_sim_df_lin <-
#   mf_sim_matrix_lin %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(mf_sim_df_lin, file = "mf_sim_df_lin.rda")
# save(mf_sim_matrix_lin, file = "mf_sim_matrix_lin.rda")
load("mf_sim_df_lin.rda")
load("mf_sim_matrix_lin.rda")

# ##TCSS
# mf_semgodata_tcss <- godata(annoDb = org.Hs.eg.db, ont = "BP", processTCSS = TRUE)
# mf_sim_df_tcss <-
#   termSim(t1 = go_mf_id, t2 = go_mf_id, semData = mf_semgodata_tcss, method = "TCSS") %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)


# ## CC
# go_cc_id <- control_dt %>% dplyr::filter(database == "GO" &
#                                            go_ontology == "CC") %>% pull(id)
# cc_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "CC")
#
# ##Wang
# cc_sim_matrix_wang <-
#   termSim(
#     t1 = go_cc_id,
#     t2 = go_cc_id,
#     semData = cc_semgodata,
#     method = "Wang"
#   )
#
# cc_sim_df_wang <-
#   cc_sim_matrix_wang %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(cc_sim_df_wang, file = "cc_sim_df_wang.rda")
# save(cc_sim_matrix_wang, file = "cc_sim_matrix_wang.rda")
load("cc_sim_df_wang.rda")
load("cc_sim_matrix_wang.rda")

# ##Resnik
# cc_sim_matrix_resnik <-
#   termSim(
#     t1 = go_cc_id,
#     t2 = go_cc_id,
#     semData = cc_semgodata,
#     method = "Resnik"
#   )
#
# cc_sim_df_resnik <-
#   cc_sim_matrix_resnik %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(cc_sim_df_resnik, file = "cc_sim_df_resnik.rda")
# save(cc_sim_matrix_resnik, file = "cc_sim_matrix_resnik.rda")
load("cc_sim_df_resnik.rda")
load("cc_sim_matrix_resnik.rda")

# ##Rel
# cc_sim_matrix_rel <-
#   termSim(
#     t1 = go_cc_id,
#     t2 = go_cc_id,
#     semData = cc_semgodata,
#     method = "Rel"
#   )
#
# cc_sim_df_rel <-
#   cc_sim_matrix_rel %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(cc_sim_df_rel, file = "cc_sim_df_rel.rda")
# save(cc_sim_matrix_rel, file = "cc_sim_matrix_rel.rda")
load("cc_sim_df_rel.rda")
load("cc_sim_matrix_rel.rda")

# ##Jiang
# cc_sim_matrix_jiang <-
#   termSim(
#     t1 = go_cc_id,
#     t2 = go_cc_id,
#     semData = cc_semgodata,
#     method = "Jiang"
#   )
#
# cc_sim_df_jiang <-
#   cc_sim_matrix_jiang %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(cc_sim_df_jiang, file = "cc_sim_df_jiang.rda")
# save(cc_sim_matrix_jiang, file = "cc_sim_matrix_jiang.rda")
load("cc_sim_df_jiang.rda")
load("cc_sim_matrix_jiang.rda")

# ##Lin
# cc_sim_matrix_lin <-
#   termSim(
#     t1 = go_cc_id,
#     t2 = go_cc_id,
#     semData = cc_semgodata,
#     method = "Lin"
#   )
#
# cc_sim_df_lin <-
#   cc_sim_matrix_lin %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)
#
# save(cc_sim_df_lin, file = "cc_sim_df_lin.rda")
# save(cc_sim_matrix_lin, file = "cc_sim_matrix_lin.rda")
load("cc_sim_df_lin.rda")
load("cc_sim_matrix_lin.rda")

# ##TCSS
# cc_semgodata_tcss <- godata(annoDb = org.Hs.eg.db, ont = "BP", processTCSS = TRUE)
# cc_sim_df_tcss <-
#   termSim(t1 = go_cc_id, t2 = go_cc_id, semData = cc_semgodata_tcss, method = "TCSS") %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)

wang_sim_df <- rbind(bp_sim_df_wang, mf_sim_df_wang, cc_sim_df_wang)
resnik_sim_df <- rbind(bp_sim_df_resnik, mf_sim_df_resnik, cc_sim_df_resnik)
rel_sim_df <- rbind(bp_sim_df_rel, mf_sim_df_rel, cc_sim_df_rel)
jiang_sim_df <- rbind(bp_sim_df_jiang, mf_sim_df_jiang, cc_sim_df_jiang)
lin_sim_df <- rbind(bp_sim_df_lin, mf_sim_df_lin, cc_sim_df_lin)

save(wang_sim_df, file = "wang_sim_df.rda")
save(resnik_sim_df, file = "resnik_sim_df.rda")
save(rel_sim_df, file = "rel_sim_df.rda")
save(jiang_sim_df, file = "jiang_sim_df.rda")
save(lin_sim_df, file = "lin_sim_df.rda")

###Wang
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  dplyr::filter(database == "GO" & go_ontology == "BP") %>%
  pull(id)

temp_data <-
  bp_sim_matrix_wang %>%
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
       filename = "wang_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(
  temp_data,
  n_neighbors = 10,
  min_dist = 0.1,
  metric = "cosine"
)

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database),
             size = 6,
             alpha = 0.8) +
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
       filename = "wang_umap_plot.pdf",
       width = 9,
       height = 7)





###jiang
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  dplyr::filter(database == "GO" & go_ontology == "BP") %>%
  pull(id)

temp_data <-
  bp_sim_matrix_jiang %>%
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
       filename = "jiang_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(
  temp_data,
  n_neighbors = 10,
  min_dist = 0.1,
  metric = "cosine"
)

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database),
             size = 6,
             alpha = 0.8) +
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
       filename = "jiang_umap_plot.pdf",
       width = 9,
       height = 7)






###lin
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  dplyr::filter(database == "GO" & go_ontology == "BP") %>%
  pull(id)

temp_data <-
  bp_sim_matrix_lin %>%
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
       filename = "lin_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(
  temp_data,
  n_neighbors = 10,
  min_dist = 0.1,
  metric = "cosine"
)

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database),
             size = 6,
             alpha = 0.8) +
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
       filename = "lin_umap_plot.pdf",
       width = 9,
       height = 7)







###rel
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  dplyr::filter(database == "GO" & go_ontology == "BP") %>%
  pull(id)

temp_data <-
  bp_sim_matrix_rel %>%
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
       filename = "rel_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(
  temp_data,
  n_neighbors = 10,
  min_dist = 0.1,
  metric = "cosine"
)

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database),
             size = 6,
             alpha = 0.8) +
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
       filename = "rel_umap_plot.pdf",
       width = 9,
       height = 7)




###resnik
##PCA analysis
##only remain the modules with more than 2 pathways
remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  dplyr::filter(database == "GO" & go_ontology == "BP") %>%
  pull(id)

temp_data <-
  bp_sim_matrix_resnik %>%
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
       filename = "resnik_pca_plot.pdf",
       width = 9,
       height = 7)

###UMAP
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
umap_res <- umap(
  temp_data,
  n_neighbors = 10,
  min_dist = 0.1,
  metric = "cosine"
)

# Create a data frame of UMAP results
umap_scores <- as.data.frame(umap_res)
colnames(umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
umap_scores <- umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
umap_plot <- ggplot(umap_scores, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = expected_module, shape = database),
             size = 6,
             alpha = 0.8) +
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
       filename = "resnik_umap_plot.pdf",
       width = 9,
       height = 7)
