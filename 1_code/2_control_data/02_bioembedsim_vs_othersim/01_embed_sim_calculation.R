library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

setwd("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding")

# 1. Get annotated gene information for pathways ====
## 22 GO terms
# go_info <-
#   control_dt |>
#   dplyr::filter(database == "GO") |>
#   dplyr::pull(id) |>
#   get_go_info()
#
# ## 14 KEGG pathways
# kegg_info <-
#   control_dt |>
#   dplyr::filter(database == "KEGG") |>
#   dplyr::pull(id) |>
#   get_kegg_pathway_info()
#
# ## 8 Ractome pathways
# reactome_info <-
#   control_dt |>
#   dplyr::filter(database == "Reactome") |>
#   dplyr::pull(id) |>
#   get_reactome_pathway_info()
#
# all_text_info <- c(go_info, kegg_info, reactome_info)
# all_combined_info <- combine_info(info = all_text_info)
#
# save(all_combined_info, file = "all_combined_info.rda")
# save(all_text_info, file = "all_text_info.rda")

load("all_combined_info.rda")
load("all_text_info.rda")

# 2. Get semantic similarity using embedding model ====
# embedding_matrix <- get_embedding_matrix(
#   text = all_combined_info,
#   api_provider = "openai",
#   text_embedding_model = "text-embedding-3-small",
#   api_key = api_key
# )
#
# save(embedding_matrix, file = "embedding_matrix.rda")
load("embedding_matrix.rda")

embedding_sim_matrix <- calculate_cosine_sim(m = embedding_matrix)
save(embedding_sim_matrix, file = "embedding_sim_matrix.rda")

embedding_sim_df <-
  as.data.frame.table(embedding_sim_matrix, responseName = "sim") |>
  dplyr::filter(Var1 != Var2) |>                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) |>
  dplyr::mutate(across(c(from, to), as.character)) |>
  dplyr::filter(from < to)

save(embedding_sim_df, file = "embedding_sim_df.rda")

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
  embedding_sim_matrix %>%
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
       filename = "biotext_embedding_pca_plot.pdf",
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
  n_neighbors = 15,
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
       filename = "biotext_embedding_umap_plot.pdf",
       width = 9,
       height = 7)
