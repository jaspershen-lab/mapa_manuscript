library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)
# load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/all_combined_info.rda")
# load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/all_text_info.rda")

# dir.create("3_data_analysis/09_select_model_param")
setwd("3_data_analysis/09_select_model_param/")

# Get embeddings ====
# openai_embedding_matrix <- get_embedding_matrix(
#   text = all_combined_info,
#   api_provider = "openai",
#   text_embedding_model = "text-embedding-3-small",
#   api_key = openai_api_key
# )
# openai_embedding_matrix |> dim()
#
# gemini_embedding_matrix <- get_embedding_matrix(
#   text = all_combined_info,
#   api_provider = "gemini",
#   text_embedding_model = "models/text-embedding-004",
#   api_key = gemini_api_key
# )
# gemini_embedding_matrix |> dim()
#
# siliconflow_embedding_matrix <- get_embedding_matrix(
#   text = all_combined_info,
#   api_provider = "siliconflow",
#   text_embedding_model = "Qwen/Qwen3-Embedding-8B",
#   api_key = siliconflow_api_key
# )
# siliconflow_embedding_matrix |> dim()
#
# save(openai_embedding_matrix, file = "openai_embedding_matrix.rda")
# save(gemini_embedding_matrix, file = "gemini_embedding_matrix.rda")
# save(siliconflow_embedding_matrix, file = "siliconflow_embedding_matrix.rda")

# Get similarities ====
load("openai_embedding_matrix.rda")
load("gemini_embedding_matrix.rda")
load("siliconflow_embedding_matrix.rda")

# openai_sim_matrix <- calculate_cosine_sim(m = openai_embedding_matrix)
# gemini_sim_matrix <- calculate_cosine_sim(m = gemini_embedding_matrix)
# siliconflow_sim_matrix <- calculate_cosine_sim(m = siliconflow_embedding_matrix)
# save(openai_sim_matrix, file = "openai_sim_matrix.rda")
# save(gemini_sim_matrix, file = "gemini_sim_matrix.rda")
# save(siliconflow_sim_matrix, file = "siliconflow_sim_matrix.rda")

load("openai_sim_matrix.rda")
load("gemini_sim_matrix.rda")
load("siliconflow_sim_matrix.rda")

remain_pathways <-
  control_dt %>%
  dplyr::count(expected_module) %>%
  dplyr::filter(n > 3)

remain_pathways <-
  control_dt %>%
  dplyr::filter(expected_module %in% remain_pathways$expected_module) %>%
  pull(id)

# openai_temp_data <-
#   openai_embedding_matrix %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "pathway") %>%
#   dplyr::filter(pathway %in% remain_pathways) %>%
#   tibble::column_to_rownames(var = "pathway")

openai_temp_data <-
  openai_sim_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::select(pathway, all_of(remain_pathways)) %>%
  dplyr::filter(pathway %in% remain_pathways) %>%
  tibble::column_to_rownames(var = "pathway") |>
  as.matrix()
openai_temp_data <- as.data.frame(1 - openai_temp_data)

colnames(openai_temp_data) <- NULL

# UMAP - openai ====
library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
set.seed(123)
openai_umap_res <- umap(
  openai_temp_data,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "cosine",
  verbose = TRUE
)

openai_umap_res <- umap(
  openai_temp_data,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean",
  verbose = TRUE
)

# Create a data frame of UMAP results
openai_umap_scores <- as.data.frame(openai_umap_res)
colnames(openai_umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
openai_umap_scores <- openai_umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- openai_umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
openai_umap_plot <- ggplot(openai_umap_scores, aes(x = UMAP1, y = UMAP2)) +
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
openai_umap_plot

ggsave(openai_umap_plot,
       filename = "openai_biotext_embedding_umap_plot.pdf",
       width = 9,
       height = 7)

# UMAP - gemini ====
gemini_temp_data <-
  gemini_sim_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::select(pathway, all_of(remain_pathways)) %>%
  dplyr::filter(pathway %in% remain_pathways) %>%
  tibble::column_to_rownames(var = "pathway") |>
  as.matrix()
gemini_temp_data <- as.data.frame(1 - gemini_temp_data)

colnames(gemini_temp_data) <- NULL

library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
set.seed(24)
gemini_umap_res <- umap(
  gemini_temp_data,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean",
  verbose = TRUE
)

# Create a data frame of UMAP results
gemini_umap_scores <- as.data.frame(gemini_umap_res)
colnames(gemini_umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
gemini_umap_scores <- gemini_umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- gemini_umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
gemini_umap_plot <- ggplot(gemini_umap_scores, aes(x = UMAP1, y = UMAP2)) +
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
gemini_umap_plot

ggsave(umap_plot,
       filename = "gemini_biotext_embedding_umap_plot.pdf",
       width = 9,
       height = 7)

# UMAP - Qwen ====
qwen_temp_data <-
  siliconflow_sim_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "pathway") %>%
  dplyr::select(pathway, all_of(remain_pathways)) %>%
  dplyr::filter(pathway %in% remain_pathways) %>%
  tibble::column_to_rownames(var = "pathway") |>
  as.matrix()

qwen_temp_data <- as.data.frame(1 - qwen_temp_data)

colnames(qwen_temp_data) <- NULL

library(uwot)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

# Run UMAP (rows = pathways)
set.seed(42)
qwen_umap_res <- umap(
  qwen_temp_data,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean",
  verbose = TRUE
)

# Create a data frame of UMAP results
qwen_umap_scores <- as.data.frame(qwen_umap_res)

colnames(qwen_umap_scores) <- c("UMAP1", "UMAP2")

# Add pathway IDs
qwen_umap_scores <- qwen_umap_scores %>%
  tibble::rownames_to_column(var = "pathway") %>%
  left_join(control_dt[, c("id", "expected_module", "database")], by = c("pathway" = "id")) %>%
  mutate(
    expected_module = factor(expected_module, levels = str_sort(
      unique(control_dt$expected_module), numeric = TRUE
    )),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  )

# Compute module label positions (centroids)
module_labels <- qwen_umap_scores %>%
  group_by(expected_module) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with ellipses and module labels
qwen_umap_plot <- ggplot(qwen_umap_scores, aes(x = UMAP1, y = UMAP2)) +
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
qwen_umap_plot

ggsave(umap_plot,
       filename = "qwen_biotext_embedding_umap_plot.pdf",
       width = 9,
       height = 7)


# comparison ====
openai_umap_scores <- openai_umap_scores |> dplyr::mutate(api_provider = "openai")
gemini_umap_scores <- gemini_umap_scores |> dplyr::mutate(api_provider = "gemini")
qwen_umap_scores <- qwen_umap_scores |> dplyr::mutate(api_provider = "siliconflow")

total_plot_data <- rbind(openai_umap_scores, gemini_umap_scores, qwen_umap_scores)

compare_plot <- ggplot(total_plot_data, aes(x = UMAP1, y = UMAP2)) +
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
  theme(legend.position = "right") +
  facet_wrap(~ api_provider)

compare_plot

ggsave(compare_plot,
       filename = "compare_biotext_embedding_umap_plot.pdf",
       width = 22,
       height = 7)

