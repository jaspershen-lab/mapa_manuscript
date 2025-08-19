library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(Cairo)

all_info <- readxl::read_xlsx("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/all_info.xlsx")
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/shared_tissues_fm/heatmap_matrix.rda")
# dir.create("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq")

## figure 6b ----
# brain_fl_meta_fm_ids <- names(heatmap_matrix_filtered["Brain (FL)",])[heatmap_matrix_filtered["Brain (FL)",] != 0]
# brain_pl_meta_fm_ids <- names(heatmap_matrix_filtered["Brain (PL)",])[heatmap_matrix_filtered["Brain (PL)",] != 0]
# common_meta_fm <- intersect(brain_fl_meta_fm_ids, brain_pl_meta_fm_ids)

brain_fl_meta_fm_ids <- names(heatmap_matrix["Brain (FL)",])[heatmap_matrix["Brain (FL)",] != 0]
brain_pl_meta_fm_ids <- names(heatmap_matrix["Brain (PL)",])[heatmap_matrix["Brain (PL)",] != 0]
common_meta_fm <- intersect(brain_fl_meta_fm_ids, brain_pl_meta_fm_ids)


brain_fl_fm_ids <- all_info |>
  filter(tissue == "Brain (FL)" & fm_module %in% common_meta_fm) |>
  pull(fm_node)
brain_fl_fm_ids <- gsub("Brain \\(FL\\)_", "", brain_fl_fm_ids)

brain_pl_fm_ids <- all_info |>
  filter(tissue == "Brain (PL)" & fm_module %in% common_meta_fm) |>
  pull(fm_node)
brain_pl_fm_ids <- gsub("Brain \\(PL\\)_", "", brain_pl_fm_ids)

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Brain__FL_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = brain_fl_fm_ids,
  llm_text = TRUE,
  text_all = TRUE
  )
plot

plot_data_dys_index <- heatmap_matrix[c("Brain (FL)"), common_meta_fm]
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

brain_fl_fm <- all_info |>
  filter(tissue == "Brain (FL)" & fm_module %in% common_meta_fm) |>
  select(fm_node, node, fm_module) |>
  rename(meta_fm = fm_module) |>
  left_join(plot_data_dys_index)

graph_data <- graph_data |>
  tidygraph::activate(nodes) |>
  mutate(
    meta_module = sapply(node, function(x) brain_fl_fm$meta_fm[grepl(x, brain_fl_fm$node)]),
    dys_index = sapply(node, function(x) brain_fl_fm$plot_data_dys_index[grepl(x, brain_fl_fm$node)])
  )

plot <-
  ggraph::ggraph(lay) +
  ggraph::geom_edge_link(
    aes(width = sim),
    color = "grey",
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = dys_index,
        size = if(analysis_type == "enrich_pathway") -log(p_adjust, 10) else abs(NES)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  # guides(fill = guide_legend(ncol = 1)) +
  # scale_fill_manual(values = colors) +
  scale_fill_gradientn(colors = c("#4576b6", "#f8f8c7", "#d93127"),
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       limits = c(-1, 1)) +
  ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(4, 10)) +
  labs(size = if(analysis_type == "enrich_pathway") "-log10(FDR adjusted P-values)" else "abs(NES)") +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot

plot <-
  plot +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = stringr::str_wrap(label, 30)),
                         size = 3,
                         repel = TRUE)

plot1 <- plot

plot1

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/figure6b_brain_fl_network_shared_with_brain_pl.pdf", width = 8, height = 5)
print(plot1)
dev.off()

## figure 6c -----
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Brain__PL_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = brain_pl_fm_ids,
  llm_text = TRUE
)
plot

plot_data_dys_index <- heatmap_matrix[c("Brain (PL)"), common_meta_fm]
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

brain_pl_fm <- all_info |>
  filter(tissue == "Brain (PL)" & fm_module %in% common_meta_fm) |>
  select(fm_node, node, fm_module) |>
  rename(meta_fm = fm_module) |>
  left_join(plot_data_dys_index)

graph_data <- graph_data |>
  tidygraph::activate(nodes) |>
  mutate(
    meta_module = sapply(node, function(x) brain_pl_fm$meta_fm[grepl(x, brain_pl_fm$node)]),
    dys_index = sapply(node, function(x) brain_pl_fm$plot_data_dys_index[grepl(x, brain_pl_fm$node)])
  )

plot <-
  ggraph::ggraph(lay) +
  ggraph::geom_edge_link(
    aes(width = sim),
    color = "grey",
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = dys_index,
        size = if(analysis_type == "enrich_pathway") -log(p_adjust, 10) else abs(NES)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  # guides(fill = guide_legend(ncol = 1)) +
  # scale_fill_manual(values = colors) +
  scale_fill_gradientn(colors = c("#4576b6", "#f8f8c7", "#d93127"),
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       limits = c(-1, 1)) +
  ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(4, 10)) +
  labs(size = if(analysis_type == "enrich_pathway") "-log10(FDR adjusted P-values)" else "abs(NES)") +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot

plot <-
  plot +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = stringr::str_wrap(label, 30)),
                         size = 3,
                         repel = TRUE)

plot2 <- plot
plot2

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/figure6c_brain_pl_network_shared_with_brain_fl.pdf", width = 8, height = 5)
print(plot2)
dev.off()

## figure 6d ----
brain_pl_meta_fm_ids <- names(heatmap_matrix["Brain (PL)",])[heatmap_matrix["Brain (PL)",] != 0]
int_duodenum_meta_fm_ids <- names(heatmap_matrix["Intestine (Duodenum)",])[heatmap_matrix["Intestine (Duodenum)",] != 0]
common_meta_fm <- intersect(brain_pl_meta_fm_ids, int_duodenum_meta_fm_ids)

brain_pl_fm_ids <- all_info |>
  filter(tissue == "Brain (PL)" & fm_module %in% common_meta_fm) |>
  pull(fm_node)
brain_pl_fm_ids <- gsub("Brain \\(PL\\)_", "", brain_pl_fm_ids)

int_duodenum_fm_ids <- all_info |>
  filter(tissue == "Intestine (Duodenum)" & fm_module %in% common_meta_fm) |>
  pull(fm_node)
int_duodenum_fm_ids <- gsub("Intestine \\(Duodenum\\)_", "", int_duodenum_fm_ids)

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Brain__PL_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = brain_pl_fm_ids,
  llm_text = TRUE,
  text_all = TRUE
)
plot

plot_data_dys_index <- heatmap_matrix[c("Brain (PL)"), common_meta_fm]
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

brain_pl_fm <- all_info |>
  filter(tissue == "Brain (PL)" & fm_module %in% common_meta_fm) |>
  select(fm_node, node, fm_module) |>
  rename(meta_fm = fm_module) |>
  left_join(plot_data_dys_index)

graph_data <- graph_data |>
  tidygraph::activate(nodes) |>
  mutate(
    meta_module = sapply(node, function(x) brain_pl_fm$meta_fm[grepl(x, brain_pl_fm$node)]),
    dys_index = sapply(node, function(x) brain_pl_fm$plot_data_dys_index[grepl(x, brain_pl_fm$node)])
  )

plot <-
  ggraph::ggraph(lay) +
  ggraph::geom_edge_link(
    aes(width = sim),
    color = "grey",
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = dys_index,
        size = if(analysis_type == "enrich_pathway") -log(p_adjust, 10) else abs(NES)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  # guides(fill = guide_legend(ncol = 1)) +
  # scale_fill_manual(values = colors) +
  scale_fill_gradientn(colors = c("#4576b6", "#f8f8c7", "#d93127"),
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       limits = c(-1, 1)) +
  ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(4, 10)) +
  labs(size = if(analysis_type == "enrich_pathway") "-log10(FDR adjusted P-values)" else "abs(NES)") +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot

plot <-
  plot +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = stringr::str_wrap(label, 30)),
                         size = 3,
                         repel = TRUE)

plot3 <- plot
plot3

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/figure6d_brain_pl_network_shared_with_int_duodenum.pdf", width = 8, height = 5)
print(plot3)
dev.off()

## figure 6e ----
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Intestine__Duodenum_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = int_duodenum_fm_ids,
  llm_text = TRUE,
  text_all = TRUE
)
plot

plot_data_dys_index <- heatmap_matrix[c("Intestine (Duodenum)"), common_meta_fm]
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

int_dm_fm <- all_info |>
  filter(tissue == "Intestine (Duodenum)" & fm_module %in% common_meta_fm) |>
  select(fm_node, node, fm_module) |>
  rename(meta_fm = fm_module) |>
  left_join(plot_data_dys_index)

graph_data <- graph_data |>
  tidygraph::activate(nodes) |>
  mutate(
    meta_module = sapply(node, function(x) int_dm_fm$meta_fm[grepl(x, int_dm_fm$node)]),
    dys_index = sapply(node, function(x) int_dm_fm$plot_data_dys_index[grepl(x, int_dm_fm$node)])
  )

plot <-
  ggraph::ggraph(lay) +
  ggraph::geom_edge_link(
    aes(width = sim),
    color = "grey",
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = dys_index,
        size = if(analysis_type == "enrich_pathway") -log(p_adjust, 10) else abs(NES)),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  # guides(fill = guide_legend(ncol = 1)) +
  # scale_fill_manual(values = colors) +
  scale_fill_gradientn(colors = c("#4576b6", "#f8f8c7", "#d93127"),
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       limits = c(-1, 1)) +
  ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(4, 10)) +
  labs(size = if(analysis_type == "enrich_pathway") "-log10(FDR adjusted P-values)" else "abs(NES)") +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.background = element_rect(fill = "transparent", color = NA)
  )
plot

plot <-
  plot +
  ggraph::geom_node_text(aes(x = x,
                             y = y,
                             label = stringr::str_wrap(label, 30)),
                         size = 3,
                         repel = TRUE)

plot4 <- plot
plot4

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/figure6e_int_duodenum_network_shared_with_brain_pl.pdf", width = 8, height = 5)
print(plot4)
dev.off()

##supplementary ====
brain_fl_meta_fm_ids <- names(heatmap_matrix_filtered["Brain (FL)",])[heatmap_matrix_filtered["Brain (FL)",] != 0]
int_duodenum_meta_fm_ids <- names(heatmap_matrix_filtered["Intestine (Duodenum)",])[heatmap_matrix_filtered["Intestine (Duodenum)",] != 0]
common_meta_fm <- c(intersect(brain_fl_meta_fm_ids, int_duodenum_meta_fm_ids), "Functional_module_132")

brain_fl_fm_ids <- all_info |>
  filter(tissue == "Brain (FL)" & fm_module %in% common_meta_fm) |>
  pull(fm_node)
brain_fl_fm_ids <- gsub("Brain \\(FL\\)_", "", brain_fl_fm_ids)

int_duodenum_fm_ids <- all_info |>
  filter(tissue == "Intestine (Duodenum)" & fm_module %in% common_meta_fm) |>
  pull(fm_node)
int_duodenum_fm_ids <- gsub("Intestine \\(Duodenum\\)_", "", int_duodenum_fm_ids)

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Brain__FL_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = brain_fl_fm_ids,
  llm_text = TRUE
)
plot

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/6d_brain_fl_fm_network.pdf", width = 5.5, height = 4)
print(plot)
dev.off()

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Intestine__Duodenum_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = int_duodenum_fm_ids,
  llm_text = TRUE
)
plot

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/6d_int_fm_network.pdf", width = 5.5, height = 4)
print(plot)
dev.off()

## FOXO signaling and longevity regulation ----
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Bronchus__Lower_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = "Functional_module_6",
  llm_text = TRUE
)
plot

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/supplementary_2_foxo_Bronchus_Lower.pdf", width = 5.5, height = 4)
print(plot)
dev.off()



