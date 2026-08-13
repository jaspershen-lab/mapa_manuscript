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

## Combined figure 6b-e with shared legend ----
library(cowplot)

# Fixed limits so all 4 panels use the same data-to-visual mapping.
# Adjust limits = c(1, 4) if your abs(NES) range differs; sim is 0-1.
shared_scale_size <- scale_size_continuous(
  range  = c(2, 10),
  limits = c(0, 10),
  breaks = c(2,4,6,8)
)
shared_scale_edge <- ggraph::scale_edge_width_continuous(
  range  = c(0.1, 1),
  limits = c(0, 1),
  breaks = c(0.25, 0.5, 0.75, 1.0)
)

p_6b <- plot1 + shared_scale_size + shared_scale_edge + theme(legend.position = "none")
p_6c <- plot2 + shared_scale_size + shared_scale_edge + theme(legend.position = "none")
p_6d <- plot3 + shared_scale_size + shared_scale_edge + theme(legend.position = "none")
p_6e <- plot4 + shared_scale_size + shared_scale_edge + theme(legend.position = "none")

shared_legend <- cowplot::get_legend(
  plot1 + shared_scale_size + shared_scale_edge + theme(legend.position = "right")
)

panels_6bcde <- cowplot::plot_grid(
  p_6b, p_6c, p_6d, p_6e,
  ncol       = 2,
  labels     = c("b", "c", "d", "e"),
  label_size = 8
)

figure_6bcde <- cowplot::plot_grid(
  panels_6bcde, shared_legend,
  rel_widths = c(1, 0.25),
  nrow       = 1
)

figure_6bcde

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/figure6bcde_combined.pdf",
         width = 13, height = 8)
print(figure_6bcde)
dev.off()

## supplementary: brain fl and duodenum ====
brain_fl_meta_fm_ids <- names(heatmap_matrix["Brain (FL)",])[heatmap_matrix["Brain (FL)",] != 0]
int_duodenum_meta_fm_ids <- names(heatmap_matrix["Intestine (Duodenum)",])[heatmap_matrix["Intestine (Duodenum)",] != 0]
common_meta_fm <- intersect(brain_fl_meta_fm_ids, int_duodenum_meta_fm_ids)

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
  llm_text = TRUE,
  text_all = TRUE
)
plot

plot_data_dys_index <- heatmap_matrix[c("Brain (FL)"), common_meta_fm]
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

int_dm_fm <- all_info |>
  filter(tissue == "Brain (FL)" & fm_module %in% common_meta_fm) |>
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

plot5 <- plot
plot5

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/supp_fig_1_brain_fl_shared_with_duodenum.pdf", width = 8, height = 5)
print(plot5)
dev.off()

## ----
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Intestine__Duodenum_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = int_duodenum_fm_ids,
  llm_text = TRUE
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

plot6 <- plot
plot6

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/supp_fig_1_duodenum_shared_with_brain_fl.pdf", width = 8, height = 5)
print(plot6)
dev.off()

## FOXO signaling and longevity regulation ----
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Liver__ML_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = "Functional_module_12",
  llm_text = TRUE
)
plot

plot_data_dys_index <- heatmap_matrix[c("Liver (ML)"), c("Functional_module_22")] |> setNames("Functional_module_22")
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

int_dm_fm <- all_info |>
  filter(tissue == "Liver (ML)" & fm_module %in% "Functional_module_22") |>
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

plot7 <- plot
plot7

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/supp_fig_2_foxo_liver_ML.pdf", width = 8, height = 5)
print(plot7)
dev.off()

##----
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Supraspinatus_tendon.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = "Functional_module_3",
  llm_text = TRUE
)
plot

plot_data_dys_index <- heatmap_matrix[c("Supraspinatus tendon"), c("Functional_module_22")] |> setNames("Functional_module_22")
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

int_dm_fm <- all_info |>
  filter(tissue == "Supraspinatus tendon" & fm_module %in% "Functional_module_22") |>
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

plot8 <- plot
plot8

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/supp_fig_2_foxo_Intestine_Ileum.pdf", width = 8, height = 5)
print(plot8)
dev.off()

##----
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Intestine__Ileum_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = "Functional_module_1",
  llm_text = TRUE
)
plot

plot_data_dys_index <- heatmap_matrix[c("Intestine (Ileum)"), c("Functional_module_22")] |> setNames("Functional_module_22")
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

int_dm_fm <- all_info |>
  filter(tissue == "Intestine (Ileum)" & fm_module %in% "Functional_module_22") |>
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

plot9 <- plot
plot9

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/supp_fig_2_foxo_Supraspinatus_tendon.pdf", width = 8, height = 5)
print(plot9)
dev.off()

##----
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_fm_res_object/llm_interpreted_fm_res_Bronchus__Lower_.rda")

plot <- plot_similarity_network(
  object = llm_interpreted_fm_res,
  database = c("go", "kegg"),
  level = "functional_module",
  module_id = "Functional_module_6",
  llm_text = TRUE
)
plot

plot_data_dys_index <- heatmap_matrix[c("Bronchus (Lower)"), c("Functional_module_22")] |> setNames("Functional_module_22")
plot_data_dys_index <- plot_data_dys_index |>
  as.data.frame() |>
  rownames_to_column(var = "meta_fm")

int_dm_fm <- all_info |>
  filter(tissue == "Bronchus (Lower)" & fm_module %in% "Functional_module_22") |>
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

plot10 <- plot
plot10

CairoPDF("3_data_analysis/05_case_study/01_monkey/bulk_rna_seq/supp_fig_2_foxo_Bronchus_Lower.pdf", width = 6, height = 5)
print(plot10)
dev.off()

##
load("3_data_analysis/05_case_study/01_monkey/liver_sn_rna_seq/liver_up_llm_interpreted_fm_res.rda")
mapa::export_functional_module(llm_interpreted_fm_res, path = "3_data_analysis/05_case_study/01_monkey/liver_sn_rna_seq/liver_cluster_u_results")

load("3_data_analysis/05_case_study/01_monkey/brain_sn_rna_seq/brain_down_llm_interpreted_fm_res.rda")
mapa::export_functional_module(llm_interpreted_fm_res, path = "3_data_analysis/05_case_study/01_monkey/brain_sn_rna_seq/brain_cluster_d_results")
load("3_data_analysis/05_case_study/01_monkey/brain_sn_rna_seq/brain_up_llm_interpreted_fm_res.rda")
mapa::export_functional_module(llm_interpreted_fm_res, path = "3_data_analysis/05_case_study/01_monkey/brain_sn_rna_seq/brain_cluster_u_results")

load("3_data_analysis/05_case_study/01_monkey/plasma_proteomics/02_proteomics_down_llm_interpreted_fm_res.rda")
mapa::export_functional_module(llm_interpreted_fm_res, path = "3_data_analysis/05_case_study/01_monkey/plasma_proteomics/02_proteomics_down_results")

export(hmdb_kegg_matched, file = "2_data/curated_pathway_data_2/hmdb_kegg_matched.xlsx")
export(hmdb_reactome_matched, file = "2_data/curated_pathway_data_2/hmdb_reactome_matched.xlsx")
export(kegg_reactome_matched, file = "2_data/curated_pathway_data_2/kegg_reactome_matched.xlsx")
