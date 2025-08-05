library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

setwd("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/")

load("all_525_result_with_module_info.rda")
load("all_variable_info.rda")

# Check whether functional modules from the same tissue are clustered in the same module
repeated_in_modules <-
  all_info |>
  group_by(tissue, fm_module) |>
  summarise(repeated_time = n(),
            .groups = 'drop')
sum(repeated_in_modules$repeated_time > 1)

# Calculate the value in the cell
unique_tissues <- unique(all_info$tissue)
all_fm_m <- unique(all_info$fm_module)

heatmap_data <- data.frame()

for (i in unique_tissues) {
  tissue_i_hp_dt <- c()
  tissue_i_info <- all_info |> filter(tissue == i)
  tissue_i_fm_m <- unique(tissue_i_info$fm_module)
  for (m in seq_along(all_fm_m)) {
    fm_m <- all_fm_m[m]
    if (fm_m %in% tissue_i_fm_m) {
      mapped_genes <- tissue_i_info$mapped_molecules[tissue_i_info$fm_module == fm_m]
      mapped_genes <- strsplit(mapped_genes, split = "/") |> unlist()
      tissue_i_mapped_genes_trend_info <- all_variable_info |>
        filter(Tissue == i & ensembl %in% mapped_genes) |>
        group_by(Cluster) |>
        summarise(num = n(), .groups = "drop")
      if (nrow(tissue_i_mapped_genes_trend_info) == 1 & "Cluster D" == tissue_i_mapped_genes_trend_info$Cluster[1]) {
        tissue_i_fm_m_value <- -1
      } else if (nrow(tissue_i_mapped_genes_trend_info) == 1 & "Cluster U" == tissue_i_mapped_genes_trend_info$Cluster[1]) {
        tissue_i_fm_m_value <- 1
      } else if (nrow(tissue_i_mapped_genes_trend_info) == 2) {
        u <- tissue_i_mapped_genes_trend_info$num[tissue_i_mapped_genes_trend_info$Cluster == "Cluster U"]
        d <- tissue_i_mapped_genes_trend_info$num[tissue_i_mapped_genes_trend_info$Cluster == "Cluster D"]
        tissue_i_fm_m_value <- (u - d)/(u + d)
      }
      tissue_i_hp_dt[m] <- tissue_i_fm_m_value
    } else {
      tissue_i_hp_dt[m] <- 0
    }
  }
  names(tissue_i_hp_dt) <- all_fm_m
  heatmap_data <- rbind(heatmap_data, data.frame(t(tissue_i_hp_dt)))
}

rownames(heatmap_data) <- unique_tissues

heatmap_matrix <- as.matrix(heatmap_data)

# Step 4: Create the heatmap
library(ComplexHeatmap)

p <-
  pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = FALSE,
  fontsize = 8,
  fontsize_row = 8,
  fontsize_col = 8,
  row_dend_side = "right",
  column_dend_side = "bottom",
  border_color = "white",
  border = TRUE
  )
p
save(heatmap_matrix, file = "shared_tissues_fm/heatmap_matrix.rda")

pdf("shared_tissues_fm/heatmap_all_56_162.pdf", width = 12, height = 8)
ComplexHeatmap::draw(p)
dev.off()

load("shared_tissues_fm/heatmap_matrix.rda")

# Remove columns with only one non-zero value
load("shared_tissues_fm/heatmap_matrix.rda")
heatmap_matrix_filtered <- heatmap_matrix[, colSums(heatmap_matrix != 0) > 3]
heatmap_matrix_filtered <- heatmap_matrix_filtered[rowSums(heatmap_matrix_filtered != 0) > 2,]

setwd("..")
tissue_sys_map <- readxl::read_xlsx("01_transcriptomics/tissue_system_map.xlsx")
annotation_row <- tissue_sys_map |>
  select(Tissue, System) |>
  filter(Tissue %in% rownames(heatmap_matrix_filtered)) |>
  column_to_rownames("Tissue")
annotation_row <- annotation_row[rownames(heatmap_matrix_filtered), , drop = FALSE]

library(ComplexHeatmap)

# library(paletteer)
# paletteer_d("rcartocolor::Temps")

# colors <- colorRampPalette(c("#CF597EFF", "#E88471FF", "#EEB479FF", "#E9E29CFF",
#             "#9CCB86FF", "#39B185FF", "#009392FF"))(10)

system_colors <- c("Adipose tissue" = "#483732", "Digestive system" = "#064c96", "Endocrine system" = "#acd386",
                   "Immune system" = "#e435a1", "Integumentary system" = "#eac324", "Muscular system" = "#5e2d8c",
                   "Nervous system" = "#10939f", "Reproductive system" = "#e37c7b", "Respiratory system" = "#e66514",
                   "Skeletal system" = "#3a9b48")

# right_anno <- rowAnnotation(
#   "ARI Score" = anno_boxplot(heatmap_matrix,
#                              gp = gpar(col = box_colors)),
#   width = unit(3, "cm")
# )

# left_anno <- HeatmapAnnotation(
#   "Aging related genes" =
# )

p <-
  pheatmap(
    heatmap_matrix_filtered,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    display_numbers = FALSE,
    fontsize = 8,
    fontsize_row = 8,
    fontsize_col = 8,
    row_dend_side = "right",
    column_dend_side = "bottom",
    border_color = "white",
    border = TRUE,
    # annotation_row = annotation_row,
    # annotation_row_side = "right"
    right_annotation = rowAnnotation(
      df = annotation_row,
      col = list(System = system_colors)
    ),
    left_annotation = row_ha,
    top_annotation = col_ha,
    heatmap_legend_param = list(
      title = "Dysregulation index"
    )
  )

lgd_gene <- Legend(
  title = "Gene number",
  labels = c("Cluster U", "Cluster D"),
  legend_gp = gpar(fill = c("#800074", "#298c8c"))
)

lgd_fm <- Legend(
  title = "Tissue number (Dysregulation index)",
  labels = c("< 0", "= 0", "> 0"),
  legend_gp = gpar(fill = c("#5e4c5f", "#999999", "#ffbb6f"))
)

draw(p, annotation_legend_list = list(lgd_gene, lgd_fm), annotation_legend_side = "right")

pdf("01_transcriptomics/shared_tissues_fm/filtered_heatmap.pdf", width = 11, height = 6)
ComplexHeatmap::draw(p)
dev.off()

### Annotation row ====
filtered_unique_tissues <- unique(rownames(heatmap_matrix_filtered))
filtered_fm_m <- unique(colnames(heatmap_matrix_filtered))
filtered_all_info <- all_info |> dplyr::filter(fm_module %in% filtered_fm_m)

# Three columns: Tissue, Number of cluster U gene and cluster D
tissue_gene_dt <- data.frame()

for (i in filtered_unique_tissues) {
  tissue_i_info <- filtered_all_info |> filter(tissue == i)
  mapped_genes <- tissue_i_info |> pull(mapped_molecules) |> strsplit(split = "/") |> unlist() |> unique()
  tissue_i_mapped_genes_trend_info <- all_variable_info |>
    filter(Tissue == i & ensembl %in% mapped_genes) |>
    group_by(Cluster) |>
    summarise(num = n(), .groups = "drop")
  if (nrow(tissue_i_mapped_genes_trend_info) == 2) {
    cluster_u_gene_num <- tissue_i_mapped_genes_trend_info$num[tissue_i_mapped_genes_trend_info$Cluster == "Cluster U"]
    cluster_d_gene_num <- tissue_i_mapped_genes_trend_info$num[tissue_i_mapped_genes_trend_info$Cluster == "Cluster D"]
  } else if (nrow(tissue_i_mapped_genes_trend_info) == 1 & tissue_i_mapped_genes_trend_info$Cluster[1] == "Cluster U") {
    cluster_u_gene_num <- tissue_i_mapped_genes_trend_info$num[1]
    cluster_d_gene_num <- 0
  } else if (nrow(tissue_i_mapped_genes_trend_info) == 1 & tissue_i_mapped_genes_trend_info$Cluster[1] == "Cluster D") {
    cluster_u_gene_num <- 0
    cluster_d_gene_num <- tissue_i_mapped_genes_trend_info$num[1]
  }

  tissue_gene_dt <- rbind(tissue_gene_dt,
                          data.frame("tissue" = i,
                                     "cluster_u" = cluster_u_gene_num,
                                     "cluster_d" = cluster_d_gene_num
                          ))
}

tissue_gene_dt <- tissue_gene_dt |> column_to_rownames("tissue")
save(tissue_gene_dt, file = "01_transcriptomics/shared_tissues_fm/tissue_gene_dt.rda")
load("01_transcriptomics/shared_tissues_fm/tissue_gene_dt.rda")
row_ha <- ComplexHeatmap::rowAnnotation(
  gene_num = anno_barplot(tissue_gene_dt,
                          gp = gpar(fill = c("cluster_u" = "#800074",
                                             "cluster_d" = "#298c8c")),
                          axis_param = list(direction = "reverse"),
                          width = unit(1.5, "cm")),
  show_annotation_name = FALSE
)

## Column annotation ====
# fm_tissue_dt <- data.frame()
#
# for (i in filtered_fm_m) {
#   fm_i_info <- filtered_all_info |> filter(fm_module == i)
#   cluster_gene <- data.frame()
#   for (tissue_i in fm_i_info$tissue) {
#     mapped_genes <- fm_i_info |>
#       dplyr::filter(tissue == tissue_i) |>
#       pull(mapped_molecules) |>
#       strsplit(split = "/") |>
#       unlist()
#     tissue_i_mapped_genes_trend_info <- all_variable_info |>
#       dplyr::filter(Tissue == tissue_i & ensembl %in% mapped_genes) |>
#       dplyr::select(Cluster, ensembl)
#     cluster_gene <- rbind(cluster_gene, tissue_i_mapped_genes_trend_info)
#   }
#   cluster_gene <- cluster_gene |> distinct(Cluster, ensembl)
#   cluster_gene_num <- cluster_gene |>
#     group_by(Cluster) |>
#     summarise(num = n(), .groups = "drop")
#   if (nrow(cluster_gene_num) == 2) {
#     cluster_u_gene_num <- cluster_gene_num$num[cluster_gene_num$Cluster == "Cluster U"]
#     cluster_d_gene_num <- cluster_gene_num$num[cluster_gene_num$Cluster == "Cluster D"]
#   } else if (nrow(cluster_gene_num) == 1 & cluster_gene_num$Cluster[1] == "Cluster U") {
#     cluster_u_gene_num <- cluster_gene_num$num[1]
#     cluster_d_gene_num <- 0
#   } else if (nrow(cluster_gene_num) == 1 & cluster_gene_num$Cluster[1] == "Cluster D") {
#     cluster_u_gene_num <- 0
#     cluster_d_gene_num <- cluster_gene_num$num[1]
#   }
#
#   fm_gene_dt <- rbind(fm_gene_dt,
#                       data.frame("fm" = i,
#                                  "cluster_u" = cluster_u_gene_num,
#                                  "cluster_d" = cluster_d_gene_num
#                           ))
# }
#
# fm_gene_dt <- fm_gene_dt |> column_to_rownames("fm")
# save(fm_gene_dt, file = "01_transcriptomics/shared_tissues_fm/fm_gene_dt.rda")

fm_tissue_num <- data.frame(
  fm = colnames(heatmap_matrix_filtered),
  negative = apply(heatmap_matrix_filtered, 2, function(x) sum(x < 0, na.rm = TRUE)),
  zero = apply(heatmap_matrix_filtered, 2, function(x) sum(x == 0, na.rm = TRUE)),
  positive = apply(heatmap_matrix_filtered, 2, function(x) sum(x > 0, na.rm = TRUE))
)
fm_tissue_num <- fm_tissue_num |> select(-fm)

col_ha <- ComplexHeatmap::HeatmapAnnotation(
  gene_num = anno_barplot(fm_tissue_num,
                          gp = gpar(fill = c("negative" = "#5e4c5f",
                                             "zero" = "#999999",
                                             "positive" = "#ffbb6f")),
                          height = unit(1.5, "cm")),
  show_annotation_name = FALSE
)

