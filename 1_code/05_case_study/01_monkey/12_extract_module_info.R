library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(Cairo)

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/llm_interpreted_res_processing_results_summary.rda")
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_all_singleton_file_info.rda")
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/error_module_not_enough_to_cluster.rda")
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/all_525_result_with_module_info.rda")
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/shared_tissues_fm/heatmap_matrix.rda")

setwd("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics")

heatmap_matrix_filtered <- heatmap_matrix[, colSums(heatmap_matrix != 0) > 3]
heatmap_matrix_filtered <- heatmap_matrix_filtered[rowSums(heatmap_matrix_filtered != 0) > 2,]

heatmap_matrix_filtered |> dim()

ht_fm <- colnames(heatmap_matrix_filtered)
ht_tissues <- rownames(heatmap_matrix_filtered)

fm_file_info <- final_results |>
  dplyr::filter(status == "success") |>
  dplyr::select(tissue, fm_file)

dir.create("shared_tissues_fm/tissue_module_networks")

for (fm in ht_fm) {
  # Get tissue fm id in the meta-fm
  tissue_fm_ids <- all_info |>
    dplyr::filter(fm_module == fm) |>
    dplyr::filter(module_size > 1) |>
    dplyr::select(tissue, fm_node, module_size, fm_module)

  # create module network for each tissue
  if (nrow(tissue_fm_ids) == 0) {
    message("For ", fm, ": It only contains singleton.")
  } else {
    dir_path <- file.path("shared_tissues_fm/tissue_module_networks", fm)
    dir.create(dir_path)

    for (t in tissue_fm_ids$tissue) {
      tissue_fm_id <- tissue_fm_ids$fm_node[tissue_fm_ids$tissue == t]
      fm_id <- str_extract(tissue_fm_id, "Functional_module_\\d+")

      file_path <- fm_file_info$fm_file[fm_file_info$tissue == t]
      load(file_path)
      plot <-
        plot_module_info(
          object = llm_interpreted_fm_res,
          level = "functional_module",
          module_id = fm_id,
          llm_text = TRUE
        )

      network <- plot[[1]]
      plot_file_name <- paste(tissue_fm_id, "network.pdf", sep = "_")
      plot_path <- file.path(dir_path, plot_file_name)
      CairoPDF(plot_path, width = 10, height = 8)
      print(network)
      dev.off()
    }
  }
}
