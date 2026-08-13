project_root <- normalizePath(".", mustWork = TRUE)

if (!file.exists(file.path(project_root, "mapa_manuscript.Rproj"))) {
  stop("Run the scripts from the mapa_manuscript project root.")
}

code_dir <- file.path(
  project_root, "1_code", "2_control_data", "07_alt_llm_model_test"
)
output_dir <- file.path(
  project_root, "3_data_analysis", "02_control_data", "07_alt_llm_model_test"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

control_data_file <- file.path(project_root, "2_data", "control_data.xlsx")
pathway_text_file <- file.path(
  project_root, "3_data_analysis", "02_control_data",
  "02_bioembedsim_vs_othersim", "biotext_embedding", "all_combined_info.rda"
)
original_embedding_file <- file.path(
  project_root, "3_data_analysis", "02_control_data",
  "02_bioembedsim_vs_othersim", "biotext_embedding", "embedding_matrix.rda"
)
qwen_embedding_cache_file <- file.path(
  project_root, "3_data_analysis", "09_select_model_param",
  "siliconflow_embedding_matrix.rda"
)
gpt_annotation_file <- file.path(
  project_root, "3_data_analysis", "4_random_module_dataset", "final_result.Rdata"
)
expert_embedding_file <- file.path(
  project_root, "3_data_analysis", "02_control_data", "05_benchmarking",
  "comparison_result", "all_result_long_data.rda"
)

qwen_embedding_model <- "Qwen/Qwen3-Embedding-8B"
qwen_llm_model <- "Qwen/Qwen3.6-35B-A3B"
annotation_embedding_model <- "text-embedding-3-small"

# MAPA 3.0.0 defaults from get_functional_modules.list().
mapa_similarity_cutoff <- 0.55
mapa_clustering_method <- "louvain"
analysis_seed <- 20260807L

first_nonempty_env <- function(names, default = "") {
  values <- Sys.getenv(names, unset = "")
  hit <- which(nzchar(values))[1L]
  if (is.na(hit)) default else unname(values[hit])
}
