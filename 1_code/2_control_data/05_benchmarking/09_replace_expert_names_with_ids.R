## Replace human expert names with the expert IDs defined in expert_name_id.xlsx.

result_dir <- file.path(
  "3_data_analysis", "02_control_data", "05_benchmarking",
  "comparison_result"
)
mapping_file <- file.path(
  "2_data", "human_expert_functional_module_annotation",
  "expert_name_id.xlsx"
)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("This script requires the readxl package.")
}
if (!file.exists(mapping_file)) {
  stop("Mapping workbook not found: ", mapping_file)
}

expert_id_table <- readxl::read_excel(mapping_file)
required_mapping_columns <- c("expert_name", "expert_id")
if (!all(required_mapping_columns %in% names(expert_id_table))) {
  stop(
    "Mapping workbook must contain: ",
    paste(required_mapping_columns, collapse = ", ")
  )
}

expert_id_by_alias <- stats::setNames(
  expert_id_table$expert_id,
  expert_id_table$expert_name
)
expert_id_by_name <- unname(expert_id_by_alias[expert_alias_by_name])
names(expert_id_by_name) <- names(expert_alias_by_name)

if (anyNA(expert_id_by_name)) {
  stop("At least one expert alias was not found in the mapping workbook.")
}

replace_values <- function(values) {
  values <- as.character(values)
  matched <- match(values, names(expert_id_by_name))
  values[!is.na(matched)] <- expert_id_by_name[matched[!is.na(matched)]]
  values
}

replace_method_values <- function(file, object_name) {
  input_env <- new.env(parent = emptyenv())
  load(file, envir = input_env)
  object <- input_env[[object_name]]
  object$method <- replace_values(object$method)
  assign(object_name, object)
  save(list = object_name, file = file)
}

replace_expert_column <- function(file, object_name) {
  input_env <- new.env(parent = emptyenv())
  load(file, envir = input_env)
  object <- input_env[[object_name]]
  if ("expert_name" %in% names(object)) {
    object$expert_name <- replace_values(object$expert_name)
    names(object)[names(object) == "expert_name"] <- "expert_id"
  } else if ("expert_id" %in% names(object)) {
    object$expert_id <- replace_values(object$expert_id)
  } else {
    stop(object_name, " contains neither expert_name nor expert_id.")
  }
  assign(object_name, object)
  save(list = object_name, file = file)
}

replace_method_values(
  file.path(result_dir, "all_result_long_data.rda"),
  "all_result_long_data"
)
replace_method_values(
  file.path(result_dir, "filtered_all_long_data.rda"),
  "filtered_all_long_data"
)
replace_method_values(
  file.path(result_dir, "filtered_all_result_long_data.rda"),
  "filtered_all_result_long_data"
)
replace_expert_column(
  file.path(result_dir, "similarity_results.rda"),
  "similarity_results"
)
replace_expert_column(
  file.path(result_dir, "combined_similarity_results.rda"),
  "combined_similarity_results"
)

message("Expert names replaced with expert IDs in five source result files.")
