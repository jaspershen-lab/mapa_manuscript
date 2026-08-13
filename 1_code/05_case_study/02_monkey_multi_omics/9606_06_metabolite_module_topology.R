library(r4projects)
library(igraph)

setwd(r4projects::get_project_wd())
rm(list = ls())

data_dir <- file.path(
  "2_data", "case_study", "02_monkey_multi_omics", "taxid_9606"
)
module_file <- file.path(
  data_dir,
  "metabolite_containing_multi_omics_module_feature_counts.csv"
)
output_file <- file.path(
  data_dir,
  paste0(
    "metabolite_containing_multi_omics_module_",
    "metabolite_degree_articulation.csv"
  )
)

tissues <- c(
  "aortic_arch", "kidney", "ovary", "spleen",
  "stomach", "thymus", "thyroid"
)
directions <- c("up", "down")

read_target_modules <- function(path) {
  input_lines <- readLines(path, warn = FALSE)
  header_line <- grep(
    "^Tissue,Direction,Module(,|$)",
    input_lines
  )

  if (length(header_line) != 1) {
    stop("Could not identify a unique table header in: ", path)
  }

  result <- read.csv(
    path,
    skip = header_line - 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # The source CSV has a trailing comma, which creates an unnamed empty column.
  result <- result[
    , !is.na(names(result)) & nzchar(names(result)),
    drop = FALSE
  ]

  required_columns <- c(
    "Tissue", "Direction", "Module", "module_content_number",
    "true_omics_num", "Metabolite"
  )
  missing_columns <- setdiff(required_columns, names(result))
  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in target-module table: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  result$Tissue <- tolower(trimws(result$Tissue))
  result$Direction <- tolower(trimws(result$Direction))
  result$Module <- trimws(result$Module)
  result$target_order <- seq_len(nrow(result))

  result
}

node_info_value <- function(node_info, field, default = NA_character_) {
  if (is.null(node_info) || is.null(node_info[[field]])) {
    return(default)
  }

  value <- node_info[[field]]
  if (length(value) == 0 || is.na(value[[1]])) {
    return(default)
  }

  as.character(value[[1]])
}

target_modules <- read_target_modules(module_file)

if (anyDuplicated(target_modules[c("Tissue", "Direction", "Module")])) {
  stop("Target-module table contains duplicated Tissue/Direction/Module rows")
}

invalid_targets <- with(
  target_modules,
  module_content_number < 3 | true_omics_num < 2 | Metabolite < 1
)
if (any(invalid_targets)) {
  stop(
    "Target-module table contains rows that do not satisfy: ",
    "module_content_number >= 3, true_omics_num >= 2, and Metabolite >= 1"
  )
}

unexpected_tissues <- setdiff(unique(target_modules$Tissue), tissues)
unexpected_directions <- setdiff(unique(target_modules$Direction), directions)
if (length(unexpected_tissues) > 0 || length(unexpected_directions) > 0) {
  stop("Target-module table contains an unexpected tissue or direction")
}

all_results <- list()
result_index <- 0L

for (tissue in tissues) {
  for (direction in directions) {
    rda_file <- file.path(
      data_dir, tissue, direction, "multi_omics_fm.rda"
    )
    if (!file.exists(rda_file)) {
      stop("Missing graph-data file: ", rda_file)
    }

    object_env <- new.env(parent = emptyenv())
    object_names <- load(rda_file, envir = object_env)
    if (!"multi_omics_fm" %in% object_names) {
      stop("Object 'multi_omics_fm' was not found in: ", rda_file)
    }

    multi_omics_fm <- object_env$multi_omics_fm
    graph_data <- multi_omics_fm$graph_data
    module_summary <- multi_omics_fm$functional_module_result

    if (!inherits(graph_data, "igraph")) {
      stop("graph_data is not an igraph object in: ", rda_file)
    }
    if (igraph::is_directed(graph_data)) {
      stop("Expected an undirected graph in: ", rda_file)
    }

    current_targets <- target_modules[
      target_modules$Tissue == tissue &
        target_modules$Direction == direction,
      ,
      drop = FALSE
    ]
    if (nrow(current_targets) == 0) {
      next
    }

    global_degree <- igraph::degree(
      graph_data,
      mode = "all",
      loops = FALSE
    )

    for (i in seq_len(nrow(current_targets))) {
      target <- current_targets[i, , drop = FALSE]
      functional_module <- target$Module[[1]]
      graph_module <- sub(
        "^Functional_module_",
        "Module_",
        functional_module
      )

      summary_row <- module_summary[
        module_summary$module == functional_module,
        ,
        drop = FALSE
      ]
      if (nrow(summary_row) != 1) {
        stop(
          "Expected exactly one functional-module summary for ",
          tissue, "/", direction, "/", functional_module
        )
      }
      if (
        summary_row$module_content_number[[1]] < 3 ||
          summary_row$true_omics_num[[1]] < 2 ||
          !isTRUE(summary_row$include_metabolites[[1]])
      ) {
        stop(
          "Module no longer satisfies the selection criteria: ",
          tissue, "/", direction, "/", functional_module
        )
      }
      if (
        summary_row$module_content_number[[1]] !=
          target$module_content_number[[1]] ||
          summary_row$true_omics_num[[1]] != target$true_omics_num[[1]]
      ) {
        stop(
          "Module metadata does not match the target-module table: ",
          tissue, "/", direction, "/", functional_module
        )
      }

      module_vertices <- which(igraph::V(graph_data)$module == graph_module)
      if (length(module_vertices) != target$module_content_number[[1]]) {
        stop(
          "Module size in graph_data does not match the target-module table: ",
          tissue, "/", direction, "/", functional_module
        )
      }

      module_graph <- igraph::induced_subgraph(
        graph_data,
        vids = module_vertices
      )
      module_degree <- igraph::degree(
        module_graph,
        mode = "all",
        loops = FALSE
      )

      articulation_positions <- as.integer(
        igraph::articulation_points(module_graph)
      )
      is_articulation <- seq_len(igraph::vcount(module_graph)) %in%
        articulation_positions

      metabolite_positions <- which(
        igraph::V(module_graph)$node_type == "metabolite"
      )
      if (length(metabolite_positions) != target$Metabolite[[1]]) {
        stop(
          "Metabolite count in graph_data does not match the target table: ",
          tissue, "/", direction, "/", functional_module
        )
      }

      metabolite_info <- igraph::V(module_graph)$node_info[
        metabolite_positions
      ]
      metabolite_name <- vapply(
        metabolite_info,
        node_info_value,
        character(1),
        field = "cpd_name"
      )
      diff_metric <- suppressWarnings(as.numeric(vapply(
        metabolite_info,
        node_info_value,
        character(1),
        field = "diff_metric"
      )))

      result_index <- result_index + 1L
      all_results[[result_index]] <- data.frame(
        Tissue = tissue,
        Direction = direction,
        Module = functional_module,
        module_content_number = target$module_content_number[[1]],
        true_omics_num = target$true_omics_num[[1]],
        metabolite_id = igraph::V(module_graph)$node_id[
          metabolite_positions
        ],
        metabolite_name = metabolite_name,
        diff_metric = diff_metric,
        degree_within_module = as.numeric(
          module_degree[metabolite_positions]
        ),
        degree_global_graph = as.numeric(
          global_degree[module_vertices[metabolite_positions]]
        ),
        is_articulation_point_within_module = is_articulation[
          metabolite_positions
        ],
        target_order = target$target_order[[1]],
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(all_results) == 0) {
  stop("No target-module metabolite results were generated")
}

metabolite_topology <- do.call(rbind, all_results)
metabolite_topology <- metabolite_topology[
  order(metabolite_topology$target_order, metabolite_topology$metabolite_id),
  ,
  drop = FALSE
]
metabolite_topology$target_order <- NULL
row.names(metabolite_topology) <- NULL

expected_rows <- sum(target_modules$Metabolite)
if (nrow(metabolite_topology) != expected_rows) {
  stop(
    "Expected ", expected_rows, " metabolite rows, generated ",
    nrow(metabolite_topology)
  )
}

observed_keys <- unique(
  metabolite_topology[c("Tissue", "Direction", "Module")]
)
target_keys <- target_modules[c("Tissue", "Direction", "Module")]
if (nrow(observed_keys) != nrow(target_keys)) {
  stop("Not every target module is represented in the output")
}

write.csv(
  metabolite_topology,
  file = output_file,
  row.names = FALSE,
  na = ""
)

message(
  "Wrote ", nrow(metabolite_topology), " metabolite rows from ",
  nrow(target_modules), " modules to: ", output_file
)
