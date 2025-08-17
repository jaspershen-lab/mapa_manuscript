library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

## Run MAPA workflow ====
# load("3_data_analysis/02_control_data/05_benchmarking/enriched_res_control_dt_1.rda")

# dir.create("2_data/reproducibility_curated_dataset_1/")
setwd("2_data/reproducibility_curated_dataset_1/")
# # Get variable info
# enriched_res_control_dt_1@enrichment_go_result <-
#   enriched_res_control_dt_1@enrichment_go_result |>
#   dplyr::mutate(ONTOLOGY = "")
#
# enriched_res_control_dt_1@enrichment_kegg_result <-
#   enriched_res_control_dt_1@enrichment_kegg_result |>
#   dplyr::mutate(category = "") |>
#   dplyr::mutate(subcategory = "")
#
# for (i in 1:nrow(enriched_res_control_dt_1@enrichment_go_result@result)) {
#   gene_id_i <- str_split(enriched_res_control_dt_1@enrichment_go_result@result$geneID[i], "/") |> unlist()
#   if (length(gene_id_i) > 20) {
#     random_gene <- sample(length(gene_id_i), 20)
#     selected_gene_id <- gene_id_i[random_gene]
#     enriched_res_control_dt_1@enrichment_go_result@result$geneID[i] <- paste(selected_gene_id, collapse = "/")
#   }
# }
#
# for (i in 1:nrow(enriched_res_control_dt_1@enrichment_kegg_result@result)) {
#   gene_id_i <- str_split(enriched_res_control_dt_1@enrichment_kegg_result@result$geneID[i], "/") |> unlist()
#   if (length(gene_id_i) > 20) {
#     random_gene <- sample(length(gene_id_i), 20)
#     selected_gene_id <- gene_id_i[random_gene]
#     enriched_res_control_dt_1@enrichment_kegg_result@result$geneID[i] <- paste(selected_gene_id, collapse = "/")
#   }
# }
#
# for (i in 1:nrow(enriched_res_control_dt_1@enrichment_reactome_result@result)) {
#   gene_id_i <- str_split(enriched_res_control_dt_1@enrichment_reactome_result@result$geneID[i], "/") |> unlist()
#   if (length(gene_id_i) > 20) {
#     random_gene <- sample(length(gene_id_i), 20)
#     selected_gene_id <- gene_id_i[random_gene]
#     enriched_res_control_dt_1@enrichment_reactome_result@result$geneID[i] <- paste(selected_gene_id, collapse = "/")
#   }
# }
#
# enriched_go_result <- enriched_res_control_dt_1@enrichment_go_result@result
# enriched_kegg_result <- enriched_res_control_dt_1@enrichment_kegg_result@result
# enriched_reactome_result <- enriched_res_control_dt_1@enrichment_reactome_result@result
#
# variable_info_go <-
#   enriched_go_result |>
#   pull(geneID) |>
#   str_split("/") |>
#   unlist()
# variable_info_kegg <-
#   enriched_kegg_result |>
#   pull(geneID) |>
#   str_split("/") |>
#   unlist()
# variable_info_reactome <-
#   enriched_reactome_result |>
#   pull(geneID) |>
#   str_split("/") |>
#   unlist()
# variable_info <- c(variable_info_go, variable_info_kegg, variable_info_reactome)
# variable_info <- unique(variable_info)
#
# variable_info <- data.frame("variable_info" = variable_info)
# variable_info <- variable_info |>
#   dplyr::rename(entrezid = variable_info)
#
# converted_variable_info <- convert_id(
#   data = variable_info,
#   query_type = "gene",
#   from_id_type = "entrezid",
#   organism = "org.Hs.eg.db"
# )
#
# variable_info <- converted_variable_info
#
# save(variable_info, file = "variable_info.rda")
#
# enriched_res_control_dt_1@variable_info <- variable_info
#
# enriched_res_control_dt_1@process_info$merge_pathways@parameter$count.cutoff.go = -1
# enriched_res_control_dt_1@process_info$merge_pathways@parameter$count.cutoff.kegg = -1
# enriched_res_control_dt_1@process_info$merge_pathways@parameter$count.cutoff.reactome = -1
#
# save(enriched_res_control_dt_1, file = "enriched_res_control_dt_1.rda")

load("enriched_res_control_dt_1.rda")

# bioembed_sim_res <-
#   get_bioembedsim(
#     object = enriched_res_control_dt_1,
#     api_provider = "openai",
#     text_embedding_model = "text-embedding-3-small",
#     api_key = api_key,
#     database = c("go", "kegg", "reactome"),
#     count.cutoff.go = 0,
#     count.cutoff.kegg = 0,
#     count.cutoff.reactome = 0
#   )
#
# fm_result <- get_functional_modules(
#   object = bioembed_sim_res,
#   sim.cutoff = 0.55,
#   cluster_method = "louvain"
# )
#
# llm_interpreted_res <-
#   llm_interpret_module(
#     object = fm_result,
#     module_content_number_cutoff = 1,
#     api_key = api_key,
#     embedding_output_dir = "embedding_output",
#     retmax = 1
#   )

# Loop 15 times
for (i in 1:15) {
  cat("Starting iteration", i, "of 15...\n")

  # Your original code block
  bioembed_sim_res <- get_bioembedsim(
    object = enriched_res_control_dt_1,
    api_provider = "openai",
    text_embedding_model = "text-embedding-3-small",
    api_key = api_key,
    database = c("go", "kegg", "reactome"),
    count.cutoff.go = 0,
    count.cutoff.kegg = 0,
    count.cutoff.reactome = 0
  )

  fm_result <- get_functional_modules(
    object = bioembed_sim_res,
    sim.cutoff = 0.55,
    cluster_method = "louvain"
  )

  llm_interpreted_res <- llm_interpret_module(
    object = fm_result,
    module_content_number_cutoff = 1,
    api_key = api_key,
    embedding_output_dir = "embedding_output",
    retmax = 1
  )

  # Save the result with incremented filename
  filename <- paste0("llm_interpreted_res_", i)
  assign(filename, llm_interpreted_res)
  save(list = filename, file = paste0(filename, ".rda"))

  cat("Completed iteration", i, "- saved as", paste0(filename, ".rda"), "\n")

  # Sleep for 1 minute (60 seconds) unless it's the last iteration
  if (i < 15) {
    cat("Sleeping for 1 minute before next iteration...\n")
    Sys.sleep(60)
  }
}

## Compare the clustering reproducibility ====
library(mclust)

# Function to extract clustering results from a single file
extract_clustering_result <- function(file_path, run_id) {
  # Create a new environment to avoid variable name conflicts
  env <- new.env()
  load(file_path, envir = env)

  var_name <- ls(env)[1]  # Assuming there's only one object in the .rda file

  # Extract clustering result
  clustering_res <- env[[var_name]]@merged_module$result_with_module %>%
    dplyr::select(node, module) %>%
    dplyr::mutate(module = as.numeric(sub("Functional_module_", "", module))) %>%
    dplyr::rename(!!paste0("run_", run_id) := module)

  return(clustering_res)
}

all_clustering_results <- list()

for (i in 1:15) {
  file_path <- paste0("llm_interpreted_res_", i, ".rda")

  if (file.exists(file_path)) {
    cat("Loading file:", file_path, "\n")
    all_clustering_results[[i]] <- extract_clustering_result(file_path, i)
  } else {
    cat("Warning: File", file_path, "not found\n")
  }
}

# Remove NULL entries if any files were missing
all_clustering_results <- all_clustering_results[!sapply(all_clustering_results, is.null)]
num_runs <- length(all_clustering_results)

cat("Successfully loaded", num_runs, "clustering results\n")

# Merge all clustering results by node
merged_results <- all_clustering_results[[1]]
for (i in 2:length(all_clustering_results)) {
  merged_results <- merge(merged_results, all_clustering_results[[i]], by = "node", all = TRUE)
}

# Check for missing values
missing_summary <- merged_results %>%
  select(-node) %>%
  summarise_all(~sum(is.na(.)))
missing_summary

# Calculate pairwise ARI between all runs
run_columns <- paste0("run_", 1:num_runs)
n_pairs <- choose(num_runs, 2)
ari_matrix <- matrix(NA, nrow = num_runs, ncol = num_runs)
rownames(ari_matrix) <- run_columns
colnames(ari_matrix) <- run_columns

# Fill diagonal with 1 (perfect agreement with itself)
diag(ari_matrix) <- 1

# Calculate pairwise ARI
pair_count <- 0

for (i in 1:(num_runs-1)) {
  for (j in (i+1):num_runs) {
    run_i <- run_columns[i]
    run_j <- run_columns[j]

    # Extract clustering assignments for both runs (remove NAs)
    complete_cases <- complete.cases(merged_results[[run_i]], merged_results[[run_j]])

    if (sum(complete_cases) > 0) {
      cluster_i <- merged_results[[run_i]][complete_cases]
      cluster_j <- merged_results[[run_j]][complete_cases]

      # Calculate ARI
      ari_value <- mclust::adjustedRandIndex(cluster_i, cluster_j)
      ari_matrix[i, j] <- ari_value
      ari_matrix[j, i] <- ari_value  # Matrix is symmetric

      pair_count <- pair_count + 1
      if (pair_count %% 10 == 0) {
        cat("Processed", pair_count, "of", n_pairs, "pairs\n")
      }
    } else {
      cat("Warning: No complete cases for runs", i, "and", j, "\n")
    }
  }
}

cat("Completed calculating", pair_count, "pairwise ARI values\n")
sum(ari_matrix == 1)

## Compare the annotation reproducibility ====

# Function to extract functional module data from a single file
extract_functional_modules <- function(file_path, run_id) {
  # Create a new environment to avoid variable name conflicts
  env <- new.env()
  load(file_path, envir = env)

  # Get the variable name
  var_name <- ls(env)[1]
  loaded_data <- env[[var_name]]

  # Extract functional modules with more than 1 node
  fm_with_llm_summary <- loaded_data@merged_module$functional_module_result %>%
    dplyr::filter(module_content_number > 1) %>%
    dplyr::select(module, node)

  # Extract LLM annotations for each module
  module_annotations <- list()

  for (module in fm_with_llm_summary$module) {
    # Extract nodes for this module
    module_nodes <- fm_with_llm_summary$node[fm_with_llm_summary$module == module] |> str_split(";") |> unlist()
    module_nodes_sorted <- sort(module_nodes)

    # Create node set signature (sorted concatenated string)
    node_set_signature <- paste(module_nodes_sorted, collapse = "|")

    # Extract LLM annotation
    if (module %in% names(loaded_data@llm_module_interpretation)) {
      llm_annotation <- paste0(
        loaded_data@llm_module_interpretation[[module]]$generated_name$module_name,
        ". ",
        loaded_data@llm_module_interpretation[[module]]$generated_name$summary
      )
    } else {
      llm_annotation <- NA
      cat("Warning: No LLM annotation found for", module, "in run", run_id, "\n")
    }

    # Store module information
    module_annotations[[length(module_annotations) + 1]] <- data.frame(
      run_id = run_id,
      functional_module_id = module,
      node_set_signature = node_set_signature,
      llm_annotation = llm_annotation,
      stringsAsFactors = FALSE
    )
  }

  # Combine all modules for this run
  if (length(module_annotations) > 0) {
    run_modules <- do.call(rbind, module_annotations)
    return(run_modules)
  } else {
    return(data.frame(
      run_id = character(0),
      functional_module_id = character(0),
      node_count = numeric(0),
      node_set_signature = character(0),
      llm_annotation = character(0),
      stringsAsFactors = FALSE
    ))
  }
}

# Load all functional module data from 15 runs
all_module_data <- list()

for (i in 1:15) {
  file_path <- paste0("llm_interpreted_res_", i, ".rda")

  if (file.exists(file_path)) {
    cat("Processing file:", file_path, "\n")
    all_module_data[[i]] <- extract_functional_modules(file_path, i)
  } else {
    cat("Warning: File", file_path, "not found\n")
  }
}

# Remove NULL entries
all_module_data <- all_module_data[!sapply(all_module_data, is.null)]

# Combine all module data
combined_modules <- do.call(rbind, all_module_data)

load("combined_modules.rda")
pairwise_comparisons <- list()
comparison_count <- 0

for (i in 1:(nrow(combined_modules) - 1)) {
  for (j in (i + 1):nrow(combined_modules)) {
    comparison_count <- comparison_count + 1

    pairwise_comparisons[[comparison_count]] <- data.frame(
      node_set_1 = combined_modules$node_set_signature[i],
      run_id_1 = combined_modules$run_id[i],
      functional_module_id_1 = combined_modules$functional_module_id[i],
      llm_annotation_1 = combined_modules$llm_annotation[i],
      node_set_2 = combined_modules$node_set_signature[j],
      run_id_2 = combined_modules$run_id[j],
      functional_module_id_2 = combined_modules$functional_module_id[j],
      llm_annotation_2 = combined_modules$llm_annotation[j],
      stringsAsFactors = FALSE
    )
  }
}


pairwise_annotation_comparisons <- do.call(rbind, pairwise_comparisons)
pairwise_annotation_comparisons <- pairwise_annotation_comparisons |> mutate(class = if_else(node_set_1 == node_set_2,
                                                                                             "same_fm",
                                                                                             "diff_fm"))
save(pairwise_annotation_comparisons, file = "pairwise_annotation_comparisons.rda")

## Get embedding
# combined_modules <- mutate(combined_modules, annotation_emb = NA)
#
# for (i in 1:nrow(combined_modules)) {
#   llm_anno <- combined_modules$llm_annotation[i]
#   anno_emb <- mapa::get_embedding(chunk = llm_anno,
#                                   api_key = api_key,
#                                   model_name = "text-embedding-3-small",
#                                   api_provider = "openai")
#   list_emb <- list(anno_emb)
#   combined_modules$annotation_emb[i] <- list_emb
# }
#
# anno_fm_map <- combined_modules
#
# save(anno_fm_map, file = "anno_fm_map.rda")

load("anno_fm_map.rda")

cosine_similarity <- function(vec1, vec2) {
  dot_product <- sum(vec1 * vec2)
  norm_vec1 <- sqrt(sum(vec1^2))
  norm_vec2 <- sqrt(sum(vec2^2))

  if (norm_vec1 == 0 | norm_vec2 == 0) {
    return(0)
  }

  return(dot_product / (norm_vec1 * norm_vec2))
}

for (i in 1:nrow(pairwise_annotation_comparisons)) {
  run_id_1 <- pairwise_annotation_comparisons$run_id_1[i]
  fm_id_1 <- pairwise_annotation_comparisons$functional_module_id_1[i]
  anno_fm_1_emb <- anno_fm_map$annotation_emb[anno_fm_map$run_id == run_id_1 &
                                              anno_fm_map$functional_module_id == fm_id_1]

  run_id_2 <- pairwise_annotation_comparisons$run_id_2[i]
  fm_id_2 <- pairwise_annotation_comparisons$functional_module_id_2[i]
  anno_fm_2_emb <- anno_fm_map$annotation_emb[anno_fm_map$run_id == run_id_2 &
                                              anno_fm_map$functional_module_id == fm_id_2]

  pairwise_annotation_comparisons$sim[i] <- cosine_similarity(anno_fm_1_emb[[1]], anno_fm_2_emb[[1]])

}

setwd(get_project_wd())
# dir.create("3_data_analysis/08_reproducibility")
pairwise_annotation_comparisons_sim_res <- pairwise_annotation_comparisons
save(pairwise_annotation_comparisons_sim_res, file = "3_data_analysis/08_reproducibility/pairwise_annotation_comparisons_sim_res.rda")




