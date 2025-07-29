library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/summary/llm_interpreted_res_processing_results_summary.rda")

setwd("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/")

# Step1: Collect text information for each functional module ====
# Initialize the combined data frame
all_info <- data.frame()

# TYPE 1: For tissues that have modules and singletons (LLM interpreted results)
llm_res_tissue <- final_results[is.na(final_results$error_message), ]

for (i in 1:nrow(llm_res_tissue)) {
  fm_file_name <- llm_res_tissue$fm_file[i]
  tissue <- llm_res_tissue$tissue[i]
  cluster <- llm_res_tissue$cluster[i]

  file_path <- file.path("llm_interpreted_res", fm_file_name)

  if (!file.exists(file_path)) {
    cat("Warning: File", fm_file_name, "not found. Skipping...\n")
    next
  }

  cat("Processing LLM results:", tissue, "-", cluster, "\n")

  load(file_path)

  object <- llm_interpreted_fm_res

  llm_summary <- object@llm_module_interpretation
  fm_result <- object@merged_module$functional_module_result

  # For FM with module_content_number > 1 (actual modules)
  for (fm_name in names(llm_summary)) {
    fm_info <- fm_result[fm_result$module == fm_name, ]

    if (nrow(fm_info) == 0) next

    llm_interpretation <- sprintf("%s: %s",
                                  llm_summary[[fm_name]]$generated_name$module_name,
                                  llm_summary[[fm_name]]$generated_name$summary)

    collected_info <- data.frame(
      "tissue" = tissue,
      "cluster" = cluster,
      "node" = fm_info$node,
      "mapped_molecules" = fm_info$geneID,
      "mapped_molecule_count" = fm_info$Count,
      "functional_module_name" = paste(c(tissue, cluster, fm_name), collapse = "_"),
      "Annotation" = llm_interpretation,
      stringsAsFactors = FALSE
    )

    all_info <- rbind(all_info, collected_info)
  }

  # For FM with module_content_number = 1 (singletons)
  singleton_module <- fm_result[fm_result$module_content_number == 1, ]

  for (fm_indx in 1:nrow(singleton_module)) {
    if (nrow(singleton_module) == 0) break

    id <- singleton_module$node[fm_indx]
    fm_info <- singleton_module[fm_indx, ]

    pathway_info <- tryCatch({
      if (grepl("^GO:", id)) {
        mapa::get_go_info(go_ids = id)
      } else if (grepl("^R-", id)) {
        mapa::get_reactome_pathway_info(reactome_ids = id)
      } else {
        mapa::get_kegg_pathway_info(kegg_ids = id)
      }
    }, error = function(e) {
      cat("Error getting pathway info for", id, ":", e$message, "\n")
      return(list(list(term_name = "Unknown", term_definition = "No definition available")))
    })

    collected_info <- data.frame(
      "tissue" = tissue,
      "cluster" = cluster,
      "node" = fm_info$node,
      "mapped_molecules" = fm_info$geneID,
      "mapped_molecule_count" = fm_info$Count,
      "functional_module_name" = paste(c(tissue, cluster, fm_info$module), collapse = "_"),
      "Annotation" = sprintf("%s: %s",
                             pathway_info[[1]]$term_name,
                             pathway_info[[1]]$term_definition),
      stringsAsFactors = FALSE
    )

    all_info <- rbind(all_info, collected_info)
  }
}

# TYPE 2: For tissues that only have singletons (no clustering possible)
fm_res_tissue <- final_results[grepl("^module_content_number_cutoff", final_results$error_message), ]

for (i in 1:nrow(fm_res_tissue)) {
  tissue_name <- fm_res_tissue$tissue[i]
  cluster_name <- fm_res_tissue$cluster[i]

  safe_tissue_name <- gsub("[^A-Za-z0-9]", "_", tissue_name)
  safe_cluster_name <- gsub("[^A-Za-z0-9]", "_", cluster_name)

  fm_filename <- paste0("fm_res_", safe_tissue_name, "_", safe_cluster_name, ".rda")
  file_path <- file.path("llm_interpreted_res", fm_filename)

  if (!file.exists(file_path)) {
    cat("Warning: File", fm_filename, "not found. Skipping...\n")
    next
  }

  cat("Processing FM results:", tissue_name, "-", cluster_name, "\n")

  load(file_path)

  object <- functional_module_res
  singleton_module <- object@merged_module$functional_module_result

  for (fm_indx in 1:nrow(singleton_module)) {
    if (nrow(singleton_module) == 0) break

    id <- singleton_module$node[fm_indx]
    fm_info <- singleton_module[fm_indx, ]

    pathway_info <- tryCatch({
      if (grepl("^GO:", id)) {
        mapa::get_go_info(go_ids = id)
      } else if (grepl("^R-", id)) {
        mapa::get_reactome_pathway_info(reactome_ids = id)
      } else {
        mapa::get_kegg_pathway_info(kegg_ids = id)
      }
    }, error = function(e) {
      cat("Error getting pathway info for", id, ":", e$message, "\n")
      return(list(list(term_name = "Unknown", term_definition = "No definition available")))
    })

    collected_info <- data.frame(
      "tissue" = tissue_name,
      "cluster" = cluster_name,
      "node" = fm_info$node,
      "mapped_molecules" = fm_info$geneID,
      "mapped_molecule_count" = fm_info$Count,
      "functional_module_name" = paste(c(tissue_name, cluster_name, fm_info$module), collapse = "_"),
      "Annotation" = sprintf("%s: %s",
                             pathway_info[[1]]$term_name,
                             pathway_info[[1]]$term_definition),
      stringsAsFactors = FALSE
    )

    all_info <- rbind(all_info, collected_info)
  }
}

# TYPE 3: For tissues that only have enriched pathways (cannot perform clustering)
no_clustering_tissue <- final_results[grepl("must have n >= 2 objects to cluster", final_results$error_message), ]

for (i in 1:nrow(no_clustering_tissue)) {
  tissue_name <- no_clustering_tissue$tissue[i]
  cluster_name <- no_clustering_tissue$cluster[i]

  safe_tissue_name <- gsub("[^A-Za-z0-9]", "_", tissue_name)
  safe_cluster_name <- gsub("[^A-Za-z0-9]", "_", cluster_name)
  filename <- paste0("embedsim_res_", safe_tissue_name, "_", safe_cluster_name, ".rda")
  file_path <- file.path("embedsim_result", filename)

  if (!file.exists(file_path)) {
    cat("Warning: File", filename, "not found. Skipping...\n")
    next
  }

  cat("Processing embedsim results:", tissue_name, "-", cluster_name, "\n")

  # Load the embedsim_res object
  load(file_path)

  # Get the single pathway ID from the 1x1 similarity matrix
  id <- rownames(embedsim_res$sim_matrix)[1]

  # Get pathway information and corresponding enrichment data
  pathway_info <- tryCatch({
    if (grepl("^GO:", id)) {
      fm_info <- embedsim_res$enriched_pathway@enrichment_go_result@result[
        embedsim_res$enriched_pathway@enrichment_go_result@result$ID == id, ]
      pathway_data <- mapa::get_go_info(go_ids = id)
      list(info = pathway_data, enrichment = fm_info)
    } else if (grepl("^R-", id)) {
      fm_info <- embedsim_res$enriched_pathway@enrichment_reactome_result@result[
        embedsim_res$enriched_pathway@enrichment_reactome_result@result$ID == id, ]
      pathway_data <- mapa::get_reactome_pathway_info(reactome_ids = id)
      list(info = pathway_data, enrichment = fm_info)
    } else {
      fm_info <- embedsim_res$enriched_pathway@enrichment_kegg_result@result[
        embedsim_res$enriched_pathway@enrichment_kegg_result@result$ID == id, ]
      pathway_data <- mapa::get_kegg_pathway_info(kegg_ids = id)
      list(info = pathway_data, enrichment = fm_info)
    }
  }, error = function(e) {
    cat("Error getting pathway info for", id, ":", e$message, "\n")
    return(list(
      info = list(list(term_name = "Unknown", term_definition = "No definition available")),
      enrichment = data.frame(geneID = "", Count = 0, stringsAsFactors = FALSE)
    ))
  })

  if (nrow(pathway_info$enrichment) == 0) {
    cat("Warning: No enrichment data found for", id, "\n")
    next
  }

  collected_info <- data.frame(
    "tissue" = tissue_name,
    "cluster" = cluster_name,
    "node" = id,
    "mapped_molecules" = pathway_info$enrichment$geneID[1],
    "mapped_molecule_count" = pathway_info$enrichment$Count[1],
    "functional_module_name" = paste(c(tissue_name, cluster_name, id), collapse = "_"),
    "Annotation" = sprintf("%s: %s",
                           pathway_info$info[[1]]$term_name,
                           pathway_info$info[[1]]$term_definition),
    stringsAsFactors = FALSE
  )

  all_info <- rbind(all_info, collected_info)
}

# Final check and summary
cat("\nData collection completed!\n")
cat("Total rows in all_info:", nrow(all_info), "\n")
cat("Unique tissues:", length(unique(all_info$tissue)), "\n")
cat("Unique clusters:", length(unique(all_info$cluster)), "\n")

# save(all_info, file = "fm_sim_result/all_info.rda")

# Step2: Get embedding matrix ====
# Initialize list to store embeddings
embedding_list <- list()
failed_embeddings <- c()

# Process each functional module
for (i in 1:nrow(all_info)) {
  fm_name <- all_info$functional_module_name[i]
  fm_text <- all_info$Annotation[i]

  cat("Processing", i, "of", nrow(all_info), ":", fm_name, "\n")

  # Skip if text is empty or NA
  if (is.na(fm_text) || fm_text == "" || fm_text == ": ") {
    cat("  Warning: Empty or invalid text for", fm_name, ". Skipping...\n")
    failed_embeddings <- c(failed_embeddings, fm_name)
    next
  }

  # Get embeddings with error handling
  tryCatch({
    embeddings <- mapa::get_embedding(chunk = fm_text,
                                      api_provider = "openai",
                                      api_key = api_key,
                                      model_name = "text-embedding-3-small")

    # Store the embedding vector
    embedding_list[[fm_name]] <- embeddings

    cat("  Successfully embedded:", fm_name, "(", length(embeddings), "dimensions )\n")

    Sys.sleep(0.1)

  }, error = function(e) {
    cat("  Error embedding", fm_name, ":", e$message, "\n")
    failed_embeddings <- c(failed_embeddings, fm_name)
  })
}

# Convert list to matrix
cat("\nConverting embeddings to matrix...\n")

if (length(embedding_list) > 0) {
  # Get the number of dimensions from the first successful embedding
  n_dims <- length(embedding_list[[1]])
  cat("Embedding dimensions:", n_dims, "\n")

  # Create matrix with functional module names as row names
  embedding_matrix <- matrix(NA,
                             nrow = length(embedding_list),
                             ncol = n_dims,
                             dimnames = list(names(embedding_list),
                                             paste0("dim_", 1:n_dims)))

  # Fill the matrix
  for (i in 1:length(embedding_list)) {
    fm_name <- names(embedding_list)[i]
    embedding_matrix[i, ] <- embedding_list[[fm_name]]
  }

  cat("Embedding matrix created successfully!\n")
  cat("Matrix dimensions:", nrow(embedding_matrix), "x", ncol(embedding_matrix), "\n")

} else {
  cat("Error: No successful embeddings were generated!\n")
  embedding_matrix <- NULL
}

# Report results
cat("\nEmbedding Summary:\n")
cat("Successful embeddings:", length(embedding_list), "\n")
cat("Failed embeddings:", length(failed_embeddings), "\n")

embedding_matrix |> dim()
# [1] 1119 1536


# Step3: Calculate cosine similarity ====
calculate_cosine_sim <- function(m){
  dot_product <- m %*% t(m)
  norm_product <- sqrt(rowSums(m^2)) %*% t(sqrt(rowSums(m^2)))
  cosine_sim <- dot_product / norm_product
  return(cosine_sim)
}

sim_matrix <- calculate_cosine_sim(embedding_matrix)
setwd(get_project_wd())
save(sim_matrix, file = "2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/sim_matrix.rda")

# Step4: Cluster functional module ====
merge_by_hierarchical <- function(sim_matrix,
                                  hclust.method,
                                  sim.cutoff) {
  cosine_dist <- 1 - sim_matrix
  ## Convert distance matrix to a 'dist' object
  cosine_dist_obj <- as.dist(cosine_dist)
  ## Perform hierarchical clustering
  hc <- hclust(cosine_dist_obj, method = hclust.method)

  clusters <- cutree(hc, h = 1 - sim.cutoff)
  cluster_result <-
    data.frame(node = hc$labels,
               module = paste("Functional_module", as.character(clusters), sep = "_"))

  return(cluster_result)
}

fm_cluster_result <- merge_by_hierarchical(sim_matrix = sim_matrix,
                                           hclust.method = "ward.D2",
                                           sim.cutoff = 0.55)

test <- as.data.frame(table(fm_cluster_result$module))
barplot(test$Freq)
sum(test$Freq == 1)

save(fm_cluster_result, file = "2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/fm_cluster_result.rda")

fm_cluster_result <- fm_cluster_result |> dplyr::rename(fm_node = node, fm_module = module)
all_info <- all_info |> dplyr::rename(fm_node = functional_module_name)
all_info <- all_info |> left_join(fm_cluster_result, by = "fm_node")
save(all_info, file = "2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/all_info.rda")

# Step5: Generate summary for the clustered functional modules ====
load("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/all_info.rda")

format_module_data_for_prompt <- function(cluster_modules) {

  # Check if input is valid
  if (nrow(cluster_modules) == 0) {
    return("No modules found in this cluster.")
  }

  formatted_modules <- c()

  for (i in 1:nrow(cluster_modules)) {
    module <- cluster_modules[i, ]

    # Extract regulation direction from cluster column
    regulation <- ifelse(grepl("Cluster U", module$cluster, ignore.case = TRUE),
                         "upregulated",
                         "downregulated")

    # Format each module
    module_text <- sprintf(
      "**Module %d**: %s\n- Tissue: %s\n- Regulation: %s\n-",
      i,
      module$Annotation,
      module$tissue,
      regulation
    )

    formatted_modules <- c(formatted_modules, module_text)
  }

  # Combine all modules
  result <- paste(formatted_modules, collapse = "\n")
  return(result)
}

generate_cluster_analysis_prompt <- function(cluster_modules, prompt_template) {

  # Format the module data
  module_data <- format_module_data_for_prompt(cluster_modules)

  # Replace placeholders in template
  complete_prompt <- gsub("\\{module_data\\}", module_data, prompt_template)

  return(complete_prompt)
}

prompt_template <- readLines("1_code/05_case_study/cluster_summary_prompt.md", warn = FALSE)
prompt_template <- paste(prompt_template, collapse = "\n")

all_annotation_result <- list()  # Initialize this before the loop

for (fm_module_name in unique(all_info$fm_module)) {
  cluster_modules <- all_info[all_info$fm_module == fm_module_name, ]
  prompt_text <- generate_cluster_analysis_prompt(cluster_modules = cluster_modules,
                                                  prompt_template = prompt_template)
  messages <- list(
    list(role = "system", content = "You are an efficient and insightful assistant to a molecular biologist."),
    list(
      role = "user",
      content = prompt_text
    )
  )

  success <- FALSE
  result <- NULL

  for (attempt in 1:3) {

    cat("Attempt", attempt, "for cluster:", fm_module_name, "\n")

    tryCatch({
      gpt_response <- mapa::gpt_api_call(messages = messages,
                                         api_key = api_key,
                                         model = "gpt-4o-mini-2024-07-18",
                                         api_provider = "openai")

      # Clean response (remove markdown)
      clean_response <- gsub("```json\\s*|```\\s*", "", gpt_response)
      clean_response <- trimws(clean_response)

      # Try to parse JSON
      result <- jsonlite::fromJSON(clean_response)

      # Check if required fields exist
      if (!is.null(result$module_name) && !is.null(result$summary)) {
        cat("✓ Success for cluster:", fm_module_name, "\n")
        success <- TRUE
        break  # Exit retry loop on success
      } else {
        cat("Missing required fields in attempt", attempt, "\n")
      }

    }, error = function(e) {
      cat("Error in attempt", attempt, ":", e$message, "\n")
    })
  }

  # Handle failure case
  if (!success) {
    cat("❌ Failed after 3 attempts for cluster:", fm_module_name, "\n")
    result <- list(
      module_name = paste("Failed_Cluster", fm_module_name),
      summary = paste("Unable to generate summary for cluster", fm_module_name, "after 3 attempts.")
    )
  }

  all_annotation_result[[fm_module_name]] <- list(
    module_name = result$module_name,
    summary = result$summary
  )
}

all_res_annotation_df <- dplyr::bind_rows(all_annotation_result, .id = "fm_module")

result_with_module <- all_info
functional_module_llm_interpreted_result <- all_res_annotation_df

save(result_with_module, file = "2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/result_with_module.rda")
save(functional_module_llm_interpreted_result, file = "2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/fm_sim_result/functional_module_llm_interpreted_result.rda")
