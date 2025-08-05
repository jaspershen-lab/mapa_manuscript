library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

setwd("2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/")
load("llm_interpreted_res_processing_results_summary.rda")
load("fm_all_singleton_file_info.rda")
load("error_module_not_enough_to_cluster.rda")

# Step1: Collect text information for each functional module ====
# Initialize the combined data frame
all_info <- data.frame()

# TYPE 1: For tissues that have modules and singletons (LLM interpreted results)
llm_res_tissue <- final_results[is.na(final_results$error_message), ]

for (i in 1:nrow(llm_res_tissue)) {
  file_path <- llm_res_tissue$fm_file[i]
  tissue <- llm_res_tissue$tissue[i]

  if (!file.exists(file_path)) {
    cat("Warning: File", file_path, "not found. Skipping...\n")
    next
  }

  cat("Processing LLM results:", tissue, "\n")

  load(file_path)

  object <- llm_interpreted_fm_res

  llm_summary <- object@llm_module_interpretation
  fm_result <- object@merged_module$functional_module_result

  # For FM with module_content_number > 1 (actual modules)
  for (fm_name in names(llm_summary)) {
    fm_info <- fm_result[fm_result$module == fm_name, ]

    if (nrow(fm_info) == 0) next

    llm_interpretation <- sprintf("%s. %s",
                                  llm_summary[[fm_name]]$generated_name$module_name,
                                  llm_summary[[fm_name]]$generated_name$summary)

    collected_info <- data.frame(
      "tissue" = tissue,
      "node" = fm_info$node,
      "mapped_molecules" = fm_info$geneID,
      "mapped_molecule_count" = fm_info$Count,
      "functional_module_name" = paste(c(tissue, fm_name), collapse = "_"),
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
      "node" = fm_info$node,
      "mapped_molecules" = fm_info$geneID,
      "mapped_molecule_count" = fm_info$Count,
      "functional_module_name" = paste(c(tissue, fm_info$module), collapse = "_"),
      "Annotation" = sprintf("%s: %s",
                             pathway_info[[1]]$term_name,
                             pathway_info[[1]]$term_definition),
      stringsAsFactors = FALSE
    )

    all_info <- rbind(all_info, collected_info)
  }
}

# TYPE 2: For tissues that only have singletons (no clustering possible)
fm_res_tissue <- final_results2

for (i in 1:nrow(fm_res_tissue)) {
  tissue_name <- fm_res_tissue$tissue[i]
  file_path <- fm_res_tissue$fm_file[i]

  if (!file.exists(file_path)) {
    cat("Warning: File", file_path, "not found. Skipping...\n")
    next
  }

  cat("Processing FM results:", tissue_name, "\n")

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
      "node" = fm_info$node,
      "mapped_molecules" = fm_info$geneID,
      "mapped_molecule_count" = fm_info$Count,
      "functional_module_name" = paste(c(tissue_name, id), collapse = "_"),
      "Annotation" = sprintf("%s. %s",
                             pathway_info[[1]]$term_name,
                             pathway_info[[1]]$term_definition),
      stringsAsFactors = FALSE
    )

    all_info <- rbind(all_info, collected_info)
  }
}

# TYPE 3: For tissues that only have enriched pathways (cannot perform clustering)
load("successful_biotext_sim_results_file_info_2.rda")
error_module_not_enough_to_cluster <-
  successful_biotext_sim_results_file_info_2 |>
  dplyr::filter(!(Tissue %in% final_results2$tissue))
save(error_module_not_enough_to_cluster, file = "error_module_not_enough_to_cluster.rda")
no_clustering_tissue <- error_module_not_enough_to_cluster

for (i in 1:nrow(no_clustering_tissue)) {
  tissue_name <- no_clustering_tissue$Tissue[i]
  file_path <- no_clustering_tissue$filename[i]

  if (!file.exists(file_path)) {
    cat("Warning: File", file_path, "not found. Skipping...\n")
    next
  }

  cat("Processing embedsim results:", tissue_name, "\n")

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
    "node" = id,
    "mapped_molecules" = pathway_info$enrichment$geneID[1],
    "mapped_molecule_count" = pathway_info$enrichment$Count[1],
    "functional_module_name" = paste(c(tissue_name, id), collapse = "_"),
    "Annotation" = sprintf("%s. %s",
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

save(all_info, file = "all_info.rda")

# Step2: Get embedding matrix ====
# Initialize list to store embeddings
load("all_info.rda")

## Check module size distribution for each tissue =====
# 2.1: Calculate module sizes
all_info <- all_info %>%
  mutate(
    module_size = str_count(node, ";") + 1  # Count semicolons + 1 for total pathways
  )

# 2.2: Convert module_size to factor with reversed order for proper stacking
all_info <- all_info %>%
  mutate(
    # Create factor with levels in descending order so size 1 appears leftmost
    module_size_factor = factor(module_size, levels = rev(sort(unique(module_size))))
  )

# 2.3: Create summary data for plotting
plot_data <- all_info %>%
  group_by(tissue, module_size) %>%
  summarise(count = n(), .groups = "drop") %>%
  # Calculate total modules per tissue for ordering
  group_by(tissue) %>%
  mutate(total_modules = sum(count)) %>%
  ungroup() %>%
  # Convert to factor and reverse order for stacking (smallest to largest from left to right)
  mutate(module_size_factor = factor(module_size, levels = rev(sort(unique(module_size)))))
library(paletteer)
paletteer_d("MetBrewer::Signac")
colors <- c("#FBE183FF", "#F4C40FFF", "#FE9B00FF",
            "#D8443CFF", "#9B3441FF", "#DE597CFF",
            "#E87B89FF", "#E6A2A6FF", "#AA7AA1FF",
            "#9F5691FF", "#633372FF", "#1F6E9CFF",
            "#2B9B81FF", "#92C051FF")
module_size_fill_color <- colorRampPalette(colors)(16)
# 2.4: Create the stacked barplot
p <- ggplot(plot_data, aes(x = reorder(tissue, total_modules), y = count, fill = module_size_factor)) +
  geom_col(position = "stack", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = module_size_fill_color, name = "Module Size\n(# pathways)") +
  labs(
    # title = "Distribution of Functional Module Sizes by Tissue",
    x = "Tissue",
    y = "Number of Modules",
    fill = "Module Size\n(# pathways)"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  # Add total count labels at the end of bars
  geom_text(
    data = plot_data %>% group_by(tissue) %>% summarise(total = sum(count)),
    aes(x = tissue, y = total, label = total, fill = NULL),
    hjust = -0.1, size = 3, inherit.aes = FALSE
  )

p

ggsave(plot = p, filename = "56_tissues_module_number_distribution.pdf",
       width = 8, height = 8)
### =====

# Cluster functional modules ====
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
save(embedding_matrix, file = "all_525_embedding_matrix.rda")

# Step3: Calculate cosine similarity ====
calculate_cosine_sim <- function(m){
  dot_product <- m %*% t(m)
  norm_product <- sqrt(rowSums(m^2)) %*% t(sqrt(rowSums(m^2)))
  cosine_sim <- dot_product / norm_product
  return(cosine_sim)
}

sim_matrix <- calculate_cosine_sim(embedding_matrix)
save(sim_matrix, file = "all_525_fm_sim_matrix.rda")

# To prevent cluster from the same tissue cluster together -> set the sim to 0 between them
load("all_525_fm_sim_matrix.rda")
tissues <- sapply(strsplit(rownames(sim_matrix), "_"), function(x) x[1])
same_tissue_mask <- outer(tissues, tissues, "==")
sim_matrix[same_tissue_mask] <- 0

filtered_sim_matrix <- sim_matrix
save(filtered_sim_matrix, file = "filtered_sim_matrix.rda")

# Step4: Cluster functional module ====
load("01_transcriptomics/filtered_sim_matrix.rda")
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

fm_cluster_result <- merge_by_hierarchical(sim_matrix = filtered_sim_matrix,
                                           hclust.method = "ward.D2",
                                           sim.cutoff = 0.55)
## Try graph-based clustering ====
# load("filtered_sim_matrix.rda")
# node_data <- all_info |> dplyr::rename(pathways_within_module = node,
#                                        node = functional_module_name)
# node_data <- node_data |> select(node, everything())
#
# edge_data <- filtered_sim_matrix %>%
#   # Set lower triangle and diagonal to NA
#   {.[lower.tri(., diag = TRUE)] <- NA; .} %>%
#   # Convert to dataframe with row/column info
#   as.data.frame() %>%
#   tibble::rownames_to_column("from") %>%
#   tidyr::pivot_longer(-from, names_to = "to", values_to = "sim") %>%
#   # Remove NA values (lower triangle and diagonal)
#   filter(!is.na(sim))
#
# filtered_edge_data <- edge_data |> filter(sim > 0.55)
#
# graph_obj <- igraph::graph_from_data_frame(filtered_edge_data, directed = FALSE,
#                                            vertices = node_data)
# comm <- igraph::cluster_louvain(graph_obj, weights = igraph::E(graph_obj)$sim)
#
# graph_based_cluster_res <- data.frame(node = node_data$node,
#                                       module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
#
# test <- data.frame(table(graph_based_cluster_res$module))
# fm_cluster_result <- graph_based_cluster_res

test <- as.data.frame(table(fm_cluster_result$module))
sum(test$Freq == 1)
barplot(test$Freq[test$Freq > 2])

# p <- ggplot(test, aes(x = Freq)) +
#   geom_histogram(binwidth = 1,
#                  fill = "steelblue",
#                  color = "white",
#                  alpha = 0.7) +
#   stat_bin(binwidth = 5,
#            geom = "text",
#            aes(label = after_stat(count)),
#            vjust = -0.5,
#            size = 3.5) +
#   labs(title = "Distribution of Functional Module Sizes",
#        x = "Module Size (Frequency)",
#        y = "Number of Modules") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
#         axis.title = element_text(size = 12),
#         axis.text = element_text(size = 10)) +
#   scale_x_continuous(breaks = seq(0, max(test$Freq), by = 1)) +
#   ylim(0, max(table(test$Freq)) * 1.1)
#
# p

# ggsave(plot = p, filename = "fm_module_size_distribution.pdf", width = 8, height = 6)
# save(fm_cluster_result, file = "all_525_fm_cluster_result_wardD2_0.55.rda")

fm_cluster_result <- fm_cluster_result |> dplyr::rename(fm_node = node, fm_module = module)
all_info <- all_info |> dplyr::rename(fm_node = functional_module_name)
all_info <- all_info |> left_join(fm_cluster_result, by = "fm_node")
save(all_info, file = "all_525_result_with_module_info.rda")

# # Get unique tissues for each functional module
# module_tissues <- all_info %>%
#   group_by(fm_module) %>%
#   summarise(tissues = list(unique(tissue)), .groups = 'drop')
#
# # Create all pairwise combinations
# module_pairs <- expand.grid(
#   module1 = module_tissues$fm_module,
#   module2 = module_tissues$fm_module,
#   stringsAsFactors = FALSE
# )
#
# jaccard_index <- function(set1, set2) {
#   intersection <- length(intersect(set1, set2))
#   union <- length(union(set1, set2))
#   if (union == 0) return(0)  # Handle edge case
#   return(intersection / union)
# }
#
# # Calculate Jaccard index for each pair
# jaccard_results <- module_pairs %>%
#   rowwise() %>%
#   mutate(
#     tissues1 = list(module_tissues$tissues[module_tissues$fm_module == module1][[1]]),
#     tissues2 = list(module_tissues$tissues[module_tissues$fm_module == module2][[1]]),
#     jaccard_index = jaccard_index(tissues1, tissues2),
#     intersection_size = length(intersect(tissues1, tissues2)),
#     union_size = length(union(tissues1, tissues2)),
#     tissues1_size = length(tissues1),
#     tissues2_size = length(tissues2)
#   ) %>%
#   select(-tissues1, -tissues2)
# jaccard_results <- jaccard_results |> filter(module1 > module2)
# jaccard_results_0.5 <- jaccard_results |> filter(union_size > 1 & jaccard_index > 0.5)
#
# save(jaccard_results, file = "shared_tissues_fm/fm_m_overlap_tissue_jc_index.rda")
# save(jaccard_results_0.5, file = "shared_tissues_fm/fm_m_overlap_tissue_jc_index_0.5.rda")
#
# same_tissue <- jaccard_results_0.5 |> filter(jaccard_index == 1)
# same_tissue <- same_tissue |> mutate(shared_tissue = NA)
# for (i in 1:nrow(same_tissue)) {
#   module2_tissue <- all_info |> filter(fm_module == same_tissue$module2[i]) |> pull(tissue) |> unique()<- all_info |> filter(fm_module == same_tissue$module1[i]) |> pull(tissue) |> unique() |> paste(collapse = ",")
#   same_tissue$shared_tissue[i] <- shared_tissues
# }
#
# save(same_tissue, file = "shared_tissues_fm/fm_m_with_same_tissues.rda")

# Step5: Generate summary for the clustered functional modules ====
load("all_525_result_with_module_info.rda")
setwd(get_project_wd())

format_module_data_for_prompt <- function(cluster_modules) {

  # Check if input is valid
  if (nrow(cluster_modules) == 0) {
    return("No modules found in this cluster.")
  }

  formatted_modules <- c()

  for (i in 1:nrow(cluster_modules)) {
    module <- cluster_modules[i, ]

    # Format each module
    module_text <- sprintf(
      "**Module %d**: %s\n- Tissue: %s.",
      i,
      module$Annotation,
      module$tissue
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
      if (!is.null(result$`meta-module_name`) && !is.null(result$summary)) {
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
      meta_module_name = paste("Failed_Cluster", fm_module_name),
      summary = paste("Unable to generate summary for cluster", fm_module_name, "after 3 attempts.")
    )
  }

  all_annotation_result[[fm_module_name]] <- list(
    meta_module_name = result$`meta-module_name`,
    summary = result$summary
  )
}

all_res_annotation_df <- dplyr::bind_rows(all_annotation_result, .id = "fm_module")
meta_module_llm_interpreted_result <- all_res_annotation_df

save(meta_module_llm_interpreted_result, file = "2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/shared_tissues_fm/meta_module_llm_interpreted_result.rda")
export(all_res_annotation_df, file = "2_data/case_study/01_monkey/tissue_specific_analysis_results/01_transcriptomics/shared_tissues_fm/meta_module_llm_interpreted_result.xlsx")
