# functional_module_res <- get_functional_modules(
#   object = embedsim_res,
#   sim.cutoff = 0.55,
#   cluster_method = "h_ward.D2"
# )
#
#
# clustering_quality <- assess_clustering_quality(functional_module_res)
# clustering_quality$size_plot

## Step4: LLM interpretation ====
# llm_interpreted_fm_res <-
#   llm_interpret_module(
#     object = functional_module_res,
#     api_provider = "openai",
#     llm_model = "gpt-4o-mini-2024-07-18",
#     embedding_model = "text-embedding-3-small",
#     api_key = api_key,
#     embedding_output_dir = "embedding_output/",
#     module_content_number_cutoff = 1,
#     orgdb = mf.orgdb
#   )
#
# save(llm_interpreted_fm_res, file = "2_data/case_study/01_monkey/05_llm_interpreted_result/01_bulk_rna_seq_up_llm_interpreted_fm_res.rda")

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Load the data
dt <- import("2_data/case_study/01_monkey/mmc2.xlsx",
             skip = 1,
             header = TRUE,
             sheet = 2)

setwd("3_data_analysis/05_case_study/01_monkey/tissue_specific/")

# Load annotation database (uncomment if not already loaded)
ah <- AnnotationHub::AnnotationHub()
mf.orgdb <- ah[["AH119900"]] # Taxonomy ID: 9541

dt <- dt |> dplyr::rename(ensembl = GeneID)

all_variable_info <- convert_id(data = dt,
                                query_type = "gene",
                                from_id_type = "ensembl",
                                ah_id = "AH119900",
                                return_orgdb = FALSE)

all_variable_info <- all_variable_info |>
  dplyr::mutate(symbol = Symbol) |>
  dplyr::select(-Symbol)
all_variable_info$symbol[all_variable_info$symbol == "NA"] <- NA


# Get unique tissues in the dataset
unique_tissues <- unique(dt$Tissue)
print(paste("Found tissues:", paste(unique_tissues, collapse = ", ")))

# Clusters to analyze
clusters_to_analyze <- c("Cluster U", "Cluster D")

# Function to run analysis for a specific tissue and cluster
run_tissue_cluster_analysis <- function(tissue_name, cluster_name) {

  cat(paste("\n=== Processing", tissue_name, "-", cluster_name, "===\n"))

  # Filter data for specific tissue and cluster
  variable_info <- all_variable_info |>
    dplyr::filter(Tissue == tissue_name & Cluster == cluster_name)

  # Check if we have data for this combination
  if(nrow(variable_info) == 0) {
    cat(paste("No data found for", tissue_name, "-", cluster_name, ". Skipping...\n"))
    return(NULL)
  }

  cat(paste("Found", nrow(variable_info), "genes for", tissue_name, "-", cluster_name, "\n"))

  tryCatch({

    # Step 1: Enrichment analysis
    cat("Step 1: Running enrichment analysis...\n")

    enrich_pathway_res <- enrich_pathway(
      variable_info = variable_info,
      query_type = "gene",
      database = c("go", "kegg"),
      save_to_local = FALSE,
      go.orgdb = mf.orgdb,
      go.keytype = "ENSEMBL",
      go.ont = "ALL",
      kegg.organism = "mcf", # Taxid:9541
      kegg.keytype = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH"
    )

    # Step 2: Similarity calculation and clustering
    cat("Step 2: Calculating embeddings and similarity...\n")

    embedsim_res <- get_bioembedsim(
      object = enrich_pathway_res,
      api_provider = "openai",
      api_key = api_key,
      text_embedding_model = "text-embedding-3-small",
      database = c("go", "kegg")
    )

    # Create safe filename (replace spaces and special characters)
    safe_tissue_name <- gsub("[^A-Za-z0-9]", "_", tissue_name)
    safe_cluster_name <- gsub("[^A-Za-z0-9]", "_", cluster_name)
    filename <- paste0("embedsim_res_", safe_tissue_name, "_", safe_cluster_name, ".rda")

    # Save results
    cat(paste("Saving results to:", filename, "\n"))
    save(embedsim_res, file = filename)

    # Return summary information
    result_summary <- list(
      tissue = tissue_name,
      cluster = cluster_name,
      n_genes = nrow(variable_info),
      filename = filename,
      success = TRUE
    )

    cat(paste("Successfully completed analysis for", tissue_name, "-", cluster_name, "\n"))
    return(result_summary)

  }, error = function(e) {
    cat(paste("Error processing", tissue_name, "-", cluster_name, ":", e$message, "\n"))
    return(list(
      tissue = tissue_name,
      cluster = cluster_name,
      error = e$message,
      success = FALSE
    ))
  })
}

# Main analysis loop
cat("Starting tissue-specific analysis...\n")
cat("=================================\n")

# Initialize results storage
all_results <- list()
counter <- 1

# Loop through all tissues and clusters
for(tissue in unique_tissues) {
  for(cluster in clusters_to_analyze) {

    result <- run_tissue_cluster_analysis(tissue, cluster)

    if(!is.null(result)) {
      all_results[[counter]] <- result
      counter <- counter + 1
    }

    # Add a small delay to avoid overwhelming APIs if needed
    Sys.sleep(1)
  }
}

# Create summary report
cat("\n=== ANALYSIS SUMMARY ===\n")
successful_analyses <- sapply(all_results, function(x) x$success)
n_successful <- sum(successful_analyses, na.rm = TRUE)
n_total <- length(all_results)

cat(paste("Total analyses attempted:", n_total, "\n"))
cat(paste("Successful analyses:", n_successful, "\n"))
cat(paste("Failed analyses:", n_total - n_successful, "\n"))

if(n_successful > 0) {
  cat("\nSuccessful analyses:\n")
  for(i in which(successful_analyses)) {
    result <- all_results[[i]]
    cat(paste("- ", result$tissue, " - ", result$cluster,
              " (", result$n_genes, " genes) -> ", result$filename, "\n"))
  }
}

error_res <- data.frame(Tissue = NA, Cluster = NA, e_message = NA)

if(any(!successful_analyses, na.rm = TRUE)) {
  cat("\nFailed analyses:\n")
  for(i in which(!successful_analyses)) {
    result <- all_results[[i]]
    # cat(paste("- ", result$tissue, " - ", result$cluster,
    #           " Error:", result$error, "\n"))
    error_res <- rbind(error_res, data.frame(Tissue = result$tissue,
                                             Cluster = result$cluster,
                                             e_message = result$error))
  }
}

# Save summary results
summary_results <- list(
  analysis_date = Sys.time(),
  total_analyses = n_total,
  successful_analyses = n_successful,
  results_details = all_results,
  unique_tissues = unique_tissues,
  clusters_analyzed = clusters_to_analyze
)

save(summary_results, file = "tissue_cluster_analysis_summary.rda")
cat("\nSummary saved to: tissue_cluster_analysis_summary.rda\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
