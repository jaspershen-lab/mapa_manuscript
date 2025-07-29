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
save(all_variable_info, file = "2_data/case_study/01_monkey/01_bulk_rna_seq/all_variable_info.rda")

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


# Check failed results ====
load("2_data/case_study/01_monkey/01_bulk_rna_seq/all_variable_info.rda")
setwd("3_data_analysis/05_case_study/01_monkey/tissue_specific/")
load("error_res.rda")

error_res_no_pathways <- error_res |> dplyr::filter(grepl("^No pathways", e_message))
error_res_check <- error_res |> dplyr::filter(!grepl("^No pathways", e_message))

error_check_all_results <- list()
counter <- 1

for(i in 1:nrow(error_res_check)) {
  tissue_name <- error_res_check$Tissue[i]
  cluster <- error_res_check$Cluster[i]

  result <- run_tissue_cluster_analysis(tissue_name = tissue_name,
                                        cluster_name = cluster)

  if(!is.null(result)) {
    error_check_all_results[[counter]] <- result
    counter <- counter + 1
  }

  Sys.sleep(1)
}

# Create summary report
cat("\n=== ANALYSIS SUMMARY ===\n")
successful_analyses <- sapply(error_check_all_results, function(x) x$success)
n_successful <- sum(successful_analyses, na.rm = TRUE)
n_total <- length(error_check_all_results)

cat(paste("Total analyses attempted:", n_total, "\n"))
cat(paste("Successful analyses:", n_successful, "\n"))
cat(paste("Failed analyses:", n_total - n_successful, "\n"))

error_no_pathway <- error_res |>
  dplyr::filter(grepl("^No pathways", e_message))

save(error_no_pathway, file = "error_no_pathway.rda")

load("tissue_cluster_analysis_summary.rda")
successful_results <- lapply(summary_results$results_details,
                             function(x) {
                               tissue <- x$tissue
                               cluster <- x$cluster

                               data.frame("tissue" = tissue, "cluster" = cluster)
                             }) %>% do.call(rbind, .)

successful_results <- anti_join(
  successful_results,
  error_no_pathway,
  by = c(tissue  = "Tissue",
         cluster = "Cluster")
)

save(successful_results, file = "successful_results.rda")

# Clustering and llm interpretation ====
# Load your successful_results dataframe

create_safe_filename <- function(tissue_name, cluster_name) {
  safe_tissue_name <- gsub("[^A-Za-z0-9]", "_", tissue_name)
  safe_cluster_name <- gsub("[^A-Za-z0-9]", "_", cluster_name)
  filename <- paste0("embedsim_res_", safe_tissue_name, "_", safe_cluster_name, ".rda")
  return(filename)
}

# Function to process a single embedsim_res file
process_single_file <- function(tissue_name, cluster_name, api_key, mf.orgdb) {

  # Create filename
  filename <- create_safe_filename(tissue_name, cluster_name)

  # Check if file exists
  if (!file.exists(filename)) {
    cat("Warning: File", filename, "not found. Skipping...\n")
    return(NULL)
  }

  cat("Processing:", tissue_name, "-", cluster_name, "\n")

  # Load the embedsim_res object
  load(filename)

  # Check if embedsim_res object exists
  if (!exists("embedsim_res")) {
    cat("Warning: embedsim_res object not found in", filename, ". Skipping...\n")
    return(NULL)
  }

  tryCatch({
    # Step 1: Get functional modules
    cat("  Getting functional modules...\n")
    functional_module_res <- get_functional_modules(
      object = embedsim_res,
      sim.cutoff = 0.55,
      cluster_method = "h_ward.D2"
    )

    # Step 2: LLM interpretation
    cat("  Performing LLM interpretation...\n")
    llm_interpreted_fm_res <- llm_interpret_module(
      object = functional_module_res,
      api_provider = "openai",
      llm_model = "gpt-4o-mini-2024-07-18",
      embedding_model = "text-embedding-3-small",
      api_key = api_key,
      embedding_output_dir = "embedding_output/",
      module_content_number_cutoff = 1,
      orgdb = mf.orgdb
    )


    safe_tissue_name <- gsub("[^A-Za-z0-9]", "_", tissue_name)
    safe_cluster_name <- gsub("[^A-Za-z0-9]", "_", cluster_name)

    # Save functional modules result
    fm_filename <- paste0("llm_interpreted_fm_res_", safe_tissue_name, "_", safe_cluster_name, ".rda")
    save(llm_interpreted_fm_res, file = fm_filename)

    cat("  Completed successfully. Results saved as:", fm_filename, "\n")

    # Return summary information
    return(list(
      tissue = tissue_name,
      cluster = cluster_name,
      status = "success",
      fm_file = fm_filename
    ))

  }, error = function(e) {
    cat("  Error processing", tissue_name, "-", cluster_name, ":", e$message, "\n\n")
    return(list(
      tissue = tissue_name,
      cluster = cluster_name,
      status = "error",
      error_message = e$message
    ))
  })
}

# Main processing function
process_all_files <- function(successful_results, api_key, mf.orgdb) {

  # Create embedding output directory if it doesn't exist
  if (!dir.exists("embedding_output/")) {
    dir.create("embedding_output/", recursive = TRUE)
  }

  # Initialize results tracking
  processing_results <- list()

  # Process each row in successful_results
  for (i in 1:nrow(successful_results)) {
    tissue_name <- successful_results$tissue[i]
    cluster_name <- successful_results$cluster[i]

    result <- process_single_file(tissue_name, cluster_name, api_key, mf.orgdb)
    processing_results[[i]] <- result
  }

  # Convert results to dataframe
  results_df <- do.call(rbind, lapply(processing_results, function(x) {
    if (is.null(x)) return(NULL)
    data.frame(
      tissue = x$tissue,
      cluster = x$cluster,
      status = x$status,
      fm_file = ifelse(is.null(x$fm_file), NA, x$fm_file),
      error_message = ifelse(is.null(x$error_message), NA, x$error_message),
      stringsAsFactors = FALSE
    )
  }))

  return(results_df)
}


# Run the processing
cat("Starting batch processing of", nrow(successful_results), "files...\n\n")
final_results <- process_all_files(successful_results, api_key, mf.orgdb)

# Print summary
cat("Processing complete!\n")
cat("Success:", sum(final_results$status == "success", na.rm = TRUE), "\n")
cat("Errors:", sum(final_results$status == "error", na.rm = TRUE), "\n")

# Save the processing results
save(final_results, file = "llm_interpreted_res_processing_results_summary.rda")

# Get functional module result for
# module_content_number_cutoff should be smaller than the maximum of
# all module content numbers in your functional module result.

error_module_num_cutoff <- final_results |> dplyr::filter(grepl("^module_content_number_cutoff", error_message))

process_single_file_2 <- function(tissue_name, cluster_name, mf.orgdb) {

  # Create filename
  filename <- create_safe_filename(tissue_name, cluster_name)

  # Check if file exists
  if (!file.exists(filename)) {
    cat("Warning: File", filename, "not found. Skipping...\n")
    return(NULL)
  }

  cat("Processing:", tissue_name, "-", cluster_name, "\n")

  # Load the embedsim_res object
  load(filename)

  # Check if embedsim_res object exists
  if (!exists("embedsim_res")) {
    cat("Warning: embedsim_res object not found in", filename, ". Skipping...\n")
    return(NULL)
  }

  tryCatch({
    # Step 1: Get functional modules
    cat("  Getting functional modules...\n")
    functional_module_res <- get_functional_modules(
      object = embedsim_res,
      sim.cutoff = 0.55,
      cluster_method = "h_ward.D2"
    )

    safe_tissue_name <- gsub("[^A-Za-z0-9]", "_", tissue_name)
    safe_cluster_name <- gsub("[^A-Za-z0-9]", "_", cluster_name)

    # Save functional modules result
    fm_filename <- paste0("fm_res_", safe_tissue_name, "_", safe_cluster_name, ".rda")
    save(functional_module_res, file = fm_filename)

    cat("  Completed successfully. Results saved as:", fm_filename, "\n")

    # Return summary information
    return(list(
      tissue = tissue_name,
      cluster = cluster_name,
      status = "success",
      fm_file = fm_filename
    ))

  }, error = function(e) {
    cat("  Error processing", tissue_name, "-", cluster_name, ":", e$message, "\n\n")
    return(list(
      tissue = tissue_name,
      cluster = cluster_name,
      status = "error",
      error_message = e$message
    ))
  })
}

process_all_files_2 <- function(error_results, mf.orgdb) {

  # Initialize results tracking
  processing_results <- list()

  # Process each row in successful_results
  for (i in 1:nrow(error_results)) {
    tissue_name <- error_results$tissue[i]
    cluster_name <- error_results$cluster[i]

    result <- process_single_file_2(tissue_name, cluster_name, mf.orgdb)
    processing_results[[i]] <- result
  }

  # Convert results to dataframe
  results_df <- do.call(rbind, lapply(processing_results, function(x) {
    if (is.null(x)) return(NULL)
    data.frame(
      tissue = x$tissue,
      cluster = x$cluster,
      status = x$status,
      fm_file = ifelse(is.null(x$fm_file), NA, x$fm_file),
      error_message = ifelse(is.null(x$error_message), NA, x$error_message),
      stringsAsFactors = FALSE
    )
  }))

  return(results_df)
}

final_results2 <- process_all_files_2(error_results = error_module_num_cutoff,
                                      mf.orgdb = mf.orgdb)

# Print summary
cat("Success:", sum(final_results2$status == "success", na.rm = TRUE), "\n")
cat("Errors:", sum(final_results2$status == "error", na.rm = TRUE), "\n")
