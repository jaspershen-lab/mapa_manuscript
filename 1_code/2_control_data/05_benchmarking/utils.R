get_kegg_pathway_id_info <- function(kegg_ids){
  chunk_size <- 10
  chunks <- split(kegg_ids, ceiling(seq_along(kegg_ids) / chunk_size))
  kegg_info <- list()

  for (i in 1:length(chunks)) {
    sub_kegg_info <-
      tryCatch({
        # Try to process the chunk as a batch first
        KEGGREST::keggGet(dbentries = chunks[[i]]) %>%
          purrr::map(
            function(x) {
              # Initialize kegg2genename as empty
              kegg2genename <- NA
              kegg2class <- NA

              # Only get annotated Entrez ID if include_gene_name is TRUE
              if ("GENE" %in% names(x)) {
                kegg2genename <-
                  x$GENE[seq(1, length(x$GENE), 2)] |>
                  paste(collapse = "/")
              }

              if ("CLASS" %in% names(x)) {
                kegg2class <- x$CLASS
              }

              # Collect base info (always included)
              all_info <- list(
                "id" = unname(x$ENTRY),
                "class" = kegg2class,
                "term_name" = sub(" - [^(]+\\([^)]+\\)$", "", x$NAME),
                "term_definition" = paste(x$DESCRIPTION, collapse = " "),
                "annotated_genes" = kegg2genename
              )

              return(all_info)
            })
      }, error = function(e) {
        # If batch processing fails, process each ID individually
        message("Batch processing failed, trying individual IDs: ", e$message)

        # Process each ID in the current chunk individually
        individual_results <- list()
        for (kegg_id in chunks[[i]]) {
          single_result <- tryCatch({
            # Get a single KEGG entry
            entry <- KEGGREST::keggGet(dbentries = kegg_id)[[1]]

            # Initialize kegg2genename as empty
            kegg2genename <- NA
            kegg2class <- NA

            # # Only get annotated Entrez ID if include_gene_name is TRUE
            if ("GENE" %in% names(entry)) {
              kegg2genename <-
                x$GENE[seq(1, length(x$GENE), 2)] |>
                paste(collapse = "/")
            }

            if ("CLASS" %in% names(x)) {
              kegg2class <- x$CLASS
            }

            # Collect base info (always included)
            all_info <- list(
              "id" = unname(x$ENTRY),
              "class" = kegg2class,
              "term_name" = sub(" - [^(]+\\([^)]+\\)$", "", x$NAME),
              "term_definition" = paste(x$DESCRIPTION, collapse = " "),
              "annotated_genes" = kegg2genename
            )

            all_info
          }, error = function(e2) {
            message("  Failed to retrieve KEGG ID: ", kegg_id, " (", e2$message, ")")
            return(NULL)
          })

          # Add non-NULL results to the list
          if (!is.null(single_result)) {
            individual_results <- c(individual_results, list(single_result))
          }

          # Brief pause to avoid overwhelming the API
          Sys.sleep(0.1)
        }

        return(individual_results)
      })

    kegg_info <- c(kegg_info, sub_kegg_info)
  }
  return(kegg_info)
}

get_go_term_info <- function(go_ids) {

  chunk_size <- 100
  chunks <- split(go_ids, ceiling(seq_along(go_ids) / chunk_size))
  go_info <- list()
  for (i in 1:length(chunks)) {
    sub_go_info <-
      mapa::quickgo_api(go_ids = chunks[[i]]) %>%
      purrr::map(
        function(x) {

          # Collect all info
          all_info <- list(
            "id" = x$id,
            "sub_ontology" = x$aspect,
            "term_name" = x$name,
            "term_definition" = x$definition$text
          )

          return(all_info)
        }
      )
    go_info <- c(go_info, sub_go_info)
  }
  return(go_info)
}
