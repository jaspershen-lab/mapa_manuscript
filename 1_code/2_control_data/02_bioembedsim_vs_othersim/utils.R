
# Extract GO term info (pathway_id, name, definition(description)) ====
get_go_info <- function(go_ids) {
  #### Generate GO terms and annotated Entrez Gene identifiers dict
  go2egs <- NULL
  eg2genename <- NULL

  # # Only retrieve gene-related information if include_gene_name is TRUE
  # if (include_gene_name) {
  #   go2egs <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go_ids, columns = "ENTREZID", keytype = "GOALL"))
  #   #### Generate entrez ID to gene name dict
  #   eg2genename <-
  #     suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go2egs$ENTREZID, columns = "GENENAME")) %>%
  #     dplyr::distinct()
  # }

  chunk_size <- 100
  chunks <- split(go_ids, ceiling(seq_along(go_ids) / chunk_size))
  go_info <- list()
  for (i in 1:length(chunks)) {
    sub_go_info <-
      quickgo_api(go_ids = chunks[[i]]) %>%
      purrr::map(
        function(x) {
          # Initialize empty gene name string
          go2genename <- ""

          # # Only process gene names if include_gene_name is TRUE
          # if (include_gene_name) {
          #   # Get GENENAME: GO ID -> ENTREZID -> GENENAME
          #   go2genename <- go2egs %>%
          #     dplyr::filter(GOALL == x$id) %>%
          #     dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
          #     dplyr::left_join(eg2genename, by = "ENTREZID") %>%
          #     dplyr::pull(GENENAME) %>%
          #     paste0(collapse = ", ")
          # }

          # # Get PMID
          # if ("xrefs" %in% names(x$definition)){
          #   # Extract 'dbId' from each element in the list
          #   pmid_list <- lapply(x$definition$xrefs,
          #                       function(r) {
          #                         if(r$dbCode == "PMID") {
          #                           pmid <- r$dbId
          #                         } else {
          #                           pmid <- NULL
          #                         }})
          #   # Convert the list of dbId values into a string
          #   pmid <- paste(unlist(pmid_list), collapse = ",")
          # } else {
          #   pmid <- ""
          # }

          # Collect all info
          all_info <- list(
            "id" = x$id,
            "term_name" = x$name,
            "term_definition" = x$definition$text
          )

          # # Only add annotated_genename field if include_gene_name is TRUE
          # if (include_gene_name) {
          #   all_info$annotated_genename <- go2genename
          # }

          return(all_info)
        }
      )
    go_info <- c(go_info, sub_go_info)
  }
  return(go_info)
}

# Get core information about a list of terms based on their ids
quickgo_api <- function(go_ids) {
  url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", paste(go_ids, collapse = ","))
  tryCatch(
    expr = {
      # Create a request object
      req <- httr2::request(url)

      resp <- req %>%
        httr2::req_headers("Accept" = "application/json") %>%
        httr2::req_retry(max_tries = 3,
                         max_seconds = 60,
                         # Condition for when to retry
                         after = \(resp) is.null(resp) && resp$status_code != 200
        ) %>%
        httr2::req_perform()

      # Parse json to get a list
      info <- httr2::resp_body_json(resp)$result
    },
    error = function(e) {
      warning("Failed to get info from QuickGO after 3 retries:", e$message)
      NULL
    })
}

# Extract KEGG info (pathway_id, name, definition(description)) ====
get_kegg_pathway_info <- function(kegg_ids){
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

              # # Only get annotated Entrez ID if include_gene_name is TRUE
              # if (include_gene_name && "GENE" %in% names(x)) {
              #   kegg2genename <-
              #     x$GENE[seq(2, length(x$GENE), 2)] %>%
              #     stringr::str_match(pattern = ";\\s*(.*?)\\s*\\[") %>%
              #     as.data.frame() %>%
              #     dplyr::pull(V2) %>%
              #     paste0(collapse = ",")
              # }

              # # Get PMID
              # if ("REFERENCE" %in% names(x)) {
              #   pmid_list <- lapply(x$REFERENCE,
              #                       function(r) {
              #                         if (grepl("^PMID", r$REFERENCE)) {
              #                           pmid <- gsub("[^:]+:", "", r$REFERENCE)
              #                         } else {
              #                           pmid <- ""
              #                         }
              #                       })
              #   pmid <- paste(unlist(pmid_list), collapse = ",")
              # } else {
              #   pmid <- ""
              # }

              # Collect base info (always included)
              all_info <- list(
                "id" = unname(x$ENTRY),
                # "term_name" = sub(" - Homo sapiens \\(human\\)$", "", x$NAME),
                "term_name" = sub(" - [^(]+\\([^)]+\\)$", "", x$NAME),
                "term_definition" = paste(x$DESCRIPTION, collapse = " ")
              )

              # # Only add gene name information if requested
              # if (include_gene_name) {
              #   all_info$annotated_genename <- kegg2genename
              # }

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

            # # Only get annotated Entrez ID if include_gene_name is TRUE
            # if (include_gene_name && "GENE" %in% names(entry)) {
            #   kegg2genename <-
            #     entry$GENE[seq(2, length(entry$GENE), 2)] %>%
            #     stringr::str_match(pattern = ";\\s*(.*?)\\s*\\[") %>%
            #     as.data.frame() %>%
            #     dplyr::pull(V2) %>%
            #     paste0(collapse = ",")
            # }

            # Collect base info (always included)
            all_info <- list(
              "id" = unname(entry$ENTRY),
              # "term_name" = sub(" - Homo sapiens \\(human\\)$", "", entry$NAME),
              "term_name" = sub(" - [^(]+\\([^)]+\\)$", "", entry$NAME),
              "term_definition" = paste(entry$DESCRIPTION, collapse = " ")
            )

            # Only add gene name information if requested
            # if (include_gene_name) {
            #   all_info$annotated_genename <- kegg2genename
            # }

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

# Extract Reactome info (pathway_id, name, definition) ====
get_reactome_pathway_info <- function(reactome_ids) {
  #### Initialize gene-related variables
  reactome2egs <- NULL
  eg2genename <- NULL

  #### Only retrieve gene-related information if include_gene_name is TRUE
  # if (include_gene_name) {
  #   #### Generate Reactome terms and annotated Entrez Gene identifiers dict
  #   suppressMessages(reactome2egs <- AnnotationDbi::select(reactome.db, keys = reactome_ids, columns = "ENTREZID", keytype = "PATHID"))
  #   #### Generate entrez ID to gene name dict
  #   eg2genename <-
  #     suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = reactome2egs$ENTREZID, columns = "GENENAME")) %>%
  #     dplyr::distinct()
  # }

  chunk_size <- 20
  chunks <- split(reactome_ids, ceiling(seq_along(reactome_ids) / chunk_size))
  reactome_info <- list()

  for (i in 1:length(chunks)) {
    reactome_api_result <- suppressMessages(rbioapi::rba_reactome_query(ids = chunks[[i]]))
    if (length(chunks[[i]]) == 1) {reactome_api_result <- list(reactome_api_result)}

    sub_reactome_info <-
      reactome_api_result %>%
      purrr::map(
        function(x) {
          # Initialize empty gene name string
          reactome2genename <- ""

          # # Only process gene names if include_gene_name is TRUE
          # if (include_gene_name) {
          #   # Get GENENAME: Reactome ID -> ENTREZID -> GENENAME
          #   reactome2genename <- reactome2egs %>%
          #     dplyr::filter(PATHID == x$stId) %>%
          #     dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
          #     dplyr::left_join(eg2genename, by = "ENTREZID") %>%
          #     dplyr::pull(GENENAME) %>%
          #     paste0(collapse = ", ")
          # }

          # # Get PMID
          # if ("literatureReference" %in% names(x)) {
          #   pmid_list <- lapply(x$literatureReference,
          #                       function(r) {
          #                         if ("pubMedIdentifier" %in% names(r)) {
          #                           pmid <- r$pubMedIdentifier
          #                         } else {
          #                           pmid <- ""
          #                         }
          #                       })
          #   pmid <- paste(unlist(pmid_list), collapse = ",")
          # } else {
          #   pmid <- ""
          # }

          # Collect base info (always included)
          all_info <- list(
            "id" = x$stId,
            "term_name" = x$displayName,
            "term_definition" = gsub("(<BR>|<br>)", "", x$summation[[1]]$text)
          )

          # # Only add gene name information if requested
          # if (include_gene_name) {
          #   all_info$annotated_genename <- reactome2genename
          # }

          return(all_info)
        }
      )
    reactome_info <- c(reactome_info, sub_reactome_info)
  }
  return(reactome_info)
}

# Combine information (term name, term description, gene names) into a string ====
combine_info <- function(info) {
  info %>%
    purrr::map(
      function(x) {
        # if (include_gene_name) {
        #   text_info <- sprintf("Pathway name: %s\nDefinition: %s\nAnnotated gene names: %s", x$term_name, x$term_definition, x$annotated_genename)
        # } else {
        #   text_info <- sprintf("%s: %s", x$term_name, x$term_definition)
        # }
        text_info <- sprintf("%s: %s", x$term_name, x$term_definition)

        text <- list(
          "id" = x$id,
          "text_info" = text_info
        )

        return(text)
      }
    )
}

# Get embeddings ====
get_embedding_matrix <-
  function(
    text,
    api_provider = c("openai", "gemini"),
    text_embedding_model,
    api_key) {

    if (missing(api_provider)) {
      stop("api_provider is required.")
    }
    api_provider <- match.arg(api_provider)

    if (missing(text_embedding_model)) {
      stop("text_embedding_model is required.")
    }

    if (missing(api_key)) {
      stop("api_key is required.")
    }

    if (api_provider == "openai") {
      embedding_matrix <-
        text %>%
        purrr::map_df(function(x) {
          token_num <-
            rtiktoken::get_token_count(text = x$text_info, model = "text-embedding-3-small")
          # max input for openai embedding model is 8191 (March 05, 2025)
          # if (include_gene_name == TRUE & token_num > 8000) {
          #   x$text_info <- sub(pattern = "\nAnnotated gene names.*$", "", x$text_info)
          #   message(paste0("Text information of ", x$id, " exceeded token limit and annotated gene names were removed."))
          # }
          embedding <-
            get_openai_embedding_internal(input_text = x$text_info, text_embedding_model = text_embedding_model, api_key = api_key) %>%
            t() %>%
            as.data.frame()
          rownames(embedding) <- x$id

          embedding
        }) %>%
        as.matrix()
    } else if (api_provider == "gemini") {
      message("Use the tokenizer(i.e. Cl100kBase) for OpenAI's embedding model to approximate the token number for Gemini model.")
      embedding_matrix <-
        text %>%
        purrr::map_df(function(x) {
          token_num <-
            rtiktoken::get_token_count(text = x$text_info, model = "text-embedding-3-small")
          if (token_num > 1500) {
            query <- paste0("According to the following description: ", x$text_info, " Please generate a short summary about 100 words to extract the core function information of the pathway described.")
            llm_summary <- gemini_api_call(
              input_text = query,
              api_key = api_key
            )
            # Considering the token limitation (2048 tokens) of text embedding model, use llm to summarize the pathway description
            message(paste0("For pathway ID: ", x$id, ", use the summary generated by gemini-1.5-flash to do text embedding since pathway information exceeded token limitation."))
            x$text_info <- paste0("Pathway name:", x$id, "\n", llm_summary)
          }

          embedding <-
            get_gemini_embedding_internal(input_text = x$text_info, text_embedding_model = text_embedding_model, api_key = api_key) %>%
            t() %>%
            as.data.frame()
          rownames(embedding) <- x$id

          embedding
        }) %>%
        as.matrix()
    }

    return(embedding_matrix)
  }


get_openai_embedding_internal <-
  function(input_text, text_embedding_model = "text-embedding-3-small", api_key) {

    url <- "https://api.openai.com/v1/embeddings"

    # Body specifying model and text
    data <- list(
      model = text_embedding_model,  # Use the specified model or default
      input = input_text  # Input text
    )

    embedding <- tryCatch(
      expr = {
        # Create a request object
        req <- httr2::request(url)

        resp <- req %>%
          httr2::req_auth_bearer_token(token = api_key) %>%
          httr2::req_body_json(data = data) %>%
          httr2::req_retry(max_tries = 3,
                           max_seconds = 60,
                           after = \(resp) is.null(resp) && resp$status_code != 200 # Condition for when to retry
          ) %>%
          httr2::req_perform()

        embedding <-
          # Parsed JSON -> an embedding, a list
          httr2::resp_body_json(resp) %>%
          {.$data[[1]]$embedding} %>%
          unlist()
      },
      error = function(e) {
        warning("Failed to get embedding after 3 retries:", e$message)
        NULL
      })

    return(embedding)
  }

get_gemini_embedding_internal <-
  function(input_text, text_embedding_model = "text-embedding-004", api_key) {

    url <- paste0("https://generativelanguage.googleapis.com/v1beta/models/", text_embedding_model, ":embedContent?key=", api_key)

    embedding <- tryCatch(
      expr = {
        # Create a request object
        req <- httr2::request(url)

        # Get response
        resp <- req %>%
          httr2::req_headers("Content-Type" = "application/json") %>%
          httr2::req_body_json(
            list(
              model = paste0("models/", text_embedding_model),
              content = list(parts = list(list(text = input_text)))
            )
          ) %>%
          httr2::req_retry(max_tries = 3,
                           max_seconds = 60,
                           after = \(resp) is.null(resp) && resp$status_code != 200 # Condition for when to retry
          ) %>%
          httr2::req_perform()

        # Parse the JSON response
        embedding <- resp %>%
          httr2::resp_body_json() %>%
          {.$embedding$values} %>%
          unlist()
      },
      error = function(e) {
        warning("Failed to get embedding after 3 retries:", e$message)
        NULL
      }
    )

    return(embedding)
  }

gemini_api_call <-
  function(input_text,
           model = "gemini-1.5-flash",
           api_key,
           temperature = 1,
           topK = 40,
           topP = 0.95,
           maxOutputTokens = 8192,
           responseMimeType = "text/plain") {
    # Construct the full URL with the API key and model as query parameters
    url <- paste0("https://generativelanguage.googleapis.com/v1beta/models/", model, ":generateContent?key=", api_key)

    # Create the request body
    request_body <- list(
      contents = list(
        list(
          role = "user",
          parts = list(list(text = input_text))
        )
      ),
      generationConfig = list(
        temperature = temperature,
        topK = topK,
        topP = topP,
        maxOutputTokens = maxOutputTokens,
        responseMimeType = responseMimeType
      )
    )

    # Create and send the request using httr2
    response <- httr2::request(url) %>%
      httr2::req_headers("Content-Type" = "application/json") %>%
      httr2::req_body_json(request_body) %>%
      httr2::req_perform()

    # Check for errors
    if (httr2::resp_is_error(response)) {
      stop(paste("API request failed with status:", resp_status(response)))
    }

    # Parse and return the JSON response
    result <- response %>% httr2::resp_body_json()

    # Extract generated text
    return(result$candidates[[1]]$content$parts[[1]]$text)
  }

# Calculate pairwise cosine similarity ====
calculate_cosine_sim <- function(m){
  dot_product <- m %*% t(m)
  norm_product <- sqrt(rowSums(m^2)) %*% t(sqrt(rowSums(m^2)))
  cosine_sim <- dot_product / norm_product
  return(cosine_sim)
}

# Calculate similarity based on gene overlap ====
term_similarity_internal <-
  function(gl,
           measure.method = c("jaccard", "dice", "overlap", "kappa"),
           remove_negative = TRUE) {

    measure.method <- match.arg(measure.method)

    all <- unique(unlist(gl))
    gl <- lapply(gl, function(x) as.numeric(factor(x, levels = all)))
    n <- length(gl)

    pathway_gene_m <- matrix(0, ncol = length(all), nrow = n)
    for(i in seq_len(n)) {
      pathway_gene_m[i, gl[[i]]] = 1
    }
    pathway_gene_m <- as(pathway_gene_m, "sparseMatrix")

    if(measure.method == "kappa") {
      mat <- kappa_dist(pathway_gene_m, remove_negative = remove_negative)
    } else if(measure.method == "overlap") {
      mat <- overlap_dist(pathway_gene_m)
    } else {
      mat <- proxyC::simil(pathway_gene_m, method = measure.method)
    }

    sim_matrix <-  as.matrix(mat)
    diag(sim_matrix) <- 1
    rownames(sim_matrix) = colnames(sim_matrix) = names(gl)

    return(sim_matrix)
  }

kappa_dist <- function(m, remove_negative = TRUE) {
  tab <- ncol(m)
  po <- proxyC::simil(m, method = "simple matching")
  m_yes <- Matrix::rowSums(m)
  m_no <- abs(Matrix::rowSums(m - 1))
  pe <- (outer(m_yes, m_yes, FUN = "*") + outer(m_no, m_no, FUN = "*"))/tab^2
  k <- (po - pe)/(1 - pe)
  if(remove_negative) k[k < 0] <- 0
  return(k)
}

overlap_dist <- function(m) {
  n = Matrix::rowSums(m)
  proxyC::simil(m, method = "dice")*outer(n, n, FUN = "+")/2/outer(n, n, pmin)
}
