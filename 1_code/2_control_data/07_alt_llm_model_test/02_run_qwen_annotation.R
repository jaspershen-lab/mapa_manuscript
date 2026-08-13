source("1_code/2_control_data/07_alt_llm_model_test/00_config.R")

required_packages <- c("httr2", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

api_key <- first_nonempty_env(c("QWEN_API_KEY", "SILICONFLOW_API_KEY"))
if (!nzchar(api_key)) {
  stop(
    "Set QWEN_API_KEY or SILICONFLOW_API_KEY before running Qwen annotation. ",
    "No API key is stored in this repository."
  )
}
configured_base_url <- first_nonempty_env("QWEN_API_BASE_URL")
if (nzchar(configured_base_url)) {
  api_base_url <- sub("/+$", "", configured_base_url)
} else {
  candidate_base_urls <- c(
    "https://api.siliconflow.cn/v1",
    "https://api.siliconflow.com/v1"
  )
  endpoint_is_authorized <- function(base_url) {
    response <- tryCatch(
      httr2::request(paste0(base_url, "/models")) |>
        httr2::req_headers(Authorization = paste("Bearer", api_key)) |>
        httr2::req_timeout(seconds = 30) |>
        httr2::req_error(is_error = function(response) FALSE) |>
        httr2::req_perform(),
      error = function(e) NULL
    )
    !is.null(response) && httr2::resp_status(response) == 200L
  }
  authorized <- vapply(candidate_base_urls, endpoint_is_authorized, logical(1L))
  if (!any(authorized)) {
    stop("The Qwen API key was rejected by both SiliconFlow endpoints.")
  }
  api_base_url <- candidate_base_urls[which(authorized)[1L]]
  message("Automatically selected Qwen API endpoint: ", api_base_url)
}
chat_endpoint <- paste0(api_base_url, "/chat/completions")
temperature <- as.numeric(first_nonempty_env("QWEN_TEMPERATURE", "0.7"))
max_tokens <- as.integer(first_nonempty_env("QWEN_MAX_TOKENS", "4096"))
if (!is.finite(temperature) || !is.finite(max_tokens)) {
  stop("QWEN_TEMPERATURE and QWEN_MAX_TOKENS must be numeric.")
}

load_one <- function(path, expected_name = NULL) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  if (!is.null(expected_name) && expected_name %in% object_names) {
    return(env[[expected_name]])
  }
  if (length(object_names) != 1L) {
    stop("Expected one object in ", path, "; found: ", paste(object_names, collapse = ", "))
  }
  env[[object_names]]
}

extract_json_object <- function(text) {
  if (is.null(text) || !length(text) || !nzchar(text)) return(NULL)
  cleaned <- trimws(as.character(text)[1L])
  cleaned <- sub("^```(?:json)?[[:space:]]*", "", cleaned, perl = TRUE)
  cleaned <- sub("[[:space:]]*```$", "", cleaned, perl = TRUE)
  first_brace <- regexpr("\\{", cleaned)[1L]
  last_brace <- max(gregexpr("\\}", cleaned)[[1L]])
  if (first_brace < 1L || last_brace < first_brace) return(NULL)
  candidate <- substr(cleaned, first_brace, last_brace)
  parsed <- tryCatch(
    jsonlite::fromJSON(candidate, simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(parsed) || is.null(parsed$module_name) || is.null(parsed$summary)) {
    return(NULL)
  }
  confidence <- suppressWarnings(as.numeric(parsed$confidence_score))[1L]
  if (!is.finite(confidence) || confidence < 0 || confidence > 1) return(NULL)
  list(
    module_name = as.character(parsed$module_name)[1L],
    summary = as.character(parsed$summary)[1L],
    confidence_score = confidence,
    json = candidate
  )
}

call_qwen <- function(messages) {
  request_body <- list(
    model = qwen_llm_model,
    messages = messages,
    max_tokens = max_tokens,
    temperature = temperature,
    n = 1
  )
  request <- httr2::request(chat_endpoint) |>
    httr2::req_method("POST") |>
    httr2::req_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json",
      Accept = "application/json"
    ) |>
    httr2::req_body_json(request_body) |>
    httr2::req_timeout(seconds = 300) |>
    httr2::req_retry(
      max_tries = 3,
      max_seconds = 600,
      is_transient = function(response) {
        httr2::resp_status(response) %in% c(408, 409, 429, 500, 502, 503, 504)
      }
    )
  response <- httr2::req_perform(request)
  response_body <- httr2::resp_body_json(response, simplifyVector = FALSE)
  content <- response_body$choices[[1L]]$message$content
  if (is.list(content) && !is.character(content)) {
    text_parts <- vapply(
      content,
      function(part) if (!is.null(part$text)) part$text else "",
      character(1L)
    )
    content <- paste(text_parts, collapse = "")
  }
  as.character(content)[1L]
}

generate_one_annotation <- function(messages) {
  raw_response <- call_qwen(messages)
  parsed <- extract_json_object(raw_response)
  if (is.null(parsed)) {
    repair_messages <- c(
      messages,
      list(
        list(role = "assistant", content = raw_response),
        list(
          role = "user",
          content = paste(
            "Return only valid JSON with exactly these fields:",
            "module_name (string), summary (string), confidence_score (number from 0 to 1)."
          )
        )
      )
    )
    raw_response <- call_qwen(repair_messages)
    parsed <- extract_json_object(raw_response)
  }
  if (is.null(parsed)) {
    stop("Qwen returned invalid annotation JSON after a format-repair attempt.")
  }
  list(
    module_name = parsed$module_name,
    summary = parsed$summary,
    confidence_score = parsed$confidence_score,
    prompt = messages,
    raw_response = raw_response,
    model = qwen_llm_model,
    api_base_url = api_base_url
  )
}

gpt_annotation_result <- load_one(gpt_annotation_file, "final_result")
if (length(gpt_annotation_result) != 16L) {
  stop("Expected 16 GPT annotations (8 real + 8 random); found ",
       length(gpt_annotation_result), ".")
}
if (any(!vapply(
  gpt_annotation_result,
  function(item) length(item$generated_name$prompt) >= 2L,
  logical(1L)
))) {
  stop("At least one GPT result does not contain the original saved prompt.")
}

cache_dir <- file.path(output_dir, "qwen_annotation_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
qwen_annotation_result <- gpt_annotation_result

for (module_name in names(gpt_annotation_result)) {
  cache_name <- paste0(gsub("[^A-Za-z0-9_-]+", "_", module_name), ".rds")
  cache_file <- file.path(cache_dir, cache_name)
  if (file.exists(cache_file)) {
    generated_name <- readRDS(cache_file)
    if (!identical(generated_name$model, qwen_llm_model)) {
      stop("Cached result was produced by a different model: ", cache_file)
    }
    message("Using cached annotation: ", module_name)
  } else {
    message("Generating Qwen annotation: ", module_name)
    generated_name <- generate_one_annotation(
      gpt_annotation_result[[module_name]]$generated_name$prompt
    )
    saveRDS(generated_name, cache_file)
  }
  qwen_annotation_result[[module_name]]$generated_name <- generated_name
  save(
    qwen_annotation_result,
    file = file.path(output_dir, "qwen_annotation_result.rda")
  )
}

qwen_annotation_table <- do.call(
  rbind,
  lapply(names(qwen_annotation_result), function(module_name) {
    generated <- qwen_annotation_result[[module_name]]$generated_name
    data.frame(
      module = module_name,
      matched_module = sub("^Random ", "", module_name),
      group = if (grepl("^Random ", module_name)) "Random" else "Real",
      module_name = generated$module_name,
      summary = generated$summary,
      confidence_score = as.numeric(generated$confidence_score),
      model = qwen_llm_model,
      stringsAsFactors = FALSE
    )
  })
)
qwen_annotation_table <- qwen_annotation_table[
  order(qwen_annotation_table$group, qwen_annotation_table$matched_module),
]
utils::write.csv(
  qwen_annotation_table,
  file.path(output_dir, "qwen_annotation_result.csv"),
  row.names = FALSE
)

run_metadata <- data.frame(
  field = c(
    "model", "api_base_url", "temperature", "max_tokens", "annotation_count",
    "prompt_source", "retrieval_source", "completed_at"
  ),
  value = c(
    qwen_llm_model, api_base_url, temperature, max_tokens,
    nrow(qwen_annotation_table),
    "saved prompts in GPT final_result.Rdata",
    "saved related_paper and pubmed_result in GPT final_result.Rdata",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  run_metadata,
  file.path(output_dir, "qwen_annotation_run_metadata.csv"),
  row.names = FALSE
)
message("Qwen annotation complete: ", nrow(qwen_annotation_table), " modules.")
