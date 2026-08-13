source("1_code/2_control_data/07_alt_llm_model_test/00_config.R")

required_packages <- c(
  "dplyr", "ggplot2", "ggrepel", "httr2", "jsonlite", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
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

flatten_annotations <- function(result, llm) {
  do.call(rbind, lapply(names(result), function(module) {
    generated <- result[[module]]$generated_name
    data.frame(
      llm = llm,
      module = module,
      matched_module = sub("^Random ", "", module),
      group = if (grepl("^Random ", module)) "Random" else "Real",
      module_name = as.character(generated$module_name)[1L],
      summary = as.character(generated$summary)[1L],
      confidence_score = as.numeric(generated$confidence_score)[1L],
      text = paste0(generated$module_name, ". ", generated$summary),
      stringsAsFactors = FALSE
    )
  }))
}

cosine_similarity <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  denominator <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  max(-1, min(1, sum(x * y) / denominator))
}

request_openai_compatible_embeddings <- function(text, model) {
  api_key <- first_nonempty_env(c("OPENAI_API_KEY", "ANNOTATION_EMBEDDING_API_KEY"))
  if (!nzchar(api_key)) {
    stop(
      "No cached annotation embeddings are available. Set OPENAI_API_KEY ",
      "to run the existing text-embedding-3-small semantic-comparison workflow."
    )
  }
  api_base_url <- sub(
    "/+$", "",
    first_nonempty_env(
      "ANNOTATION_EMBEDDING_API_BASE_URL", default = "https://api.openai.com/v1"
    )
  )
  response <- httr2::request(paste0(api_base_url, "/embeddings")) |>
    httr2::req_method("POST") |>
    httr2::req_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json",
      Accept = "application/json"
    ) |>
    httr2::req_body_json(list(model = model, input = unname(text))) |>
    httr2::req_timeout(seconds = 300) |>
    httr2::req_retry(
      max_tries = 3,
      max_seconds = 600,
      is_transient = function(response) {
        httr2::resp_status(response) %in% c(408, 409, 429, 500, 502, 503, 504)
      }
    ) |>
    httr2::req_perform()
  body <- httr2::resp_body_json(response, simplifyVector = FALSE)
  ordered <- body$data[order(vapply(body$data, `[[`, integer(1L), "index"))]
  vectors <- lapply(ordered, `[[`, "embedding")
  if (length(vectors) != length(text) || length(unique(lengths(vectors))) != 1L) {
    stop("Embedding API returned an unexpected number or dimension of vectors.")
  }
  matrix <- do.call(rbind, vectors)
  rownames(matrix) <- names(text)
  matrix
}

qwen_annotation_file <- file.path(output_dir, "qwen_annotation_result.rda")
qwen_result <- load_one(qwen_annotation_file, "qwen_annotation_result")
gpt_result <- load_one(gpt_annotation_file, "final_result")

annotation_data <- dplyr::bind_rows(
  flatten_annotations(gpt_result, "GPT-4o-mini"),
  flatten_annotations(qwen_result, "Qwen3.6-35B-A3B")
) |>
  dplyr::mutate(
    embedding_id = paste(llm, module, sep = "::"),
    group = factor(group, levels = c("Random", "Real"))
  )
if (anyNA(annotation_data$confidence_score)) {
  stop("At least one LLM confidence score is missing or non-numeric.")
}

embedding_cache_file <- file.path(output_dir, "annotation_embeddings.rds")
embedding_text <- stats::setNames(annotation_data$text, annotation_data$embedding_id)
use_cache <- FALSE
if (file.exists(embedding_cache_file)) {
  cache <- readRDS(embedding_cache_file)
  use_cache <- identical(cache$model, annotation_embedding_model) &&
    identical(cache$text, embedding_text)
}
if (use_cache) {
  annotation_embedding_matrix <- cache$matrix
} else {
  annotation_embedding_matrix <- request_openai_compatible_embeddings(
    embedding_text, annotation_embedding_model
  )
  saveRDS(
    list(
      model = annotation_embedding_model,
      text = embedding_text,
      matrix = annotation_embedding_matrix,
      generated_at = Sys.time()
    ),
    embedding_cache_file
  )
}

qwen_vs_gpt_semantic_similarity <- lapply(
  unique(annotation_data$module),
  function(module) {
    rows <- annotation_data[annotation_data$module == module, ]
    qwen_id <- rows$embedding_id[rows$llm == "Qwen3.6-35B-A3B"]
    gpt_id <- rows$embedding_id[rows$llm == "GPT-4o-mini"]
    data.frame(
      module = module,
      matched_module = sub("^Random ", "", module),
      group = if (grepl("^Random ", module)) "Random" else "Real",
      qwen_module_name = rows$module_name[rows$llm == "Qwen3.6-35B-A3B"],
      gpt_module_name = rows$module_name[rows$llm == "GPT-4o-mini"],
      cosine_similarity = cosine_similarity(
        annotation_embedding_matrix[qwen_id, ],
        annotation_embedding_matrix[gpt_id, ]
      ),
      stringsAsFactors = FALSE
    )
  }
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(group = factor(group, levels = c("Random", "Real")))

qwen_vs_gpt_semantic_summary <- qwen_vs_gpt_semantic_similarity |>
  dplyr::group_by(group) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(cosine_similarity),
    sd = stats::sd(cosine_similarity),
    median = stats::median(cosine_similarity),
    q1 = as.numeric(stats::quantile(cosine_similarity, 0.25)),
    q3 = as.numeric(stats::quantile(cosine_similarity, 0.75)),
    .groups = "drop"
  ) |>
  dplyr::bind_rows(
    qwen_vs_gpt_semantic_similarity |>
      dplyr::summarise(
        group = "Overall",
        n = dplyr::n(),
        mean = mean(cosine_similarity),
        sd = stats::sd(cosine_similarity),
        median = stats::median(cosine_similarity),
        q1 = as.numeric(stats::quantile(cosine_similarity, 0.25)),
        q3 = as.numeric(stats::quantile(cosine_similarity, 0.75))
      )
  )

all_result_long_data <- load_one(expert_embedding_file, "all_result_long_data")
tool_methods <- c(
  "mapa_cluster_label", "apear_cluster_label", "paver_cluster_label",
  "enrich_plot_cluster_label"
)
expert_data <- all_result_long_data |>
  dplyr::filter(!method %in% tool_methods) |>
  dplyr::transmute(
    cluster = as.integer(cluster),
    expert_id = as.character(method),
    expert_label = as.character(combined_label),
    expert_embedding = combined_emb,
    involves_not_available = trimws(as.character(label)) == "Not Available"
  ) |>
  dplyr::arrange(cluster, expert_id)
if (any(lengths(expert_data$expert_embedding) == 0L)) {
  stop("At least one expert annotation is missing its cached embedding.")
}
if (length(unique(lengths(expert_data$expert_embedding))) != 1L) {
  stop("Expert embeddings have inconsistent dimensions.")
}
if (ncol(annotation_embedding_matrix) != lengths(expert_data$expert_embedding)[1L]) {
  stop(
    "New annotation embeddings and cached expert embeddings have different dimensions. ",
    "Use the existing text-embedding-3-small workflow."
  )
}

real_llm_data <- annotation_data |>
  dplyr::filter(group == "Real") |>
  dplyr::mutate(llm_cluster = as.integer(sub("Functional_module_", "", module)))
llm_expert_cosine_similarity <- real_llm_data |>
  dplyr::select(llm, module, llm_cluster, embedding_id) |>
  dplyr::inner_join(
    expert_data |>
      dplyr::rename(expert_cluster = cluster),
    by = c("llm_cluster" = "expert_cluster"),
    relationship = "many-to-many"
  ) |>
  dplyr::transmute(
    llm,
    module,
    cluster = llm_cluster,
    expert_id,
    expert_label,
    involves_not_available,
    cosine_similarity = mapply(
      function(embedding_id, expert_embedding) {
        cosine_similarity(
          annotation_embedding_matrix[embedding_id, ], expert_embedding
        )
      },
      embedding_id,
      expert_embedding
    )
  )

summarise_expert_consistency <- function(data, status) {
  data |>
    dplyr::group_by(llm) |>
    dplyr::summarise(
      status = status,
      n = dplyr::n(),
      mean = mean(cosine_similarity),
      sd = stats::sd(cosine_similarity),
      median = stats::median(cosine_similarity),
      q1 = as.numeric(stats::quantile(cosine_similarity, 0.25)),
      q3 = as.numeric(stats::quantile(cosine_similarity, 0.75)),
      .groups = "drop"
    )
}
llm_expert_consistency_summary <- dplyr::bind_rows(
  summarise_expert_consistency(llm_expert_cosine_similarity, "All expert annotations"),
  summarise_expert_consistency(
    dplyr::filter(llm_expert_cosine_similarity, !involves_not_available),
    "Excluding 'Not Available'"
  )
)

paired_expert_data <- llm_expert_cosine_similarity |>
  dplyr::select(llm, cluster, expert_id, cosine_similarity) |>
  tidyr::pivot_wider(names_from = llm, values_from = cosine_similarity)
expert_paired_wilcox <- stats::wilcox.test(
  paired_expert_data[["Qwen3.6-35B-A3B"]],
  paired_expert_data[["GPT-4o-mini"]],
  paired = TRUE,
  exact = FALSE
)
expert_panel_paired_wilcox <- paired_expert_data |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    paired_wilcox_p_value = stats::wilcox.test(
      `Qwen3.6-35B-A3B`, `GPT-4o-mini`, paired = TRUE, exact = FALSE
    )$p.value,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    significance = dplyr::case_when(
      paired_wilcox_p_value < 0.001 ~ "***",
      paired_wilcox_p_value < 0.01 ~ "**",
      paired_wilcox_p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

confidence_pairs <- annotation_data |>
  dplyr::select(llm, matched_module, group, confidence_score) |>
  tidyr::pivot_wider(names_from = group, values_from = confidence_score)
confidence_real_vs_random_summary <- confidence_pairs |>
  dplyr::group_by(llm) |>
  dplyr::group_modify(function(data, key) {
    paired_t <- stats::t.test(data$Real, data$Random, paired = TRUE)
    paired_wilcox <- stats::wilcox.test(
      data$Real, data$Random, paired = TRUE, exact = FALSE
    )
    paired_wilcox_symbol <- if (paired_wilcox$p.value < 0.001) {
      "***"
    } else if (paired_wilcox$p.value < 0.01) {
      "**"
    } else if (paired_wilcox$p.value < 0.05) {
      "*"
    } else {
      "ns"
    }
    rank_sum <- stats::wilcox.test(
      data$Real, data$Random, paired = FALSE, exact = FALSE
    )
    rank_sum_symbol <- if (rank_sum$p.value < 0.001) {
      "***"
    } else if (rank_sum$p.value < 0.01) {
      "**"
    } else if (rank_sum$p.value < 0.05) {
      "*"
    } else {
      "ns"
    }
    data.frame(
      n_pairs = nrow(data),
      real_mean = mean(data$Real),
      real_sd = stats::sd(data$Real),
      random_mean = mean(data$Random),
      random_sd = stats::sd(data$Random),
      mean_paired_difference = mean(data$Real - data$Random),
      paired_t_p_value = paired_t$p.value,
      paired_wilcox_p_value = paired_wilcox$p.value,
      paired_wilcox_significance = paired_wilcox_symbol,
      wilcoxon_rank_sum_p_value = rank_sum$p.value,
      rank_sum_significance = rank_sum_symbol
    )
  }) |>
  dplyr::ungroup()

write_result <- function(data, filename) {
  utils::write.csv(data, file.path(output_dir, filename), row.names = FALSE)
}
write_result(qwen_vs_gpt_semantic_similarity, "qwen_vs_gpt_semantic_similarity.csv")
write_result(qwen_vs_gpt_semantic_summary, "qwen_vs_gpt_semantic_summary.csv")
write_result(llm_expert_cosine_similarity, "llm_expert_cosine_similarity.csv")
write_result(llm_expert_consistency_summary, "llm_expert_consistency_summary.csv")
write_result(
  expert_panel_paired_wilcox,
  "llm_expert_paired_wilcoxon_by_module.csv"
)
write_result(annotation_data, "llm_annotation_and_confidence_scores.csv")
write_result(confidence_real_vs_random_summary, "confidence_real_vs_random_summary.csv")

annotation_evaluation <- list(
  qwen_vs_gpt_semantic_similarity = qwen_vs_gpt_semantic_similarity,
  qwen_vs_gpt_semantic_summary = qwen_vs_gpt_semantic_summary,
  llm_expert_cosine_similarity = llm_expert_cosine_similarity,
  llm_expert_consistency_summary = llm_expert_consistency_summary,
  expert_panel_paired_wilcox = expert_panel_paired_wilcox,
  confidence_scores = annotation_data,
  confidence_real_vs_random_summary = confidence_real_vs_random_summary,
  expert_paired_wilcox_p_value = expert_paired_wilcox$p.value,
  annotation_embedding_model = annotation_embedding_model
)
save(
  annotation_evaluation,
  file = file.path(output_dir, "annotation_evaluation.rda")
)

semantic_plot <- ggplot2::ggplot(
  qwen_vs_gpt_semantic_similarity,
  ggplot2::aes(x = group, y = cosine_similarity, fill = group)
) +
  ggplot2::geom_boxplot(width = 0.5, alpha = 0.7, outlier.shape = NA) +
  ggplot2::geom_jitter(
    width = 0.12, shape = 21, size = 2.7, colour = "black", alpha = 0.8
  ) +
  ggplot2::scale_fill_manual(values = c("Random" = "#084d68", "Real" = "#e53ba4")) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(
    x = NULL,
    y = "Qwen versus GPT cosine similarity"
  )
ggplot2::ggsave(
  file.path(output_dir, "qwen_vs_gpt_semantic_similarity.pdf"),
  semantic_plot, width = 6, height = 6
)

expert_plot_data <- llm_expert_cosine_similarity |>
  dplyr::mutate(
    llm = factor(llm, levels = c("GPT-4o-mini", "Qwen3.6-35B-A3B")),
    cluster = factor(cluster, levels = sort(unique(cluster)))
  )
expert_plot_stats <- expert_panel_paired_wilcox |>
  dplyr::mutate(
    cluster = factor(cluster, levels = levels(expert_plot_data$cluster)),
    x_start = 1,
    x_end = 2,
    x_label = 1.5,
    y_line = 0.94,
    y_label = 0.975
  )
expert_plot <- ggplot2::ggplot(
  expert_plot_data,
  ggplot2::aes(x = llm, y = cosine_similarity, fill = llm)
) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  ggplot2::geom_jitter(
    data = dplyr::filter(expert_plot_data, !involves_not_available),
    width = 0.12, shape = 21, size = 1.6, alpha = 0.75, colour = "black"
  ) +
  ggplot2::geom_point(
    data = dplyr::filter(expert_plot_data, involves_not_available),
    ggplot2::aes(
      x = llm,
      y = cosine_similarity,
      colour = "Involves expert 'Not Available' annotation"
    ),
    inherit.aes = FALSE,
    shape = 4, size = 2.5, stroke = 0.9
  ) +
  ggplot2::geom_segment(
    data = expert_plot_stats,
    ggplot2::aes(x = x_start, xend = x_end, y = y_line, yend = y_line),
    inherit.aes = FALSE,
    colour = "black", linewidth = 0.7
  ) +
  ggplot2::geom_text(
    data = expert_plot_stats,
    ggplot2::aes(x = x_label, y = y_label, label = significance),
    inherit.aes = FALSE,
    size = 4.2, fontface = "bold"
  ) +
  ggplot2::facet_wrap(
    ~cluster,
    labeller = ggplot2::labeller(cluster = function(x) paste("Functional Module", x)),
    nrow = 2
  ) +
  ggplot2::scale_fill_manual(
    values = c("GPT-4o-mini" = "#9ECAE1", "Qwen3.6-35B-A3B" = "#4C78A8")
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c("Involves expert 'Not Available' annotation" = "#D55E00")
  ) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 8)
  ) +
  ggplot2::labs(x = NULL, y = "Cosine similarity to expert annotation")
ggplot2::ggsave(
  file.path(output_dir, "llm_expert_consistency.pdf"),
  expert_plot, width = 9, height = 6
)

module_size_data <- load_one(
  file.path(
    project_root, "3_data_analysis", "4_random_module_dataset",
    "formated_result.Rdata"
  ),
  "formated_result"
) |>
  dplyr::transmute(
    module = as.character(module),
    module_size = as.numeric(module_content_number)
  )
confidence_plot_data <- annotation_data |>
  dplyr::filter(llm == "Qwen3.6-35B-A3B") |>
  dplyr::left_join(module_size_data, by = "module") |>
  dplyr::mutate(
    group = factor(group, levels = c("Random", "Real")),
    cluster_number = as.integer(sub("Functional_module_", "", matched_module))
  )
if (anyNA(confidence_plot_data$module_size)) {
  stop("At least one module size is missing from formated_result.Rdata.")
}
cluster_offsets <- stats::setNames(
  seq(-0.09, 0.09, length.out = length(unique(confidence_plot_data$cluster_number))),
  sort(unique(confidence_plot_data$cluster_number))
)
confidence_plot_data <- confidence_plot_data |>
  dplyr::mutate(
    x_base = ifelse(group == "Random", 1, 2),
    x_jittered = x_base + unname(cluster_offsets[as.character(cluster_number)])
  )
qwen_paired_wilcox <- confidence_real_vs_random_summary |>
  dplyr::filter(llm == "Qwen3.6-35B-A3B")
confidence_plot <- ggplot2::ggplot(
  confidence_plot_data,
  ggplot2::aes(y = confidence_score)
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(x = x_base, fill = group, group = group),
    width = 0.45, alpha = 0.5, outlier.shape = NA,
    colour = "#4D4D4D", linewidth = 0.8
  ) +
  ggplot2::geom_line(
    ggplot2::aes(x = x_jittered, group = matched_module),
    colour = "#4D4D4D", alpha = 0.75, linewidth = 0.55
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = x_jittered, fill = group, size = module_size
    ),
    shape = 21, colour = "#333333", stroke = 0.7, alpha = 0.95
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(
      x = x_jittered, label = cluster_number
    ),
    seed = analysis_seed,
    size = 4,
    colour = "black",
    box.padding = 0.18,
    point.padding = 0.12,
    min.segment.length = Inf,
    segment.colour = NA,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::annotate(
    "segment", x = 1, xend = 2, y = 1.012, yend = 1.012,
    colour = "black", linewidth = 0.8
  ) +
  ggplot2::annotate(
    "text", x = 1.5, y = 1.022,
    label = qwen_paired_wilcox$paired_wilcox_significance,
    size = 6, fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(values = c("Random" = "#084d68", "Real" = "#e53ba4")) +
  ggplot2::scale_size_continuous(
    name = "Module size", limits = c(3, 7), breaks = 3:7, range = c(2.8, 6.5)
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(1, 2), labels = c("Random", "Real"), limits = c(0.65, 2.35)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0.6, 1.0, 0.1),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::coord_cartesian(ylim = c(0.53, 1.035), clip = "off") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(size = 13),
    legend.position = c(0.82, 0.25),
    legend.title = ggplot2::element_text(face = "bold"),
    legend.background = ggplot2::element_rect(fill = "transparent", colour = NA)
  ) +
  ggplot2::guides(fill = "none") +
  ggplot2::labs(x = NULL, y = "Confidence Score")
ggplot2::ggsave(
  file.path(output_dir, "confidence_real_vs_random.pdf"),
  confidence_plot, width = 7, height = 6
)

qwen_semantic_overall <- dplyr::filter(
  qwen_vs_gpt_semantic_summary, as.character(group) == "Overall"
)
qwen_expert_reporting <- dplyr::filter(
  llm_expert_consistency_summary,
  llm == "Qwen3.6-35B-A3B", status == "All expert annotations"
)
qwen_confidence <- dplyr::filter(
  confidence_real_vs_random_summary, llm == "Qwen3.6-35B-A3B"
)
summary_lines <- c(
  sprintf("Qwen LLM model: %s", qwen_llm_model),
  sprintf("Semantic embedding model: %s", annotation_embedding_model),
  sprintf(
    "Qwen versus GPT semantic cosine similarity: median %.3f (IQR %.3f-%.3f), n=%d",
    qwen_semantic_overall$median, qwen_semantic_overall$q1,
    qwen_semantic_overall$q3, qwen_semantic_overall$n
  ),
  sprintf(
    "Qwen versus experts: median %.3f (IQR %.3f-%.3f), n=%d; includes 'Not Available'",
    qwen_expert_reporting$median, qwen_expert_reporting$q1,
    qwen_expert_reporting$q3, qwen_expert_reporting$n
  ),
  sprintf(
    "Paired Qwen versus GPT expert-consistency Wilcoxon p-value: %.6g",
    expert_paired_wilcox$p.value
  ),
  sprintf(
    "Qwen confidence (real): %.3f +/- %.3f; random: %.3f +/- %.3f",
    qwen_confidence$real_mean, qwen_confidence$real_sd,
    qwen_confidence$random_mean, qwen_confidence$random_sd
  ),
  sprintf(
    "Qwen real-random paired difference: %.3f; paired t p=%.6g; paired Wilcoxon p=%.6g",
    qwen_confidence$mean_paired_difference, qwen_confidence$paired_t_p_value,
    qwen_confidence$paired_wilcox_p_value
  )
)
writeLines(summary_lines, file.path(output_dir, "annotation_evaluation_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
