library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

setwd("2_data/case_study/02_monkey_multi_omics/")

tissues <- c("aortic_arch", "kidney", "ovary", "spleen",
             "stomach", "thymus", "thyroid")
directions <- c("up", "down")

split_module_ids <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(x))) {
    return(character())
  }

  ids <- trimws(strsplit(as.character(x), "/", fixed = TRUE)[[1]])
  unique(ids[nzchar(ids)])
}

load_input_gene_symbols <- function(rda_path) {
  object_env <- new.env(parent = emptyenv())
  object_names <- load(rda_path, envir = object_env)
  if (length(object_names) != 1) {
    stop("Expected one object in ", rda_path)
  }

  variable_info <- methods::slot(object_env[[object_names[[1]]]], "variable_info")
  symbols <- trimws(as.character(variable_info$symbol))
  unique(symbols[!is.na(symbols) & nzchar(symbols)])
}

# Collect feature counts for metabolite-containing modules with
# module_content_number >= 3 and true_omics_num == 3 ----
all_data <- lapply(tissues, function(tissue) {
  lapply(directions, function(dir) {
    result_dir <- file.path("taxid_9606", tissue, dir)
    rda_path <- file.path(result_dir, "multi_omics_fm.rda")
    rna_path <- file.path(result_dir, paste0(dir, "_rna_enriched_pathways.rda"))
    prot_path <- file.path(result_dir, paste0(dir, "_prot_enriched_pathways.rda"))
    if (!all(file.exists(c(rda_path, rna_path, prot_path)))) return(NULL)

    result_env <- new.env(parent = emptyenv())
    load(rda_path, envir = result_env)
    fm_result <- result_env$multi_omics_fm$functional_module_result

    rna_symbols <- load_input_gene_symbols(rna_path)
    prot_symbols <- load_input_gene_symbols(prot_path)

    target_module_info <- fm_result |>
      dplyr::filter(
        module_content_number >= 3,
        true_omics_num == 3,
        include_metabolites
      ) |>
      dplyr::select(module, module_content_number, genes, metabolites)
    if (nrow(target_module_info) == 0) return(NULL)

    purrr::map_dfr(seq_len(nrow(target_module_info)), function(i) {
      genes <- split_module_ids(target_module_info$genes[[i]])
      metabolites <- split_module_ids(target_module_info$metabolites[[i]])

      tibble::tibble(
        module = target_module_info$module[[i]],
        module_content_number = target_module_info$module_content_number[[i]],
        tissue = tissue,
        direction = dir,
        Metabolome = length(metabolites),
        Proteome = sum(genes %in% prot_symbols),
        Transcriptome = sum(genes %in% rna_symbols)
      )
    })
  })
})

# Collect confidence scores from llm_interpreted_object ----
conf_data <- lapply(tissues, function(tissue) {
  lapply(directions, function(dir) {
    llm_path <- file.path("taxid_9606", tissue, dir, "llm_interpreted_object.rda")
    if (!file.exists(llm_path)) return(NULL)

    llm_env <- new.env(parent = emptyenv())
    load(llm_path, envir = llm_env)
    interp <- llm_env$llm_interpreted_object$llm_module_interpretation
    if (is.null(interp) || length(interp) == 0) return(NULL)

    scores <- lapply(names(interp), function(mod_id) {
      cs <- interp[[mod_id]]$generated_name$confidence_score
      if (is.null(cs)) return(NULL)
      data.frame(module = mod_id, tissue = tissue, direction = dir,
                 confidence_score = as.numeric(cs),
                 stringsAsFactors = FALSE)
    })
    bind_rows(Filter(Negate(is.null), scores))
  })
})


conf_df <- bind_rows(Filter(Negate(is.null), unlist(conf_data, recursive = FALSE)))

module_summary <- bind_rows(Filter(Negate(is.null), unlist(all_data, recursive = FALSE))) |>
  mutate(
    total_count = Metabolome + Proteome + Transcriptome,
    y_label = paste0(gsub("_", " ", tissue), " (", direction, ") | ", module)
  ) |>
  left_join(conf_df, by = c("module", "tissue", "direction")) |>
  arrange(factor(tissue, levels = tissues),
          desc(module_content_number), direction, module)

if (nrow(module_summary) != 21) {
  stop("Expected 21 target modules, found ", nrow(module_summary))
}

table_data <- module_summary |>
  transmute(
    y_label,
    module,
    tissue,
    regulation_direction = direction,
    total_count,
    Metabolome,
    Proteome,
    Transcriptome,
    confidence_score
  )

write.csv(
  table_data,
  file = "taxid_9606/Supplementary_Data_2_Table_S4_multi_omics_modules.csv",
  row.names = FALSE,
  na = ""
)

plot_data <- module_summary |>
  tidyr::pivot_longer(
    cols = c(Metabolome, Proteome, Transcriptome),
    names_to = "omic",
    values_to = "count"
  ) |>
  mutate(
    omic    = factor(omic,
                     levels = c("Transcriptome", "Proteome", "Metabolome"))
  )

# Order y-axis: tissues alphabetically, then within each tissue sort by
# module_content_number descending (largest at top, smallest at bottom)
y_order <- module_summary$y_label
plot_data$y_label <- factor(plot_data$y_label, levels = rev(y_order))

# Colors for omic types
omic_colors <- c(
  "Transcriptome" = "#fae69e",
  "Proteome"      = "#f2b56f",
  "Metabolome"    = "#71b7ed"
)

plot_data_filtered <- plot_data

# Per-module total bar length and confidence score for point layer ----
point_data <- plot_data_filtered |>
  group_by(y_label, module, tissue, direction) |>
  summarise(total_count = sum(count), .groups = "drop") |>
  left_join(conf_df, by = c("module", "tissue", "direction"))

# Scale factor: map confidence_score [0,1] onto the primary x-axis range
max_count <- max(point_data$total_count, na.rm = TRUE)
scale_fac <- max_count

point_data <- point_data |>
  mutate(cs_scaled = confidence_score * scale_fac)

p <- ggplot(plot_data_filtered,
            aes(x = count, y = y_label, fill = omic)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.25) +
  geom_point(
    data = point_data |> filter(!is.na(cs_scaled)),
    aes(x = cs_scaled, y = y_label),
    inherit.aes = FALSE,
    shape = 23, fill = "grey", color = "black", size = 1.5, stroke = 0.5
  ) +
  geom_vline(xintercept = 0.6 * scale_fac, color = "grey", linetype = 2) +
  scale_fill_manual(values = omic_colors, name = "Omics") +
  scale_x_continuous(
    name = "Feature number",
    sec.axis = sec_axis(
      transform = ~ . / scale_fac,
      name   = "Confidence score",
      breaks = seq(0, 1, by = 0.2)
    )
  ) +
  labs(y = "Module") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y        = element_text(size = 8),
    axis.title.x.top   = element_text(size = 10),
    axis.text.x.top    = element_text(size = 8),
    legend.position    = "left"
  )

p

ggsave(
  filename = file.path("taxid_9606/multi_omics_tissue_feature_count_updated.pdf"),
  plot     = p,
  width    = 5,
  height   = 6.8
)

plot_data_filtered |>
  dplyr::distinct(y_label, module, tissue, direction, module_content_number) |>
  dplyr::arrange(factor(tissue, levels = tissues),
                 dplyr::desc(module_content_number), direction, module)
