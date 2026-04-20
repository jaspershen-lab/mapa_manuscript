library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

setwd("2_data/case_study/02_monkey_multi_omics/")

tissues <- c("aortic_arch", "thymus", "spleen", "thyroid")
directions <- c("up", "down")

# Collect feature counts for modules with multi_omics_num == 3 ----
all_data <- lapply(tissues, function(tissue) {
  lapply(directions, function(dir) {
    rda_path <- file.path("taxid_9606", tissue, dir, "multi_omics_fm.rda")
    if (!file.exists(rda_path)) return(NULL)

    load(rda_path)

    fm_result <- multi_omics_fm$functional_module_result
    result_node <- multi_omics_fm$result_with_module

    m3_modules <- fm_result$module[fm_result$multi_omics_num == 3]
    if (length(m3_modules) == 0) return(NULL)

    # For each node in m3 modules, extract omic source:
    #   - gene with dt_src containing "T" -> Transcriptome
    #   - gene with dt_src containing "P" -> Proteome
    #   - metabolite                       -> Metabolome
    #   - pathway nodes are excluded (no single omic source)
    m3_nodes <- result_node %>% filter(module %in% m3_modules)

    omic_rows <- lapply(seq_len(nrow(m3_nodes)), function(i) {
      row       <- m3_nodes[i, ]
      node_info <- row$node_info[[1]]

      if (row$node_type == "metabolite") {
        return(data.frame(module = row$module, omic = "Metabolome",
                          stringsAsFactors = FALSE))
      }
      if (row$node_type == "gene" && is.list(node_info) && !is.null(node_info$dt_src)) {
        src <- node_info$dt_src
        out <- data.frame(module = character(), omic = character(),
                          stringsAsFactors = FALSE)
        if (grepl("T", src))
          out <- rbind(out, data.frame(module = row$module, omic = "Transcriptome",
                                       stringsAsFactors = FALSE))
        if (grepl("P", src))
          out <- rbind(out, data.frame(module = row$module, omic = "Proteome",
                                       stringsAsFactors = FALSE))
        return(out)
      }
      NULL  # pathway nodes excluded
    })

    counts <- bind_rows(omic_rows) %>%
      group_by(module, omic) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(tissue = tissue, direction = dir)
    counts
  })
})

# Collect confidence scores from llm_interpreted_object ----
conf_data <- lapply(tissues, function(tissue) {
  lapply(directions, function(dir) {
    llm_path <- file.path("taxid_9606", tissue, dir, "llm_interpreted_object.rda")
    if (!file.exists(llm_path)) return(NULL)

    load(llm_path)  # loads llm_interpreted_object

    interp <- llm_interpreted_object$llm_module_interpretation
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

plot_data <- bind_rows(Filter(Negate(is.null), unlist(all_data, recursive = FALSE)))

# Build a y-axis label: "tissue (direction) | module_id"
plot_data <- plot_data %>%
  mutate(
    y_label = paste0(gsub("_", " ", tissue), " (", direction, ") | ", module),
    omic    = factor(omic,
                     levels = c("Transcriptome", "Proteome", "Metabolome"))
  )

# Order y-axis: tissue order as defined in `tissues`, then within each tissue
# sort by total feature count descending (largest at top, smallest at bottom)
total_counts <- plot_data %>%
  group_by(y_label) %>%
  summarise(total_count = sum(count), .groups = "drop")

plot_data <- plot_data %>%
  mutate(tissue_f = factor(tissue, levels = tissues)) %>%
  left_join(total_counts, by = "y_label") %>%
  arrange(tissue_f, desc(total_count))

y_order <- unique(plot_data$y_label)
plot_data$y_label <- factor(plot_data$y_label, levels = rev(y_order))

# Colors for omic types
omic_colors <- c(
  "Transcriptome" = "#fae69e",
  "Proteome"      = "#f2b56f",
  "Metabolome"    = "#71b7ed"
)

multi_omics_fm_id <- plot_data |>
  group_by(y_label) |>
  summarise(freq = n()) |>
  ungroup() |>
  filter(freq > 1) |>
  pull(y_label) |>
  unique()

plot_data_filtered <- plot_data |> filter(y_label %in% multi_omics_fm_id)

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
  filename = file.path("taxid_9606/multi_omics_tissue_feature_count.pdf"),
  plot     = p,
  width    = 5,
  height   = 7.2
)



plot_data_filtered |>
  dplyr::distinct(module, tissue) |>
  dplyr::group_by(tissue) |>
  summarise(freq = n()) |>
  ungroup()
