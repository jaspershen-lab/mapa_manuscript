library(r4projects)
library(mapa)

setwd(r4projects::get_project_wd())
rm(list = ls())

data_dir <- file.path(
  "2_data", "case_study", "02_monkey_multi_omics", "taxid_9606"
)
output_dir <- file.path(
  data_dir,
  "articulation_point_metabolite_multi_omics_networks"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

topology_csv <- file.path(
  data_dir,
  paste0(
    "metabolite_containing_multi_omics_module_",
    "metabolite_degree_articulation.csv"
  )
)
topology_xlsx <- sub("\\.csv$", ".xlsx", topology_csv)

if (file.exists(topology_xlsx)) {
  metabolite_topology <- as.data.frame(
    readxl::read_excel(topology_xlsx),
    stringsAsFactors = FALSE
  )
} else if (file.exists(topology_csv)) {
  metabolite_topology <- read.csv(
    topology_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
} else {
  stop(
    "Could not find the metabolite topology table as either XLSX or CSV"
  )
}

required_columns <- c(
  "Tissue", "Direction", "Module", "module_content_number",
  "metabolite_id", "metabolite_name",
  "is_articulation_point_within_module"
)
missing_columns <- setdiff(required_columns, names(metabolite_topology))
if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in metabolite topology table: ",
    paste(missing_columns, collapse = ", ")
  )
}

is_articulation <- metabolite_topology[[
  "is_articulation_point_within_module"
]]
if (is.logical(is_articulation)) {
  is_articulation[is.na(is_articulation)] <- FALSE
} else {
  is_articulation <- toupper(trimws(as.character(is_articulation))) == "TRUE"
  is_articulation[is.na(is_articulation)] <- FALSE
}

articulation_metabolites <- metabolite_topology[
  is_articulation,
  ,
  drop = FALSE
]
if (nrow(articulation_metabolites) == 0) {
  stop("No articulation-point metabolites were found")
}

module_keys <- unique(
  articulation_metabolites[c("Tissue", "Direction", "Module")]
)
module_keys <- module_keys[
  order(module_keys$Tissue, module_keys$Direction, module_keys$Module),
  ,
  drop = FALSE
]
row.names(module_keys) <- NULL

if (nrow(module_keys) != 13) {
  stop(
    "Expected 13 modules containing articulation-point metabolites, found ",
    nrow(module_keys)
  )
}

highlight_colour <- "#D62728"
node_size <- 5
manifest <- vector("list", nrow(module_keys))

for (i in seq_len(nrow(module_keys))) {
  tissue <- module_keys$Tissue[[i]]
  direction <- module_keys$Direction[[i]]
  module_id <- module_keys$Module[[i]]

  current_articulation <- articulation_metabolites[
    articulation_metabolites$Tissue == tissue &
      articulation_metabolites$Direction == direction &
      articulation_metabolites$Module == module_id,
    ,
    drop = FALSE
  ]
  articulation_ids <- sort(unique(current_articulation$metabolite_id))

  rda_file <- file.path(
    data_dir, tissue, direction, "multi_omics_fm.rda"
  )
  if (!file.exists(rda_file)) {
    stop("Missing MAPA result file: ", rda_file)
  }

  object_env <- new.env(parent = emptyenv())
  object_names <- load(rda_file, envir = object_env)
  if (!"multi_omics_fm" %in% object_names) {
    stop("Object 'multi_omics_fm' was not found in: ", rda_file)
  }
  multi_omics_fm <- object_env$multi_omics_fm

  module_nodes <- multi_omics_fm$result_with_module[
    multi_omics_fm$result_with_module$module == module_id,
    ,
    drop = FALSE
  ]
  if (nrow(module_nodes) == 0) {
    stop("Module was not found in MAPA result: ", module_id)
  }
  if (!all(articulation_ids %in% module_nodes$node_id)) {
    stop(
      "At least one articulation-point metabolite is absent from ",
      tissue, "/", direction, "/", module_id
    )
  }
  highlighted_nodes <- module_nodes[
    module_nodes$node_id %in% articulation_ids &
      module_nodes$node_type == "metabolite",
    ,
    drop = FALSE
  ]
  if (nrow(highlighted_nodes) != length(articulation_ids)) {
    stop(
      "Highlighted node identity/type check failed for ",
      tissue, "/", direction, "/", module_id
    )
  }

  set.seed(1000 + i)
  network_plot <- mapa::plot_multi_omics_module_info(
    merge_result = multi_omics_fm,
    module_id = module_id,
    node_size = node_size,
    label_size = 2.6,
    show_rwr_edge = TRUE,
    show_labels = TRUE,
    title = paste(tissue, direction, module_id, sep = " | ")
  )

  # MAPA uses black node outlines. Add a transparent overlay so that only
  # articulation-point metabolites receive a clearly visible red outline,
  # while retaining MAPA's original metabolite fill colour.
  network_plot <- network_plot +
    ggraph::geom_node_point(
      data = function(layout_data) {
        layout_data[
          layout_data$node_type == "metabolite" &
            layout_data$node_id %in% articulation_ids,
          ,
          drop = FALSE
        ]
      },
      shape = 24,
      size = node_size,
      fill = NA,
      colour = highlight_colour,
      stroke = 1.4,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      caption = paste0(
        "Red outline: articulation-point metabolite",
        if (length(articulation_ids) > 1) "s" else "",
        " (", paste(articulation_ids, collapse = ", "), ")"
      )
    ) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        colour = highlight_colour,
        face = "bold",
        hjust = 0,
        size = 9,
        margin = ggplot2::margin(t = 8)
      )
    )

  module_size <- nrow(module_nodes)
  plot_width <- max(7, min(10, 6 + sqrt(module_size) / 2))
  plot_height <- max(6, min(8.5, plot_width * 0.82))
  pdf_name <- paste(
    tissue,
    direction,
    module_id,
    "articulation_metabolite_network.pdf",
    sep = "_"
  )
  pdf_file <- file.path(output_dir, pdf_name)

  ggplot2::ggsave(
    filename = pdf_file,
    plot = network_plot,
    device = grDevices::cairo_pdf,
    width = plot_width,
    height = plot_height,
    units = "in",
    bg = "white"
  )

  manifest[[i]] <- data.frame(
    Tissue = tissue,
    Direction = direction,
    Module = module_id,
    module_content_number = module_size,
    articulation_metabolite_count = length(articulation_ids),
    articulation_metabolite_ids = paste(articulation_ids, collapse = "/"),
    pdf_file = pdf_name,
    stringsAsFactors = FALSE
  )

  message(
    "Wrote ", i, "/", nrow(module_keys), ": ", pdf_file
  )
}

manifest <- do.call(rbind, manifest)
write.csv(
  manifest,
  file = file.path(output_dir, "network_manifest.csv"),
  row.names = FALSE,
  na = ""
)

message(
  "Completed ", nrow(manifest), " MAPA module-network PDFs in: ",
  output_dir
)
