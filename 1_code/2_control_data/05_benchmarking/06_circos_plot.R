library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

setwd("3_data_analysis/02_control_data/05_benchmarking")

load("comparison_result/all_clustering_res.rda")

# Load required libraries
library(circlize)
library(dplyr)
library(RColorBrewer)

df <- all_clustering_res

# Function to create chord diagram from ground truth to others
create_chord_diagram <- function(df) {
  # Create adjacency matrix only from ground truth to others
  methods <- c("ground_truth_label", "mapa_label", "paver_label", "apear_cluster", "enrichplot_cluster")
  method_names <- c("Ground_Truth", "MAPA", "PAVER", "aPEAR", "enrichplot")

  # Get unique clusters for ground truth
  gt_clusters <- sort(unique(df$ground_truth_label))

  # Create a list to store matrices for each comparison
  all_data <- data.frame()

  # Create data for chord diagram from ground truth to each other method
  for(i in 2:length(methods)) {
    method <- methods[i]
    method_name <- method_names[i]

    # Create contingency table
    connections <- df %>%
      group_by(ground_truth_label, !!sym(method)) %>%
      summarise(count = n(), .groups = 'drop')

    # Add method labels
    connections$from <- paste0("GT_", connections$ground_truth_label)
    connections$to <- paste0(method_name, "_", connections[[method]])

    all_data <- rbind(all_data,
                      data.frame(from = connections$from,
                                 to = connections$to,
                                 value = connections$count))
  }

  # Convert to matrix format for chordDiagram
  # Get all unique sectors
  all_sectors <- unique(c(all_data$from, all_data$to))

  # Create adjacency matrix
  adj_matrix <- matrix(0, nrow = length(all_sectors), ncol = length(all_sectors))
  rownames(adj_matrix) <- all_sectors
  colnames(adj_matrix) <- all_sectors

  # Fill the matrix
  for(i in 1:nrow(all_data)) {
    adj_matrix[all_data$from[i], all_data$to[i]] <- all_data$value[i]
  }

  # Create colors for sectors
  sector_colors <- character(length(all_sectors))
  names(sector_colors) <- all_sectors

  # Color ground truth sectors in grey
  gt_sectors <- grep("^GT_", all_sectors, value = TRUE)
  sector_colors[gt_sectors] <- "grey"

  # Color other method sectors
  mapa_sectors <- grep("^MAPA", all_sectors, value = TRUE)
  if(length(mapa_sectors) > 0) sector_colors[mapa_sectors] <- "#DD5129FF"

  pave_sectors <- grep("^PAVER", all_sectors, value = TRUE)
  if(length(pave_sectors) > 0) sector_colors[pave_sectors] <- "#43B284FF"

  apea_sectors <- grep("^aPEAR", all_sectors, value = TRUE)
  if(length(apea_sectors) > 0) sector_colors[apea_sectors] <- "#0F7BA2FF"

  enri_sectors <- grep("^enrichplot", all_sectors, value = TRUE)
  if(length(enri_sectors) > 0) sector_colors[enri_sectors] <- "#FAB255FF"

  # Create link colors based on destination method
  link_colors <- matrix("grey", nrow = nrow(adj_matrix), ncol = ncol(adj_matrix))
  rownames(link_colors) <- rownames(adj_matrix)
  colnames(link_colors) <- colnames(adj_matrix)

  # Set colors for links going to each method
  for(i in 1:nrow(adj_matrix)) {
    for(j in 1:ncol(adj_matrix)) {
      if(adj_matrix[i,j] > 0) {
        col_name <- colnames(adj_matrix)[j]
        if(grepl("^MAPA", col_name)) {
          link_colors[i,j] <- adjustcolor("#DD5129FF", alpha.f = 0.6)
        } else if(grepl("^PAVER", col_name)) {
          link_colors[i,j] <- adjustcolor("#43B284FF", alpha.f = 0.6)
        } else if(grepl("^aPEAR", col_name)) {
          link_colors[i,j] <- adjustcolor("#0F7BA2FF", alpha.f = 0.6)
        } else if(grepl("^enrichplot", col_name)) {
          link_colors[i,j] <- adjustcolor("#FAB255FF", alpha.f = 0.6)
        }
      }
    }
  }

  # Create chord diagram
  circos.clear()
  chordDiagram(adj_matrix,
               grid.col = sector_colors,
               col = link_colors,
               transparency = 0,
               annotationTrack = "grid",
               preAllocateTracks = 1,
               directional = 1,  # Make it directional from ground truth
               direction.type = "arrows",
               link.arr.type = "big.arrow")

  # Add labels
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
  }, bg.border = NA)

  # Add legend
  legend_labels <- c("Ground Truth", "MAPA", "PAVER", "aPEAR", "enrichplot")
  legend_colors <- c("grey", "#DD5129FF", "#43B284FF", "#0F7BA2FF", "#FAB255FF")

  # Add legend to the plot
  legend("bottomright",
         legend = legend_labels,
         fill = legend_colors,
         cex = 0.8,
         title = "Methods",
         title.cex = 0.9,
         bg = "white",
         box.col = "black")
}

# Create the chord diagram
create_chord_diagram(df)

# Reset circos parameters
circos.clear()

pdf("comparison_result/chord_diagram.pdf", width = 8, height = 6)
create_chord_diagram(df)
dev.off()
