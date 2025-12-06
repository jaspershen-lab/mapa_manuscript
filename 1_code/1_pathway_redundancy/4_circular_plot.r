library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(circlize)

load("2_data/pathway_redundancy/all_pairs.rda")

# Filter data for sim >= 0.5
filtered_data <- all_pairs %>%
  filter(sim >= 0.5)

# Get unique pair types ordered by group
pair_types <- unique(filtered_data$pair_type)

# Create histogram data for each pair type
hist_data <- lapply(pair_types, function(pt) {
  sims <- filtered_data %>%
    filter(pair_type == pt) %>%
    pull(sim)

  # Create histogram with breaks from 0.5 to 1.0
  h <- hist(sims, breaks = seq(0.5, 1.0, by = 0.01), plot = FALSE)

  data.frame(
    pair_type = pt,
    bin_start = h$breaks[-length(h$breaks)],
    bin_end = h$breaks[-1],
    freq = h$counts,
    log_freq = log2(h$counts + 1)  # Add 1 to avoid log(0)
  )
})

hist_data <- do.call(rbind, hist_data)

# Convert pair_type to factor with the desired order
pair_type_order <- c("BP_BP", "MF_MF", "CC_CC",
                     "KEGG_KEGG", "Reactome_Reactome",
                     "BP_KEGG", "BP_Reactome", "KEGG_Reactome")

hist_data$pair_type <- factor(hist_data$pair_type, levels = pair_type_order)

# Define colors for each sector
sector_colors <- c("#503e80", "#e66e53", "#efa361", "#ebc465",
                   "#8ab17c", "#279e8d", "#23756c", "#264655")

names(sector_colors) <- pair_type_order

# Initialize circos plot =====
circos.clear()
gap_sizes <- rep(0, length(pair_type_order))
gap_sizes[length(gap_sizes)] <- 4

circos.par(start.degree = 90, gap.degree = gap_sizes)

# Initialize sectors
circos.initialize(factors = hist_data$pair_type,
                  xlim = c(0.5, 1.0))

# Create track with histograms
circos.track(
  factors = hist_data$pair_type,
  ylim = c(0, max(hist_data$log_freq, na.rm = TRUE)),
  panel.fun = function(x, y) {
    sector_data <- hist_data[hist_data$pair_type == CELL_META$sector.index, ]

    sector_col <- sector_colors[CELL_META$sector.index]

    # Draw histogram bars
    for (i in 1:nrow(sector_data)) {
      circos.rect(
        xleft = sector_data$bin_start[i],
        xright = sector_data$bin_end[i],
        ybottom = 0,
        ytop = sector_data$log_freq[i],
        col = sector_col,
        border = NA
        # lwd = 0.2
      )
    }

    # Add axes
    circos.axis(h = "top", major.at = seq(0.5, 1.0, by = 0.1),
                labels.cex = 0.5)
    # circos.yaxis(side = "left", labels.cex = 0.5)
    # Add y-axis only for the first sector
    if (CELL_META$sector.index == pair_type_order[1]) {
      circos.yaxis(side = "left", labels.cex = 0.5)
    }

  },
  bg.border = NA,
  track.height = 0.6
)

# Add sector label
circos.track(
  factors = hist_data$pair_type,
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    circos.text(
      x = mean(CELL_META$xlim),
      y = 1,
      labels = CELL_META$sector.index,
      facing = "bending.outside",
      niceFacing = FALSE,
      cex = 0.6,
      font = 2
    )
  },
  bg.border = NA,
  track.height = 0.1
)

# Clear circos parameters
circos.clear()
