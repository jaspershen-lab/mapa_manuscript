library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/all_gene_list.rda")

library(VennDiagram)
library(grid)

A <- all_gene_list[["hsa04724"]]
B <- all_gene_list[["GO:0098978"]]
C <- all_gene_list[["GO:0098985"]]

# Calculate overlaps
area1 <- length(A)
area2 <- length(B)
area3 <- length(C)
n12   <- length(intersect(A, B))
n23   <- length(intersect(B, C))
n13   <- length(intersect(A, C))
n123  <- length(Reduce(intersect, list(A, B, C)))

# Draw the triple Venn diagram
venn_plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12   = n12,
  n23   = n23,
  n13   = n13,
  n123  = n123,
  # category = c("hsa03030 \nDNA replication", "R-HSA-68952 \nDNA replication initiation", "R-HSA-69306 \nDNA replication"),

  # Specify fill colors and transparency
  fill  = c("#4f94d7", "#d15356", "#ebb869"),
  alpha = c(0.5, 0.5, 0.5),

  # Specify outline (border) colors, line style, and line width
  col   = c("#4f94d7", "#d15356", "#ebb869"),
  lty   = "solid",
  lwd   = 2
)

# Render the diagram
grid.newpage()
grid.draw(venn_plot)

ggsave(plot = venn_plot,
       filename = "3_data_analysis/02_control_data/03_venn_plot_glutamatergic_synapse.pdf",
       width = 6,
       height = 6)
