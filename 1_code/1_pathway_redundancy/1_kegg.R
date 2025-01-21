library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')


dir.create(
  "3_data_analysis/1_pathway_redundancy/1_kegg",
  showWarnings = FALSE,
  recursive = TRUE
)

setwd("3_data_analysis/1_pathway_redundancy/1_kegg")

library(KEGGgraph)
library(KEGGREST)
library(tidyverse)

load("pathway.rda")

hsa_pathway <-
  pathway
hsa_pathway@compound_list[[10]]$KEGG.ID

pathway <-
  hsa_pathway@compound_list %>%
  lapply(function(x)
    x$KEGG.ID)

names(pathway) <-
  hsa_pathway@pathway_id

# Create a binary incidence matrix
unique_metabolites <- unique(unlist(pathway))
incidence_matrix <- sapply(pathway, function(pathway) {
  as.integer(unique_metabolites %in% pathway)
})

# Compute intersection and union using matrix multiplication
intersection_matrix <- t(incidence_matrix) %*% incidence_matrix
union_matrix <- outer(rowSums(incidence_matrix), rowSums(incidence_matrix), FUN = "+") - intersection_matrix

# Compute Jaccard index
jaccard_matrix <- intersection_matrix / union_matrix

# Display the result
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(kegg_reference_pathway)
jaccard_matrix
