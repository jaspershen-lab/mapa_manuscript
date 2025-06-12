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

pathway <-
  hsa_pathway@gene_list %>%
  lapply(function(x)
    x$KEGG.ID)

names(pathway) <-
  hsa_pathway@pathway_id

##remove NULL pathways
remove_index <-
  lapply(pathway, function(x) {
    is.null(x)
  }) %>%
  unlist() %>%
  which()

if (length(remove_index) > 0) {
  pathway <- pathway[-remove_index]
}

# Function to calculate Jaccard index
jaccard_index <- function(a, b) {
  intersect_len <- length(intersect(a, b))
  union_len <- length(union(a, b))
  if (union_len == 0)
    return(0)
  return(intersect_len / union_len)
}

# Get pathway names
pathway_names <- names(pathway)

# Initialize matrix
n <- length(pathway_names)
jaccard_matrix <- matrix(0, nrow = n, ncol = n)
rownames(jaccard_matrix) <- pathway_names
colnames(jaccard_matrix) <- pathway_names

# Compute pairwise Jaccard indices
for (i in 1:n) {
  cat(i, " ")
  for (j in i:n) {
    ji <- jaccard_index(pathway[[i]], pathway[[j]])
    jaccard_matrix[i, j] <- ji
    jaccard_matrix[j, i] <- ji  # symmetric
  }
}

jaccard_matrix[upper.tri(jaccard_matrix)] %>%
  data.frame(jaccard_index = .) %>%
  ggplot() +
  geom_histogram(aes(x = jaccard_index), bins = 50)


