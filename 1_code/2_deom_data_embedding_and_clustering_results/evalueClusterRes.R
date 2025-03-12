# Jaccard
# Create example data (replace with your actual dataframe)
pathway_pairs <- jaccard_index

# Get all unique pathway names
for (i in 1:nrow(pathway_pairs)) {
  new_sub_row <- data.frame("name1" = pathway_pairs$name2[i], "name2" = pathway_pairs$name1[i], "value" = pathway_pairs$value[i])
  pathway_pairs <- rbind(pathway_pairs, new_sub_row)
}

# Convert to a wide matrix
sim_matrix <- pathway_pairs %>%
  tidyr::pivot_wider(
    names_from = name1,
    values_from = value,
    values_fill = NA  # Missing pairs will be NA
  ) %>%
  tibble::column_to_rownames("name2") %>%
  as.matrix()
sim_matrix <- sim_matrix[rownames(sim_matrix), rownames(sim_matrix)]
diag(sim_matrix) <- 1



# Use cosine similarity ====
# p.adjust.cutoff.go = 0.05
# p.adjust.cutoff.kegg = 0.05
# p.adjust.cutoff.reactome = 0.05
# count.cutoff.go = 5
# count.cutoff.kegg = 5
# count.cutoff.reactome = 5
# sim.cutoff = 0.5

## Create tidygraph object
### Collect edge data
edge_data <-
  sim_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "name1") %>%
  tidyr::pivot_longer(
    cols = -name1,
    names_to = "name2",
    values_to = "sim"
  ) %>%
  dplyr::filter(name1 != name2) %>%  # exclude self-edges
  dplyr::filter(name1 > name2) %>%   # exclude duplicates
  dplyr::rename(from = name1, to = name2) %>%
  dplyr::filter(sim > sim.cutoff)

### Collect node data
result <- data.frame()

if (!is.null(object@enrichment_go_result)) {
  result <-
    object@enrichment_go_result@result %>%
    dplyr::filter(p.adjust < p.adjust.cutoff.go) %>%
    dplyr::filter(Count > count.cutoff.go) %>%
    dplyr::select(-ONTOLOGY) %>%
    dplyr::mutate(database = "GO") %>%
    rbind(result)
}

if (!is.null(object@enrichment_kegg_result)) {
  result <-
    object@enrichment_kegg_result@result %>%
    dplyr::filter(p.adjust < p.adjust.cutoff.kegg) %>%
    dplyr::filter(Count > count.cutoff.kegg) %>%
    dplyr::select(-c(category, subcategory)) %>%
    dplyr::mutate(database = "KEGG") %>%
    rbind(result)
}

if (!is.null(object@enrichment_reactome_result)) {
  result <-
    object@enrichment_reactome_result@result %>%
    dplyr::filter(p.adjust < p.adjust.cutoff.reactome) %>%
    dplyr::filter(Count > count.cutoff.reactome) %>%
    dplyr::mutate(database = "Reactome") %>%
    rbind(result)
}

node_data <-
  result %>%
  dplyr::rename(node = ID)

### Create tidygraph object
graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

## edge_betweenness clustering based on sim ====
subnetwork <-
  suppressWarnings(igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(igraph::edge_attr(graph_data, "sim"))))
#### Assign functional module label for pathways
cluster <-
  paste("Functional_module", as.character(igraph::membership(subnetwork)), sep = "_")

# Use cosine dissimilairty (1 - cosine) ====
mcl_res <- MCL::mcl(
  sim_matrix,
  addLoops = FALSE,
  inflation = 2.5,
  expansion = 2,
  max.iter = 500,
  allow1 = TRUE
)

binary_cut_res <- simplifyEnrichment::binary_cut(mat = sim_matrix, cutoff = 0.6)
simplifyEnrichment::plot_binary_cut(mat = sim_matrix, cutoff = 0.5)

#  Use embedding to perform clustering ====



