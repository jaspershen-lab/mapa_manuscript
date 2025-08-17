library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

dir.create(
  "3_data_analysis/1_pathway_redundancy/2_go",
  showWarnings = FALSE,
  recursive = TRUE
)

setwd("3_data_analysis/1_pathway_redundancy/2_go")

library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)


# # Get GO term annotations
# go_data <- as.data.frame(org.Hs.egGO)
#
# # Filter only "BP", "MF", or "CC"
# bp_data <-
#   go_data %>%
#   dplyr::filter(Ontology == "BP")
#
# mf_data <-
#   go_data %>%
#   dplyr::filter(Ontology == "MF")
#
# cc_data <-
#   go_data %>%
#   dplyr::filter(Ontology == "CC")
#
#
# # Helper to build list: GO ID -> list of Entrez gene IDs
# make_go_list <- function(df) {
#   split(df$gene_id, df$go_id)
# }
#
# go2gene_bp <- make_go_list(bp_data)
# go2gene_mf <- make_go_list(mf_data)
# go2gene_cc <- make_go_list(cc_data)
#
#
# # Optional: filter out GO terms with <2 genes
# go2gene_bp <- go2gene_bp[sapply(go2gene_bp, length) >= 10]
# go2gene_mf <- go2gene_mf[sapply(go2gene_mf, length) >= 10]
# go2gene_cc <- go2gene_cc[sapply(go2gene_cc, length) >= 10]


library(GOSemSim)

#####BP
# Use Biological Process (BP)
hsGO_BP <- godata('org.Hs.eg.db', ont = "BP", computeIC = TRUE)

# Option 1: All-vs-all matrix
go_terms <-
  names(go2gene_bp)

# go_terms_bp <- unique(hsGO_BP@geneAnno$GO)
## NOTE: "Wang" is IC-based methods. For IC-based methods, information of GO term is species specific.
sem_matrix_bp <- mgoSim(
  go_terms_bp,
  go_terms_bp,
  semData = hsGO_BP,
  measure = "Wang",
  combine= NULL
)

#
# save(sem_matrix_bp, file = "sem_matrix_bp.rda")

###MF

hsGO_MF <- godata('org.Hs.eg.db', ont = "MF", computeIC = TRUE)

go_terms <-
  names(go2gene_mf)

sem_matrix_mf <- mgoSim(
  go_terms,
  go_terms,
  semData = hsGO_MF,
  measure = "Wang",
  combine=NULL
)

save(sem_matrix_mf, file = "sem_matrix_mf.rda")

##CC
hsGO_CC <- godata('org.Hs.eg.db', ont = "CC", computeIC = TRUE)

go_terms <- names(go2gene_cc)

sem_matrix_cc <- mgoSim(
  go_terms,
  go_terms,
  semData = hsGO_CC,
  measure = "Wang",
  combine=NULL
)

save(sem_matrix_cc, file = "sem_matrix_cc.rda")



load("sem_matrix_bp.rda")
load("sem_matrix_mf.rda")
load("sem_matrix_cc.rda")

###BP
library(ggforce)
plot <-
  sem_matrix_bp[upper.tri(sem_matrix_bp)] %>%
  data.frame(value = .) %>%
  ggplot() +
  geom_histogram(aes(x = value),
                 color = "black",
                 fill = database_color["GO"],
                 bins = 50) +
  theme_bw() +
  labs(x = "Semantic Similarity", y = "Frequency") +
  ggforce::facet_zoom(
    xlim = c(0.4, 1),
    ylim = c(0, 15000),
    zoom.size = 0.5,
    show.area = TRUE
  )
plot
ggsave(plot,
       filename = "sem_matrix_bp_histogram.pdf",
       width = 10,
       height = 5)

##MF
plot <-
  sem_matrix_mf[upper.tri(sem_matrix_mf)] %>%
  data.frame(value = .) %>%
  ggplot() +
  geom_histogram(aes(x = value),
                 color = "black",
                 fill = "#1f77b4",
                 bins = 50) +
  theme_bw() +
  labs(x = "Semantic Similarity", y = "Frequency") +
  ggforce::facet_zoom(
    xlim = c(0.4, 1),
    ylim = c(0, 35000),
    zoom.size = 0.5,
    show.area = TRUE
  )

plot
ggsave(plot,
       filename = "sem_matrix_mf_histogram.pdf",
       width = 10,
       height = 5)

##CC
plot <-
  sem_matrix_cc[upper.tri(sem_matrix_cc)] %>%
  data.frame(value = .) %>%
  ggplot() +
  geom_histogram(aes(x = value),
                 color = "black",
                 fill = "#1f77b4",
                 bins = 50) +
  theme_bw() +
  labs(x = "Semantic Similarity", y = "Frequency") +
  ggforce::facet_zoom(
    xlim = c(0.4, 1),
    ylim = c(0, 25000),
    zoom.size = 0.5,
    show.area = TRUE
  )

plot

ggsave(plot,
       filename = "sem_matrix_cc_histogram.pdf",
       width = 10,
       height = 5)


####network
##BP
# Get upper triangle indices (to avoid duplicates and self-comparisons)
upper_idx <- which(upper.tri(sem_matrix_bp), arr.ind = TRUE)

# Construct data frame
sem_matrix_bp_df <- data.frame(
  GO_term1 = rownames(sem_matrix_bp)[upper_idx[, 1]],
  GO_term2 = colnames(sem_matrix_bp)[upper_idx[, 2]],
  similarity = sem_matrix_bp[upper_idx]
)

sem_matrix_bp_df <-
  sem_matrix_bp_df %>%
  dplyr::filter(similarity > 0.9)

dim(sem_matrix_bp_df)

edge_data <-
  sem_matrix_bp_df

node_data <-
  data.frame(GO_term = unique(c(edge_data$GO_term1, edge_data$GO_term2)), stringsAsFactors = FALSE)

# Create tidygraph object
library(tidygraph)
library(ggraph)

tidygraph_obj <- tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
  dplyr::mutate(degree = centrality_degree())


# Plot the graph
plot <-
  ggraph(tidygraph_obj, layout = "stress") +
  geom_edge_fan(aes(edge_width = similarity), show.legend = TRUE) +
  geom_node_point(
    aes(size = degree),
    fill = database_color["GO"],
    shape = 21,
    color = "black"
  ) +
  # geom_node_text(aes(label = GO_term), repel = TRUE, size = 3) +
  theme_void() +
  # labs(title = "GO Term Similarity Network (BP)", subtitle = "Edges represent high semantic similarity") +
  scale_edge_width(range = c(0.1, 2))
plot
ggsave(plot,
       filename = "sem_matrix_bp_network.pdf",
       width = 10,
       height = 5)

###detect modules
library(igraph)

tidygraph_obj_ig <- as.igraph(tidygraph_obj)
modules <- cluster_edge_betweenness(tidygraph_obj_ig, weights = 1 - E(tidygraph_obj_ig)$similarity)

modules

# Add module membership to node data
node_data$module <- membership(modules)

# Create a new tidygraph object with module information
tidygraph_obj <- tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
  dplyr::mutate(degree = centrality_degree(),
                module = node_data$module)
library(viridis)

plot2 <-
  ggraph(tidygraph_obj, layout = "stress") +
  geom_edge_fan(aes(edge_width = similarity), show.legend = TRUE) +
  geom_node_point(
    aes(size = degree,
        fill = as.factor(module)),
    shape = 21,
    color = "black"
  ) +
  geom_node_text(aes(label = GO_term), repel = TRUE, size = 3) +
  theme_void() +
  # labs(title = "GO Term Similarity Network (BP)", subtitle = "Edges represent high semantic similarity") +
  scale_edge_width(range = c(0.1, 2)) +
  scale_fill_manual(values = viridis::viridis(length(unique(node_data$module)))) +
  theme(legend.position = "bottom")
plot2
ggsave(plot2,
       filename = "sem_matrix_bp_network2.pdf",
       width = 10,
       height = 15)

plot3 <-
  ggraph(tidygraph_obj, layout = "stress") +
  geom_edge_fan(aes(edge_width = similarity),
                show.legend = TRUE) +
  geom_node_point(
    aes(size = degree,
        fill = as.factor(module)),
    shape = 21,
    color = "black",
    show.legend = FALSE
  ) +
  geom_node_text(aes(label = GO_term), repel = TRUE, size = 3, show.legend = FALSE) +
  theme_void() +
  # labs(title = "GO Term Similarity Network (BP)", subtitle = "Edges represent high semantic similarity") +
  scale_edge_width(range = c(0.1, 2))
  # scale_fill_manual(values = viridis::viridis(length(unique(node_data$module))))
plot3
ggsave(plot3,
       filename = "sem_matrix_bp_network3.pdf",
       width = 15,
       height = 5)

###only for moudle 1
module1_terms <- node_data$GO_term[node_data$module == 1]

tidygraph_obj_module1 <-
tidygraph_obj %>%
  activate(nodes) %>%
  dplyr::filter(module == 1)

plot_module1 <-
  ggraph(tidygraph_obj_module1, layout = "stress") +
  geom_edge_fan(aes(edge_width = similarity), show.legend = TRUE) +
  geom_node_point(
    aes(size = degree),
    fill = database_color["GO"],
    shape = 21,
    color = "black"
  ) +
  geom_node_text(aes(label = GO_term), repel = TRUE, size = 3) +
  theme_void() +
  labs(title = "GO Term Similarity Network (Module 1)",
       subtitle = paste("Contains", length(module1_terms), "terms")) +
  scale_edge_width(range = c(0.1, 2))

plot_module1
ggsave(plot_module1,
       filename = "sem_matrix_bp_network_module1.pdf",
       width = 5,
       height = 5)


## Extract module info
library(GO.db)

# GOTERM[[node_data$GO_term[1]]]@Term
node_data_with_term_name <-
  node_data |> dplyr::mutate(term_name = Term(GO_term))

node_data_with_term_name <- node_data_with_term_name |> dplyr::arrange(module)
node_data_with_term_name <- node_data_with_term_name |> dplyr::select(module, GO_term, term_name)

save(node_data_with_term_name, file = "node_data_with_term_name.rda")
rio::export(node_data_with_term_name, file = "go_module_info.xlsx")
