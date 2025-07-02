library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

dir.create(
  "3_data_analysis/1_pathway_redundancy/3_metabolite_across_db",
  showWarnings = FALSE,
  recursive = TRUE
)

load("3_data_analysis/03_curated_pathway_dataset2/hmdb_kegg_matched_with_sim.rda")
load("3_data_analysis/03_curated_pathway_dataset2/hmdb_reactome_matched_with_sim.rda")
load("3_data_analysis/03_curated_pathway_dataset2/kegg_reactome_matched_with_sim.rda")

setwd("3_data_analysis/1_pathway_redundancy/3_metabolite_across_db")

hmdb_kegg_id <- hmdb_kegg_matched_with_sim |> select(hmdb_pathway_id, kegg_pathway_id)
hmdb_reactome_id <- hmdb_reactome_matched_with_sim |> select(hmdb_pathway_id, reactome_pathway_id)
kegg_reactome_id <- kegg_reactome_matched_with_sim |> select(kegg_pathway_id, reactome_pathway_id)

all_overlap_id <- left_join(hmdb_reactome_id, hmdb_kegg_id, by = "hmdb_pathway_id")
all_overlap_id <- all_overlap_id |> filter(!is.na(kegg_pathway_id))

graph_data <- data.frame()

for (i in 1:nrow(all_overlap_id)) {
  row_ids <- all_overlap_id[i,]
  hmdb_kegg_sim <- hmdb_kegg_matched_with_sim |>
    filter(hmdb_pathway_id == all_overlap_id[i,1] & kegg_pathway_id == all_overlap_id[i,3]) |>
    rename(pathway_name_1 = hmdb_pathway_name,
           pathway_name_2 = kegg_pathway_name,
           id_1 = hmdb_pathway_id,
           id_2 = kegg_pathway_id)
  hmdb_reactome_sim <- hmdb_reactome_matched_with_sim |>
    filter(hmdb_pathway_id == all_overlap_id[i,1] & reactome_pathway_id == all_overlap_id[i,2]) |>
    rename(pathway_name_1 = hmdb_pathway_name,
           pathway_name_2 = reactome_pathway_name,
           id_1 = hmdb_pathway_id,
           id_2 = reactome_pathway_id)
  graph_data <- rbind(graph_data, hmdb_kegg_sim, hmdb_reactome_sim)
}

hmdb_kegg_matched_sim_info <- hmdb_kegg_matched_with_sim |>
  rename(id_1 = hmdb_pathway_id,
         id_2 = kegg_pathway_id,
         pathway_name_1 = hmdb_pathway_name,
         pathway_name_2 = kegg_pathway_name)

hmdb_reactome_matched_sim_info <- hmdb_reactome_matched_with_sim |>
  rename(id_1 = hmdb_pathway_id,
         id_2 = reactome_pathway_id,
         pathway_name_1 = hmdb_pathway_name,
         pathway_name_2 = reactome_pathway_name)

kegg_reactome_matched_sim_info <- kegg_reactome_matched_with_sim |>
  rename(id_1 = kegg_pathway_id,
         id_2 = reactome_pathway_id,
         pathway_name_1 = kegg_pathway_name,
         pathway_name_2 = reactome_pathway_name)

graph_data <- graph_data |>
  bind_rows(hmdb_kegg_matched_sim_info, hmdb_reactome_matched_sim_info, kegg_reactome_matched_sim_info)

graph_data <- distinct_all(graph_data, .keep_all = TRUE)

edge_data <- graph_data |> select(id_1, id_2, jw_dist)
node_data <- data.frame(pathway_id = unique(c(graph_data$id_1, graph_data$id_2)),
                        stringsAsFactors = FALSE) |>
  mutate(database = case_when(
    grepl("^SMP", pathway_id) ~ "SMPDB",
    grepl("^hsa", pathway_id) ~ "KEGG",
    grepl("^R-HSA", pathway_id) ~ "Reactome"
  ))

smpdb_pathway_name_id <- bind_rows(
  data.frame(pathway_name = hmdb_kegg_matched_sim_info$pathway_name_1, pathway_id = hmdb_kegg_matched_sim_info$id_1),
  data.frame(pathway_name = hmdb_reactome_matched_sim_info$pathway_name_1, pathway_id = hmdb_reactome_matched_sim_info$id_1)
) |>
  distinct_all(.keep_all = TRUE)

kegg_pathway_name_id <- bind_rows(
  data.frame(pathway_name = hmdb_kegg_matched_sim_info$pathway_name_2, pathway_id = hmdb_kegg_matched_sim_info$id_2),
  data.frame(pathway_name = kegg_reactome_matched_sim_info$pathway_name_1, pathway_id = kegg_reactome_matched_sim_info$id_1)
) |>
  distinct_all(.keep_all = TRUE)

reactome_pathway_name_id <- bind_rows(
  data.frame(pathway_name = hmdb_reactome_matched_sim_info$pathway_name_2, pathway_id = hmdb_reactome_matched_sim_info$id_2),
  data.frame(pathway_name = kegg_reactome_matched_sim_info$pathway_name_2, pathway_id = kegg_reactome_matched_sim_info$id_2)
) |>
  distinct_all(.keep_all = TRUE)

name_lookup <- bind_rows(smpdb_pathway_name_id, kegg_pathway_name_id, reactome_pathway_name_id)

node_data <- node_data |>
  left_join(name_lookup, by = "pathway_id")

# Create tidygraph object
library(tidygraph)
library(ggraph)

tidygraph_obj <- tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
  dplyr::mutate(degree = centrality_degree())


# Plot the graph

plot <-
  # ggraph(tidygraph_obj, layout = "fr") +
  ggraph(tidygraph_obj, layout = "stress") +
  geom_edge_link(aes(edge_width = 1-jw_dist), show.legend = TRUE) +
  geom_node_point(
    aes(color = database,
        size = degree),
    shape = 20
  ) +
  geom_node_text(aes(label = pathway_name), repel = TRUE, size = 3) +
  theme_void() +
  scale_color_manual(values = database_color) +
  # labs(title = "GO Term Similarity Network (BP)", subtitle = "Edges represent high semantic similarity") +
  scale_edge_width(range = c(0.5, 1)) +
  labs(edge_width = "1 - Jaro-Winkler Distance")


plot

ggsave(plot,
       filename = "metabolite_across_db_redundancy.pdf",
       width = 10,
       height = 5)

# only show "Pentose phosphate pathway" related label
node_data_2 <- node_data |>
  mutate(pathway_name = if_else(grepl("Pentose phosphate pathway", pathway_name, ignore.case = TRUE), pathway_name, NA))

tidygraph_obj_2 <- tbl_graph(nodes = node_data_2,
                           edges = edge_data,
                           directed = FALSE) %>%
  dplyr::mutate(degree = centrality_degree())


# Plot the graph
plot2 <-
  # ggraph(tidygraph_obj, layout = "fr") +
  ggraph(tidygraph_obj_2, layout = "stress") +
  geom_edge_link(aes(edge_width = 1-jw_dist), show.legend = TRUE) +
  geom_node_point(
    aes(color = database,
        size = degree),
    shape = 20
  ) +
  geom_node_text(aes(label = pathway_name), repel = TRUE, size = 3) +
  theme_void() +
  scale_color_manual(values = database_color) +
  # labs(title = "GO Term Similarity Network (BP)", subtitle = "Edges represent high semantic similarity") +
  scale_edge_width(range = c(0.5, 1)) +
  labs(edge_width = "1 - Jaro-Winkler Distance")


plot2

ggsave(plot2,
       filename = "metabolite_across_db_redundancy_label_pentose_phosphate.pdf",
       width = 10,
       height = 5)


