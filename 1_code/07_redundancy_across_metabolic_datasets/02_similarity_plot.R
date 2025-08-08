library(r4projects)
setwd(get_project_wd())
rm(list = ls())

# load the similarity table
load("2_data/smpdb_kegg_similarity.rda")
load("2_data/reactome_kegg_similarity.rda")
load("2_data/smpdb_reactome_similarity.rda")

library(ggplot2)

# reactome and kegg
plot <-
  reactome_kegg_result %>%
  dplyr::filter(!is.na(similarity)) %>%
  ggplot(aes(x = similarity)) +
  geom_histogram(
    color = "black",
    fill = "#FF8C00",
    bins = 50
  ) +
  theme_bw() +
  labs(x = "Jaccard Similarity", y = "Frequency") +
  ggforce::facet_zoom(
    xlim = c(0.3, 1),
    ylim = c(0, 1250),
    zoom.size = 0.5,
    show.area = TRUE
  )

print(plot)

ggsave(plot,
       filename = "3_data_analysis/06_metabolic_pathway_redundancy/reactome_kegg.pdf",
       width = 10,
       height = 5)

# smpdb and kegg
plot <-
  smpdb_kegg_result %>%
  dplyr::filter(!is.na(similarity)) %>%
  ggplot(aes(x = similarity)) +
  geom_histogram(
    color = "black",
    fill = "#FF8C00",
    bins = 50
  ) +
  theme_bw() +
  labs(x = "Jaccard Similarity", y = "Frequency") +
  ggforce::facet_zoom(
    xlim = c(0.3, 1),
    ylim = c(0, 400),
    zoom.size = 0.5,
    show.area = TRUE
  )

print(plot)

ggsave(plot,
       filename = "3_data_analysis/06_metabolic_pathway_redundancy/smpdb_kegg.pdf",
       width = 10,
       height = 5)

# smpdb and reactome
plot <-
  smpdb_reactome_result %>%
  dplyr::filter(!is.na(similarity)) %>%
  ggplot(aes(x = similarity)) +
  geom_histogram(
    color = "black",
    fill = "#FF8C00",
    bins = 50
  ) +
  theme_bw() +
  labs(x = "Jaccard Similarity", y = "Frequency") +
  ggforce::facet_zoom(
    xlim = c(0.3, 1),
    ylim = c(0, 2000),
    zoom.size = 0.5,
    show.area = TRUE
  )

print(plot)

ggsave(plot,
       filename = "3_data_analysis/06_metabolic_pathway_redundancy/smpdb_reactome.pdf",
       width = 10,
       height = 5)
