library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(tidyverse)
library(ggvenn)
library(glue)

# load three metabolic pathway datasets
load("2_data/hmdb_pathway.rda")
load("2_data/kegg_hsa_pathway.rda")
load("2_data/reactome_human_pathway.rda")

smpdb_database <- hmdb_pathway
kegg_database <- kegg_hsa_pathway
reactome_database <- reactome_human_pathway_final

# remove pathways with no HMDB IDs
smpdb_database_filtered <- smpdb_database[!is.na(smpdb_database$HMDB_ID), ]
kegg_database_filtered <- kegg_database[!is.na(kegg_database$HMDB_ID), ]
reactome_database_filtered <- reactome_database[!is.na(reactome_database$HMDB_ID), ]

# function to split hmdb id
split_hmdb <- function(hmdb_str) {
  hmdb_str %>%
    str_split("\\{\\}") %>%
    unlist() %>%
    unique()
}

# get the list of HMDB IDs for each pathway
smpdb_list <- smpdb_database_filtered %>%
  mutate(HMDB_vec = map(HMDB_ID, split_hmdb)) %>%
  transmute(
    pathway_id_smpdb = pathway_id,
    pathway_name = pathway_name,
    HMDB_vec_smpdb  = HMDB_vec
  )


kegg_list <- kegg_database_filtered %>%
  mutate(HMDB_vec = map(HMDB_ID, split_hmdb)) %>%
  transmute(
    pathway_id_kegg = pathway_id,
    pathway_name = pathway_name,
    HMDB_vec_kegg  = HMDB_vec
  )

reactome_list <- reactome_database_filtered %>%
  mutate(HMDB_vec = map(HMDB_ID, split_hmdb)) %>%
  transmute(
    pathway_id_reactome = pathway_id,
    pathway_name = pathway_name,
    HMDB_vec_reactome  = HMDB_vec
  )

# get glycolysis from smpdb and kegg
smpdb_glycolysis <- smpdb_list %>%
  filter(pathway_id_smpdb == "SMP0000040")

kegg_glycolysis <- kegg_list %>%
  filter(pathway_id_kegg == "hsa00010")

## draw venn diagram for glycolysis
# get the unique HMDB IDs for SMPDB and KEGG glycolysis pathways
s_set <- smpdb_glycolysis$HMDB_vec_smpdb %>% unlist() %>% unique()
k_set <- kegg_glycolysis$HMDB_vec_kegg %>% unlist() %>% unique()

# draw venn plot
p_venn <- ggvenn(
  data = list(SMPDB = s_set, KEGG = k_set),
  fill_color = c("#7998ad", "#fd7541"),
  stroke_size = 0.6,
  set_name_size = 4,
  text_size = 4,
  show_percentage = FALSE
) +
  ggtitle("glycolysis in SMPDB and KEGG") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_venn)

ggsave(
  plot = p_venn,
  filename = "3_data_analysis/06_metabolic_pathway_redundancy/glycolysis.pdf",
  width = 6,
  height = 6
)

# get fatty acid metabolism from smpdb and reactome
smpdb_fatty<- smpdb_list %>%
  filter(pathway_id_smpdb == "SMP0000051")

reactome_fatty <- reactome_list %>%
  filter(pathway_id_reactome == "R-HSA-8978868")

## draw venn diagram for fatty acid metabolism
# get the unique HMDB IDs for SMPDB and Reactome fatty acid metabolism pathways
s_set <- smpdb_fatty$HMDB_vec_smpdb %>% unlist() %>% unique()
r_set <- reactome_fatty$HMDB_vec_reactome %>% unlist() %>% unique()

# draw venn plot
p_venn <- ggvenn(
  data = list(SMPDB = s_set, Reactome = r_set),
  fill_color = c("#7998ad", "#23b9c7"),
  stroke_size = 0.6,
  set_name_size = 4,
  text_size = 4,
  show_percentage = FALSE
) +
  ggtitle("fatty acid metabolism in SMPDB and Reactome") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_venn)

ggsave(
  plot = p_venn,
  filename = "3_data_analysis/06_metabolic_pathway_redundancy/fatty_acid_metabolism.pdf",
  width = 6,
  height = 6
)




