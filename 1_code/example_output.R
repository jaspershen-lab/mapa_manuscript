library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# dir.create("3_data_analysis/07_example_output")

## MAPA workflow
setwd("3_data_analysis/07_example_output/")

example_dt <- read_csv("example_ora_gene_data.csv")

variable_info <- convert_id(
  data = example_dt,
  query_type = "gene",
  from_id_type = "ensembl",
  organism = "org.Hs.eg.db"
)

enriched_pathways_res <- enrich_pathway(
  variable_info = variable_info,
  query_type = "gene",
  database = c("go", "kegg", "reactome"),
  go.orgdb = "org.Hs.eg.db",
  go.keytype = "ENTREZID",
  kegg.organism = "hsa",
  kegg.keytype = "kegg",
  reactome.organism = "human"
)

enriched_pathways_res

openai_semantic_sim_matrix <-
  get_bioembedsim(object = enriched_pathways_res,
                  api_provider = "openai",
                  text_embedding_model = "text-embedding-3-small",
                  api_key = api_key,
                  database = c("go", "kegg", "reactome"),
                  save_to_local = FALSE)

functional_modules_res <- get_functional_modules(
  object = openai_semantic_sim_matrix,
  sim.cutoff = 0.7,
  cluster_method = "louvain"
)

save(openai_semantic_sim_matrix, file = "openai_semantic_sim_matrix.rda")
save(functional_modules_res, file = "functional_modules_res.rda")

asse_report <- assess_clustering_quality(
  object = functional_modules_res
)

asse_report$evaluation_plot
asse_report$size_plot

dir.create("embedding_output")

llm_interpreted_res <-
  llm_interpret_module(
    object = functional_modules_res,
    module_content_number_cutoff = 2,
    api_key = api_key,
    embedding_output_dir = "embedding_output"
  )

save(llm_interpreted_res, file = "llm_interpreted_res.rda")

## pathway enrichment barplot ====
plot <-
  plot_pathway_bar(
    object = llm_interpreted_res,
    top_n = 10,
    y_label_width = 30,
    level = "functional_module",
    line_type = "straight",
    database = c("go", "kegg", "reactome"),
    llm_text = TRUE
  )
plot

ggsave(plot = plot, filename = "enrich_barplot_fm.pdf",
       width = 8, height = 6)

## module info plot ====
plot <-
  plot_module_info(
    object = llm_interpreted_res,
    level = "functional_module",
    module_id = "Functional_module_8",
    llm_text = TRUE
  )

plot[[1]]
plot[[2]]
plot[[3]]

library(Cairo)
CairoPDF("fm_8_network.pdf", width = 8, height = 6)
plot[[1]]
dev.off()

CairoPDF("fm_8_enriched_pathways.pdf", width = 8, height = 6)
plot[[2]]
dev.off()

CairoPDF("fm_8_wordcloud.pdf", width = 8, height = 6)
plot[[3]]
dev.off()

## similarity network ====
plot <-
  plot_similarity_network(
    object = llm_interpreted_res,
    level = "functional_module",
    database = c("go", "kegg", "reactome"),
    degree_cutoff = 3,
    llm_text = TRUE
  )

plot

library(Cairo)
CairoPDF("similarity_network.pdf", width = 11, height = 8)
plot
dev.off()

## multipartite plot with expression level heatmap ====
load("input_data.rda")
expression_data <- input_data |>
  dplyr::select(-variable_id, -score, -fdr, -class,
                -symbol, -uniprot, -entrezid) |>
  dplyr::rename(id = ensembl)

selected_module_ids <- llm_interpreted_res@merged_module$functional_module_result |> filter(module_content_number == 4 |module_content_number == 5) |> pull(module)
module_names <- llm_interpreted_res@merged_module$functional_module_result |> filter(module %in% selected_module_ids) |> dplyr::select(module, llm_module_name)

plot <-
  plot_relationship_heatmap(
    object = llm_interpreted_res,
    level = "pathway",
    expression_data = expression_data,
    module_ids = selected_module_ids,
    llm_text = TRUE
  )

plot

library(Cairo)
CairoPDF("relationship_expression_level_hp.pdf", width = 10, height = 8)
plot
dev.off()

## multipartite ====
load("llm_interpreted_res.rda")
selected_module_ids <- llm_interpreted_res@merged_module$functional_module_result |> filter(module_content_number == 4 |module_content_number == 5) |> pull(module)
llm_interpreted_res@merged_module$functional_module_result <-
  llm_interpreted_res@merged_module$functional_module_result |>
  dplyr::filter(module %in% selected_module_ids)

plot <-
  plot_relationship_network(
    object = llm_interpreted_res,
    include_functional_modules = TRUE,
    include_modules = FALSE,
    include_pathways = TRUE,
    include_molecules = TRUE,
    include_variables = FALSE,
    llm_text = TRUE
  )

plot

library(Cairo)
CairoPDF("relationship_multipartite.pdf", width = 8, height = 8)
plot
dev.off()

## export result report ====
report_functional_module(
  object = llm_interpreted_res,
  path = ".",
  type = "md"
)
