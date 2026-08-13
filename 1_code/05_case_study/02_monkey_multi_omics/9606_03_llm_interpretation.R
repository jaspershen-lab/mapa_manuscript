library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
library(org.Hs.eg.db)
source("1_code/100_tools.R")

setwd("2_data/case_study/02_monkey_multi_omics/")

# Thymus ====
load("taxid_9606/thymus/up/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_46",
                 show_rwr_edge = TRUE,
                 llm_text = TRUE)

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_15",
                 show_rwr_edge = TRUE,
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/thymus/up/llm_interpreted_object.rda")

# Thymus (down) ====
load("taxid_9606/thymus/down/multi_omics_fm.rda")

# multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
#   dplyr::filter(multi_omics_num == 3 & module %in% c("Functional_module_5", "Functional_module_2", "Functional_module_66"))

# multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
#   dplyr::filter(multi_omics_num == 3 & module == "Functional_module_5")
#
# plot_module_info(multi_omics_fm,
#                  module_id = "Functional_module_66",
#                  show_rwr_edge = TRUE)

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |> filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               phenotype = "aging",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 15,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = TRUE,
                                               api_provider = "openai")

save(llm_interpreted_object, file = "taxid_9606/thymus/down/llm_interpreted_object.rda")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_5",
                 show_rwr_edge = TRUE,
                 llm_text = TRUE)

# Spleen ====
load("taxid_9606/spleen/up/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               phenotype = "aging",
                                               local_corpus_dir = NULL,
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 15,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_65",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/spleen/up/llm_interpreted_object.rda")

# Spleen (down) ====
load("taxid_9606/spleen/down/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_65",
                 show_rwr_edge = TRUE,
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/spleen/down/llm_interpreted_object.rda")

# Ovary ====
load("taxid_9606/ovary/up/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_10",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/ovary/up/llm_interpreted_object.rda")

# Ovary (down) ====
load("taxid_9606/ovary/down/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_10",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/ovary/down/llm_interpreted_object.rda")

# Stomach ====
load("taxid_9606/stomach/up/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_1",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/stomach/up/llm_interpreted_object.rda")

# Stomach (down) ====
load("taxid_9606/stomach/down/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(
    module == "Functional_module_6",
    module_content_number >= 3,
    true_omics_num == 3,
    include_metabolites
  )

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_6",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/stomach/down/llm_interpreted_object.rda")

# Kidney ====
load("taxid_9606/kidney/up/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_1",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/kidney/up/llm_interpreted_object.rda")

# Kidney (down) ====
load("taxid_9606/kidney/down/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_1",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/kidney/down/llm_interpreted_object.rda")

# Aortic arch ====
load("taxid_9606/aortic_arch/up/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 15,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_67",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/aortic_arch/up/llm_interpreted_object.rda")

# Aortic arch (down) ====
load("taxid_9606/aortic_arch/down/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_67",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/aortic_arch/down/llm_interpreted_object.rda")

# Thyroid ====
load("taxid_9606/thyroid/up/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(multi_omics_num == 3)

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 3,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

set.seed(1)

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_24",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/thyroid/up/llm_interpreted_object.rda")

# Thyroid (down) ====
load("taxid_9606/thyroid/down/multi_omics_fm.rda")

multi_omics_fm$functional_module_result <- multi_omics_fm$functional_module_result |>
  dplyr::filter(
    module == "Functional_module_18",
    module_content_number >= 3,
    true_omics_num == 3,
    include_metabolites
  )

llm_interpreted_object <- llm_interpret_module(object = multi_omics_fm,
                                               module_content_number_cutoff = 2,
                                               llm_model = "gpt-4o-mini-2024-07-18",
                                               embedding_model = "text-embedding-3-small",
                                               api_key = api_key,
                                               embedding_output_dir = "embedding_output/",
                                               local_corpus_dir = NULL,
                                               phenotype = "aging",
                                               chunk_size = 5,
                                               years = 5,
                                               retmax = 20,
                                               similarity_filter_num = 20,
                                               GPT_filter_num = 5,
                                               orgdb = org.Hs.eg.db,
                                               output_prompt = FALSE,
                                               api_provider = "openai")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_18",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/thyroid/down/llm_interpreted_object.rda")

# result summary ====
setwd(r4projects::get_project_wd())
rm(list = ls())
setwd("2_data/case_study/02_monkey_multi_omics/taxid_9606/")

# tissues <- c("thymus", "spleen", "ovary",
#              # "stomach", "kidney", "thyroid",
#              "aortic_arch")
# directions <- c("up", "down")
#
# llm_summary <-
#   purrr::map_dfr(tissues, function(tissue) {
#     purrr::map_dfr(directions, function(direction) {
#       rda_path <- file.path(tissue, direction, "llm_interpreted_object.rda")
#       if (!file.exists(rda_path)) return(NULL)
#
#       load(rda_path)  # loads llm_interpreted_object
#
#       interp <- llm_interpreted_object$llm_module_interpretation
#
#       purrr::map_dfr(names(interp), function(module_id) {
#         gn <- interp[[module_id]]$generated_name
#         tibble::tibble(
#           tissue = tissue,
#           direction = direction,
#           module_id = module_id,
#           module_name = if (!is.null(gn$module_name)) gn$module_name else NA_character_,
#           module_summary  = if (!is.null(gn$summary)) gn$summary else NA_character_
#         )
#       })
#     })
#   })
#
# llm_summary$tissue |> unique()
#
# rio::export(llm_summary, file = "llm_summary.xlsx")

## Thymus down ====
load("thymus/down/llm_interpreted_object.rda")

# thymus down fm_2
set.seed(1)

plot <-
plot_module_info(object = llm_interpreted_object,
                 module_id = "Functional_module_2",
                 show_rwr_edge = TRUE,
                 llm_text = TRUE,
                 node_size = 6
                 )
plot

library(Cairo)
CairoPDF("thymus/down/thymus_down_fm_2.pdf", width = 8, height = 7)
plot
dev.off()

# thymus down fm_5
set.seed(1)

plot <-
  plot_module_info(object = llm_interpreted_object,
                   module_id = "Functional_module_5",
                   show_rwr_edge = TRUE,
                   node_size = 6,
                   label_size = 3.5
  )
plot

library(Cairo)
CairoPDF("thymus/down/thymus_down_fm_5.pdf", width = 8, height = 7)
plot
dev.off()

# thymus down fm_66
set.seed(39)

plot <-
  plot_module_info(object = llm_interpreted_object,
                   module_id = "Functional_module_66",
                   show_rwr_edge = TRUE,
                   node_size = 6,
                   label_size = 3.5
  )
plot

library(Cairo)
CairoPDF("thymus/down/thymus_down_fm_66.pdf", width = 8, height = 7)
plot
dev.off()

## Spleen up ====
load("spleen/up/llm_interpreted_object.rda")

# spleen up fm_65
set.seed(13)

plot <-
  plot_module_info(object = llm_interpreted_object,
                   module_id = "Functional_module_65",
                   show_rwr_edge = TRUE,
                   llm_text = TRUE,
                   node_size = 6,
                   label_size = 3.5
  )
plot

library(Cairo)
CairoPDF("spleen/up/spleen_up_fm_65.pdf", width = 8, height = 7)
plot
dev.off()

## Aortic arch up ====
load("aortic_arch/up/llm_interpreted_object.rda")

# aortic arch up fm_67
set.seed(12)

plot <-
  plot_module_info(object = llm_interpreted_object,
                   module_id = "Functional_module_67",
                   show_rwr_edge = TRUE,
                   llm_text = TRUE,
                   node_size = 6,
                   label_size = 3.5
  )
plot

library(Cairo)
CairoPDF("aortic_arch/up/aortic_arch_up_fm_67.pdf", width = 8, height = 7)
plot
dev.off()

# thyroid up fm_24
load("thyroid/up/llm_interpreted_object.rda")

set.seed(7)

plot <-
  plot_module_info(object = llm_interpreted_object,
                   module_id = "Functional_module_24",
                   show_rwr_edge = TRUE,
                   llm_text = TRUE,
                   node_size = 6,
                   label_size = 3.5
  )
plot

library(Cairo)
CairoPDF("thyroid/up/thyroid_up_fm_24.pdf", width = 8, height = 7)
plot
dev.off()

###
llm_interpreted_object$functional_module_result |>
  filter(module == "Functional_module_67") |>
  pull(pathways) |>
  stringr::str_replace_all("/", ",")

llm_interpreted_object$functional_module_result |>
  filter(module == "Functional_module_67") |>
  pull(genes) |>
  stringr::str_replace_all("/", ",")

plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_16",
                 show_rwr_edge = TRUE,
                 llm_text = TRUE)
plot_module_info(llm_interpreted_object,
                 module_id = "Functional_module_50",
                 llm_text = TRUE)

llm_interpreted_object$llm_module_interpretation$Functional_module_53$generated_name

test_beta <- test[, c(1, 8:37)]
test_p <- test[, c(1,38:67)]
