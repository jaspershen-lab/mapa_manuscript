library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
library(org.Hs.eg.db)
source("1_code/100_tools.R")

setwd("2_data/case_study/02_monkey_multi_omics/")

# Thymus ====
load("taxid_9606/thymus/up/multi_omics_fm.rda")

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
                 module_id = "Functional_module_68",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/thymus/up/llm_interpreted_object.rda")

# Thymus (down) ====
load("taxid_9606/thymus/down/multi_omics_fm.rda")

multi_omics_fm$functional_module_result |> dplyr::filter(module_content_number > 3) |> nrow()
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
                 module_id = "Functional_module_7",
                 llm_text = TRUE)

save(llm_interpreted_object, file = "taxid_9606/thymus/down/llm_interpreted_object.rda")

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
                 module_id = "Functional_module_68",
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
                 module_id = "Functional_module_16",
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

save(llm_interpreted_object, file = "taxid_9606/ovary/up/llm_interpreted_object.rda").

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
