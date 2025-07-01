library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/4_random_module_dataset/final_result.Rdata")
load("3_data_analysis/4_random_module_dataset/formated_result.Rdata")

load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_df.rda"
)

names(final_result)

dir.create("3_data_analysis/4_random_module_dataset/intra_module_similarity")
setwd("3_data_analysis/4_random_module_dataset/intra_module_similarity")

temp_data <-
  1:nrow(formated_result) %>%
  purrr::map(function(i) {
    x <-
      formated_result$pathway_id[i]
    x <-
      stringr::str_split(x, ";") %>%
      unlist()
    sim <-
      embedding_sim_df %>%
      dplyr::filter(from %in% x & to %in% x) %>%
      pull(sim)
    data.frame(module = formated_result$module[i], sim = sim)
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

temp_data <-
  temp_data %>%
  dplyr::mutate(class = case_when(stringr::str_detect(module, "Random") ~ "Random", TRUE ~ "Real")) %>%
  dplyr::mutate(module = stringr::str_replace(module, "Random ", "")) %>%
  dplyr::mutate(module = stringr::str_replace(module, "Functional_module_", "")) %>%
  dplyr::mutate(module =  factor(module, levels = stringr::str_sort(unique(module), numeric = TRUE)))


temp_data %>%
  ggplot(aes(x = module, y = sim)) +
  geom_boxplot(
    aes(fill = class),
    alpha = 0.5,
    outlier.shape = NA,
    width = 0.5
  ) +
  geom_jitter(aes(fill = class),
              shape = 21,
              alpha = 1,
              width = 0.1) +
  scale_fill_manual(values = real_random_module_color) +
  labs(x = "Module", y = "Biotext embedding similarity") +
  theme_bw() +
  theme(legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1))


library(ggplot2)
library(ggpubr)

plot <-
  temp_data %>%
  ggplot(aes(x = module, y = sim, fill = class)) +
  geom_boxplot(
    alpha = 0.5,
    outlier.shape = NA,
    width = 0.5,
    position = position_dodge(width = 0.6)
  ) +
  geom_jitter(
    shape = 21,
    alpha = 1,
    size = 2,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6)
  ) +
  scale_fill_manual(values = real_random_module_color) +
  labs(x = "Module", y = "Biotext embedding similarity") +
  theme_bw() +
  theme(legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1)) +
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     hide.ns = TRUE)

plot
ggsave(plot = plot,
       filename = "intra_module_similarity.pdf",
       height = 6,
       width = 8)
