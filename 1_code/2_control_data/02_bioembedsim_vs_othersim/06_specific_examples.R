library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

library(reshape2)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/overlap_similarity/all_gene_list.rda")

load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_df.rda"
)

###GO similarity
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/GO_terms/wang_sim_df.rda"
)
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/GO_terms/resnik_sim_df.rda"
)
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/GO_terms/rel_sim_df.rda"
)
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/GO_terms/jiang_sim_df.rda"
)
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/GO_terms/lin_sim_df.rda"
)

###component based
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/overlap_similarity/op_sim_df.rda"
)
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/overlap_similarity/kappa_sim_df.rda"
)
load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/overlap_similarity/dice_sim_df.rda"
)

load(
  "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/overlap_similarity/jaccard_sim_df.rda"
)

dir.create("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/specific_examples")
setwd("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/specific_examples")

control_dt %>%
  dplyr::filter(id %in% c("GO:0098978", "GO:0098985", "hsa04724"))

control_dt %>%
  dplyr::filter(expected_module == "Functional_module_4")

embedding_data <-
  embedding_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )

op_data <-
  op_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )

kappa_data <-
  kappa_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )

dice_data <-
  dice_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )
wang_data <-
  wang_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )
resnik_data <-
  resnik_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )
rel_data <-
  rel_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )
jiang_data <-
  jiang_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )
lin_data <-
  lin_sim_df %>%
  dplyr::filter(
    from %in% c("GO:0098978", "GO:0098985", "hsa04724") &
      to %in% c("GO:0098978", "GO:0098985", "hsa04724")
  )

temp_data <-
  embedding_data %>%
  dplyr::left_join(jaccard_sim_df, by = c("from", "to")) %>%
  dplyr::rename(embedding = sim.x, jaccard = sim.y) %>%
  dplyr::left_join(op_data, by = c("from", "to")) %>%
  dplyr::rename(op = sim) %>%
  dplyr::left_join(kappa_data, by = c("from", "to")) %>%
  dplyr::rename(kappa = sim) %>%
  dplyr::left_join(dice_data, by = c("from", "to")) %>%
  dplyr::rename(dice = sim) %>%
  dplyr::left_join(wang_data, by = c("from", "to")) %>%
  dplyr::rename(wang = sim) %>%
  dplyr::left_join(resnik_data, by = c("from", "to")) %>%
  dplyr::rename(resnik = sim) %>%
  dplyr::left_join(rel_data, by = c("from", "to")) %>%
  dplyr::rename(rel = sim) %>%
  dplyr::left_join(jiang_data, by = c("from", "to")) %>%
  dplyr::rename(jiang = sim) %>%
  dplyr::left_join(lin_data, by = c("from", "to")) %>%
  dplyr::rename(lin = sim)

###venn diagram
library(VennDiagram)
venn_data <- list(
  "GO:0098978" = all_gene_list[["GO:0098978"]],
  "GO:0098985" = all_gene_list[["GO:0098985"]],
  "hsa04724" = all_gene_list[["hsa04724"]]
)

venn_plot <- venn.diagram(
  x = venn_data,
  filename = NULL,
  fill = c("#4f94d7", "#d15356", "#ebb869"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontface = "bold",
  # main = "Venn Diagram of Gene Overlaps",
  main.cex = 2,
  main.fontface = "bold",
  margin = 0.1
)

pdf(file = "venn_plot.pdf", width = 8, height = 8)
grid.draw(venn_plot)
dev.off()


####barplot
plot <-
temp_data %>%
  dplyr::select(from, to , embedding, jaccard, op, wang) %>%
  tidyr::pivot_longer(
    cols = -c(from, to),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(variable_id = paste(from, to, sep = "_")) %>%
  ggplot(aes(x = value, y = variable_id)) +
  geom_bar(aes(fill = metric),
           stat = "identity", position = "dodge",
           show.legend = FALSE,
           color = "black") +
  facet_wrap(~metric, scales = "fixed", nrow = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Similarity score", y = "") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  ggsci::scale_fill_bmj()

ggsave(plot = plot, filename = "similarity_barplot.pdf", width = 8, height = 6)
