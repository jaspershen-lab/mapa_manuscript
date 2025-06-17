library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

setwd("3_data_analysis/02_control_data/01_expected_module_size_dictribution/")

plot <-
  control_dt |>
    dplyr::pull(expected_module) |>
    as.factor() |>
    table() |>
    as.data.frame() |>
    ggplot(aes(x = reorder(Var1, Freq), y = Freq)) +
    geom_col() +
    theme_bw() +
    theme(
      # axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 1)
      axis.text.x = element_text(angle = 60, hjust = 1)
    ) +
    xlab("Expected Functional Module") +
    ylab("Size of Expected Functional Module")

plot


plot <-
  control_dt %>%
  dplyr::mutate(expected_module = stringr::str_replace_all(expected_module, "Functional_module_", "")) %>%
  dplyr::mutate(
    expected_module = factor(expected_module, levels = stringr::str_sort(unique(expected_module), numeric = TRUE)),
    database = factor(database, levels = c("GO", "KEGG", "Reactome"))
  ) %>%
  ggplot() +
  geom_bar(aes(expected_module, fill = database),
           color = "black") +
  theme_bw() +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1)
  ) +
  xlab("") +
  ylab("Size of expected module") +
  scale_fill_manual(values = database_color)

plot

ggsave(plot = plot,
       filename = "expected_module_size_dictribution.pdf",
       height = 5,
       width = 8
       )
