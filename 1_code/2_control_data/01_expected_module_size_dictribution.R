library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

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
ggsave(plot = plot,
       filename = "3_data_analysis/02_control_data/01_expected_module_size_dictribution/expected_module_size_dictribution.pdf",
       height = 8,
       width = 8
       )
