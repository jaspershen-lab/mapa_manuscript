## Benchmarking on the informative GO set benchmark (k = 20)
## Step 6: Histogram of the ground-truth module size distribution
## (57 modules, 313 member GO terms)

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

ground_truth_dt <- readr::read_csv("2_data/informative_go_set/k_20/ground_truth.csv",
                                   show_col_types = FALSE)

module_size_dt <-
  ground_truth_dt |>
  dplyr::count(module_id, name = "module_size")

## Sanity checks against the dataset manifest (57 modules, 313 members)
stopifnot(nrow(module_size_dt) == 57,
          sum(module_size_dt$module_size) == 313)

readr::write_csv(module_size_dt,
                 "3_data_analysis/11_informative_go_benchmarking/module_size_dt.csv")

plot <-
  module_size_dt |>
  ggplot(aes(x = module_size)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey80") +
  scale_x_continuous(breaks = seq(4, max(module_size_dt$module_size), 2)) +
  theme_bw() +
  xlab("Module Size (Number of GO Terms)") +
  ylab("Number of Modules")

plot

ggsave(plot = plot,
       filename = "3_data_analysis/11_informative_go_benchmarking/module_size_distribution.pdf",
       width = 6, height = 4)
