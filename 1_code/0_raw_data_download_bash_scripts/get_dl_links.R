library(tidyverse)

dled_files_name <- read_csv("2_data/metadata/downloaded_file_links/CRA015128_dled_file_names.txt", col_names = "filename")

all_links <- read_csv("2_data/metadata/links_CRA015128_multi-tissue_metformin.txt", col_names = FALSE)

all_files <- all_links %>%
  mutate(filename = basename(X1)) %>%
  rename(link = "X1")

to_dl_links <-
  dplyr::anti_join(all_files, dled_files_name, by = "filename")

to_dl_links <- to_dl_links$link

n <- length(to_dl_links)
num_parts <- 2
split_size <- ceiling(n / num_parts)

split_vector <- lapply(seq(1, n, by = split_size), function(i) {
  to_dl_links[i:min(i + split_size - 1, n)]
})

for (i in 1:length(split_vector)) {
  # m <- m + length(split_vector[[i]])
  writeLines(split_vector[[i]], paste0("2_data/metadata/to_dl_links/CRA015128_to_dl_links_", i, ".txt"))
}


