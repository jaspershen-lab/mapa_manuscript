library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

setwd("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/")

## All pathways ====
load("combined_sim_df.rda")

temp_data <-
  combined_sim_df %>%
  dplyr::left_join(control_dt[, c("id", "database", "name", "expected_module")], by = c("from" = "id")) %>%
  dplyr::left_join(control_dt[, c("id", "database", "name", "expected_module")], by = c("to" = "id")) %>%
  dplyr::mutate(
    class = case_when(
      expected_module.x == expected_module.y ~ "same_module",
      expected_module.x != expected_module.y ~ "different_module"
    )
  ) %>%
  dplyr::mutate(
    class2 = case_when(
      database.x == database.y ~ "intra_database",
      database.x != database.y ~ "inter_database"
    )
  ) %>%
  dplyr::mutate(class = factor(class, levels = c("same_module", "different_module")))

# Calculate Cliff's delta for each metric
library(effectsize)

similarity_data <- temp_data |>
  dplyr::select(from, to, sim_embedding, sim_jaccard, sim_op, sim_kappa, sim_dice, class) %>%
  tidyr::pivot_longer(
    cols = c(sim_embedding, sim_jaccard, sim_op, sim_kappa, sim_dice),
    names_to = "metric",
    values_to = "similarity",
    names_prefix = "sim_"
  )

metrics <- unique(similarity_data$metric)

cliffs_results <- data.frame(
  metric = character(),
  cliffs_delta = numeric(),
  stringsAsFactors = FALSE
)

for (metric_name in metrics) {
  # Get data for current metric
  metric_data <- similarity_data[similarity_data$metric == metric_name, ]
  intra_vals <- metric_data$similarity[metric_data$class == "same_module"]
  inter_vals <- metric_data$similarity[metric_data$class == "different_module"]

  # Calculate Cliff's delta using effsize package
  cliff_result <- cliff.delta(intra_vals, inter_vals)
  delta_value <- cliff_result$estimate

  cliffs_results <- rbind(cliffs_results, data.frame(
    metric = metric_name,
    cliffs_delta = delta_value
  ))
}

save(cliffs_results, file = "all_pathways_cliffs_results.rda")

## calculate IQR
iqr_results <- data.frame(
  metric = character(),
  intra_q1 = numeric(),
  intra_q3 = numeric(),
  inter_q1 = numeric(),
  inter_q3 = numeric(),
  stringsAsFactors = FALSE
)

for (metric_name in metrics) {
  # Get data for current metric
  metric_data <- similarity_data[similarity_data$metric == metric_name, ]
  intra_vals <- metric_data$similarity[metric_data$class == "same_module"]
  intra_q1 <- quantile(intra_vals, 0.25)
  intra_q3 <- quantile(intra_vals, 0.75)

  inter_vals <- metric_data$similarity[metric_data$class == "different_module"]
  inter_q1 <- quantile(inter_vals, 0.25)
  inter_q3 <- quantile(inter_vals, 0.75)

  iqr_results <- rbind(iqr_results, data.frame(
    metric = metric_name,
    intra_q1 = intra_q1,
    intra_q3 = intra_q3,
    inter_q1 = inter_q1,
    inter_q3 = inter_q3
  ))
}

export(iqr_results, file = "iqr_results.xlsx")

## GO terms ====
load("GO_terms/go_combined_sim_df.rda")

temp_data <-
  go_combined_sim_df %>%
  dplyr::left_join(control_dt[, c("id", "database", "name", "expected_module")], by = c("from" = "id")) %>%
  dplyr::left_join(control_dt[, c("id", "database", "name", "expected_module")], by = c("to" = "id")) %>%
  dplyr::mutate(
    class = case_when(
      expected_module.x == expected_module.y ~ "same_module",
      expected_module.x != expected_module.y ~ "different_module"
    )
  ) %>%
  dplyr::mutate(
    class2 = case_when(
      database.x == database.y ~ "intra_database",
      database.x != database.y ~ "inter_database"
    )
  ) %>%
  dplyr::mutate(class = factor(class, levels = c("same_module", "different_module")))


# Calculate Cliff's delta for each metric
library(effectsize)

similarity_data <- temp_data |>
  dplyr::select(from, to, sim_embedding, sim_jaccard, sim_op, sim_kappa, sim_dice,
                sim_wang, sim_resnik, sim_rel, sim_jiang, sim_lin, class) %>%
  tidyr::pivot_longer(
    cols = c(sim_embedding, sim_jaccard, sim_op, sim_kappa, sim_dice,
             sim_wang, sim_resnik, sim_rel, sim_jiang, sim_lin),
    names_to = "metric",
    values_to = "similarity",
    names_prefix = "sim_"
  )

metrics <- unique(similarity_data$metric)

cliffs_results <- data.frame(
  metric = character(),
  cliffs_delta = numeric(),
  stringsAsFactors = FALSE
)

for (metric_name in metrics) {
  # Get data for current metric
  metric_data <- similarity_data[similarity_data$metric == metric_name, ]
  intra_vals <- metric_data$similarity[metric_data$class == "same_module"]
  intra_vals <- intra_vals[!is.na(intra_vals)]
  inter_vals <- metric_data$similarity[metric_data$class == "different_module"]
  inter_vals <- inter_vals[!is.na(inter_vals)]

  # Calculate Cliff's delta using effsize package
  cliff_result <- cliff.delta(intra_vals, inter_vals)
  delta_value <- cliff_result$estimate

  cliffs_results <- rbind(cliffs_results, data.frame(
    metric = metric_name,
    cliffs_delta = delta_value
  ))
}

save(cliffs_results, file = "GO_terms/go_terms_cliffs_results.rda")

for (metric_name in metrics) {
  metric_data <- similarity_data[similarity_data$metric == metric_name, ]
  na_ratio <- sum(is.na(metric_data$similarity))/nrow(metric_data)
  text <- sprintf("metric %s: %f NA values", metric_name, na_ratio)
  print(text)
}

# Remove the row with NA
temp_data_remove_na <- temp_data[!is.na(temp_data$sim_wang),]

# Calculate Cliff's delta for each metric
library(effectsize)

similarity_data <- temp_data_remove_na |>
  dplyr::select(from, to, sim_embedding, sim_jaccard, sim_op, sim_kappa, sim_dice,
                sim_wang, sim_resnik, sim_rel, sim_jiang, sim_lin, class) %>%
  tidyr::pivot_longer(
    cols = c(sim_embedding, sim_jaccard, sim_op, sim_kappa, sim_dice,
             sim_wang, sim_resnik, sim_rel, sim_jiang, sim_lin),
    names_to = "metric",
    values_to = "similarity",
    names_prefix = "sim_"
  )

metrics <- unique(similarity_data$metric)

cliffs_results <- data.frame(
  metric = character(),
  cliffs_delta = numeric(),
  stringsAsFactors = FALSE
)

for (metric_name in metrics) {
  # Get data for current metric
  metric_data <- similarity_data[similarity_data$metric == metric_name, ]
  intra_vals <- metric_data$similarity[metric_data$class == "same_module"]
  inter_vals <- metric_data$similarity[metric_data$class == "different_module"]

  # Calculate Cliff's delta using effsize package
  cliff_result <- cliff.delta(intra_vals, inter_vals)
  delta_value <- cliff_result$estimate

  cliffs_results <- rbind(cliffs_results, data.frame(
    metric = metric_name,
    cliffs_delta = delta_value
  ))
}

save(cliffs_results, file = "GO_terms/go_term_cliffs_results_remove_na.rda")
