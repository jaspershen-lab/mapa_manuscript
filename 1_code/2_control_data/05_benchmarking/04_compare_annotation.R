library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

setwd("3_data_analysis/02_control_data/05_benchmarking")
load("enriched_result.rda")
load("ground_truth_dt.rda")

pathway_info_with_cluster_res <- enriched_result@result |>
  select(ID, Description) |>
  rename(id = ID) |>
  left_join(ground_truth_dt, by = "id")

# based on previous results -> make sure that mapa's clustering is reliable
# -> here the annotation based on mapa's clustering result (human expert is the same)

# enrichplot ====
get_wordcloud <- function(cluster, pathway_info_with_cluster_res, nWords){
  words <- pathway_info_with_cluster_res$Description %>%
    gsub(" in ", " ", .) %>%
    gsub(" [0-9]+ ", " ", .) %>%
    gsub("^[0-9]+ ", "", .) %>%
    gsub(" [0-9]+$", "", .) %>%
    gsub(" [A-Za-z] ", " ", .) %>%
    gsub("^[A-Za-z] ", "", .) %>%
    gsub(" [A-Za-z]$", "", .) %>%
    gsub(" / ", " ", .) %>%
    gsub(" and ", " ", .) %>%
    gsub(" of ", " ", .) %>%
    gsub(",", " ", .) %>%
    gsub(" - ", " ", .)
  net_tot <- length(words)

  clusters <- unique(pathway_info_with_cluster_res$ground_truth_label)
  words_i <- words[which(pathway_info_with_cluster_res$ground_truth_label == cluster)]

  sel_tot <- length(words_i)
  sel_w <- get_word_freq(words_i)
  net_w_all <- get_word_freq(words)
  net_w <- net_w_all[names(sel_w)]
  tag_size <- (sel_w/sel_tot)/(net_w/net_tot)
  tag_size <- tag_size[order(tag_size, decreasing = TRUE)]
  nWords <- min(nWords, length(tag_size))
  tag <- names(tag_size[seq_len(nWords)])

  # Order of words
  dada <- strsplit(words_i, " ")
  len <- vapply(dada, length, FUN.VALUE=1)
  rank <- NULL
  for(i in seq_len(sel_tot)) {
    rank <- c(rank, seq_len(len[i]))
  }

  word_data <- data.frame(word = unlist(dada), rank = rank)
  word_rank1 <- stats::aggregate(rank ~ word, data = word_data, sum)
  rownames(word_rank1) <- word_rank1[, 1]

  word_rank1 <- word_rank1[names(sel_w), ]
  # Get an average ranking order
  word_rank1[, 2] <- word_rank1[, 2]/as.numeric(sel_w)
  tag_order <- word_rank1[tag, ]
  tag_order <- tag_order[order(tag_order[, 2]), ]
  tag_clu_i <- paste(tag_order$word, collapse=" ")
  tag_clu_i
}

get_word_freq <- function(wordd){
  dada <- strsplit(wordd, " ")
  didi <- table(unlist(dada))
  didi <- didi[order(didi, decreasing = TRUE)]
  # Get the number of each word
  word_name <- names(didi)
  fun_num_w <- function(ww){
    sum(vapply(dada, function(w){ww %in% w}, FUN.VALUE = 1))
  }
  word_num <- vapply(word_name, fun_num_w, FUN.VALUE = 1)
  word_w <- word_num[order(word_num, decreasing = TRUE)]
}

enrichplot_annotation_res <- data.frame()

for (i in unique(pathway_info_with_cluster_res$ground_truth_label)) {
  cluster_label_i <- get_wordcloud(cluster = i,
                                   pathway_info_with_cluster_res = pathway_info_with_cluster_res,
                                   nWords = 4)
  enrichplot_annotation_res <- rbind(enrichplot_annotation_res,
                                     data.frame(cluster = i,
                                                enrich_plot_cluster_label = cluster_label_i))
}

# aPEAR ====
# sim: a similarity matrix used to detect the clusters
# clusters: a vector of clusters, with names indicating the pathway name
setwd(get_project_wd())
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_sim_matrix.rda")

rowname_order <- rownames(embedding_sim_matrix)
reordered_name <- pathway_info_with_cluster_res$Description[match(rowname_order, pathway_info_with_cluster_res$id)]
rownames(embedding_sim_matrix) <- reordered_name
colnames(embedding_sim_matrix) <- reordered_name

named_clusters <- pathway_info_with_cluster_res$ground_truth_label
names(named_clusters) <- pathway_info_with_cluster_res$Description

library(foreach)
library(tibble)
library(data.table)

cluster_size <- data.frame(table(named_clusters))
clusters_rm_singleton <- as.numeric(cluster_size$named_clusters[cluster_size$Freq != 1])
filtered_clusters <- named_clusters[named_clusters %in% clusters_rm_singleton]


mapClusterNames <- function(scores, clusters) {
  clusterNames <- foreach(cluster = unique(clusters), .combine = rbind) %do% {
    name <- scores[ names(scores) %in% names(clusters[ clusters == cluster ]) ] %>%
      which.max %>%
      names

    data.table(ClusterID = cluster, Name = name)
  }

  as.data.table(clusters, keep.rownames = TRUE) %>%
    merge(clusterNames, by.x = 'clusters', by.y = 'ClusterID') %>%
    .[ , list(rn, Name) ] %>%
    deframe
}

clusterNamesPagerank <- function(sim, clusters) {
  stopifnot(rownames(sim) == colnames(sim))
  stopifnot(nrow(sim) == ncol(sim))

  paths <- rownames(sim)
  edges <- list()
  counter <- 1

  for (i in 1:(nrow(sim) - 1)) {
    for (j in (i + 1):ncol(sim)) {
      value <- sim[ i, j ]

      clusteri <- clusters[ paths[ i ] ]
      clusterj <- clusters[ paths[ j ] ]
      if (!anyNA(c(clusteri, clusterj)) && clusteri == clusterj) {
        edges[[counter]] <- data.table(from = paths[ i ], to = paths[ j ], weight = value)
        counter <- counter + 1
      }
    }
  }

  edges <- rbindlist(edges)

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  scores <- igraph::page.rank(g)$vector

  mapClusterNames(scores, clusters)
}

apear_annotation_res_without_singleton <- clusterNamesPagerank(sim = embedding_sim_matrix,
                                                               clusters = filtered_clusters)
# For singleton, aPEAR does not give annotation, just use the original pathway name
apear_annotation_res <- data.frame(Description = names(apear_annotation_res_without_singleton),
                                   apear_cluster_label = unname(apear_annotation_res_without_singleton))
apear_annotation_res <- apear_annotation_res |>
  right_join(pathway_info_with_cluster_res, by = "Description")
apear_annotation_res <- apear_annotation_res |>
  mutate(apear_cluster_label = if_else(is.na(apear_cluster_label), Description, apear_cluster_label)) |>
  select(ground_truth_label, apear_cluster_label) |>
  distinct(ground_truth_label, .keep_all = TRUE) |>
  rename(cluster = ground_truth_label)

# PAVER ====
# Average the embeddings within each cluster
setwd(get_project_wd())
load("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/biotext_embedding/embedding_matrix.rda")

rownames(embedding_matrix) <- rownames(embedding_sim_matrix)
paver_embedding_df <- embedding_matrix |> as.data.frame()
paver_embedding_df <- paver_embedding_df |> mutate(UniqueID = rownames(paver_embedding_df)) |> select(UniqueID, everything())

clustering <- named_clusters |>
  tibble::enframe(name = "UniqueID", value = "Cluster") |>
  dplyr::mutate(Cluster = as.factor(.data$Cluster)) |>
  dplyr::inner_join(paver_embedding_df, by = "UniqueID")

avg_cluster_embeddings <- clustering |>
  dplyr::group_by(.data$Cluster) |>
  dplyr::summarise(dplyr::across(!.data$UniqueID, mean))

cosine_dissimilarity <- function(mat, root = FALSE) {
  sim <- mat / sqrt(rowSums(mat * mat))
  sim <- sim %*% t(sim)

  if (root == T) {
    D_sim <- stats::as.dist(suppressWarnings(sqrt(1 - sim)))
  } else {
    D_sim <- stats::as.dist(1 - sim)
  }

  D_sim[D_sim < 0] <- 0 # Precision errors can lead to negative or NA distances
  D_sim[is.na(D_sim)] <- 0

  D_sim
}

# Identify pathways most similar to the average of each cluster i.e., it's "theme"
nearestpathways <- purrr::map_chr(1:nlevels(clustering$Cluster), function(i) {
  # Grab all embeddings within the cluster
  cluster_embeddings <- clustering %>%
    dplyr::filter(.data$Cluster == i) %>%
    dplyr::select(-.data$Cluster) %>%
    tibble::column_to_rownames("UniqueID")

  # Grab the average embedding of the cluster
  avg_cluster_embedding <- avg_cluster_embeddings %>%
    dplyr::filter(.data$Cluster == i) %>%
    dplyr::select(-.data$Cluster) %>%
    as.numeric()

  # Calculate cosine dissimilarity for each embedding in the cluster
  dissimilarities <- cluster_embeddings %>%
    dplyr::group_nest(dplyr::row_number()) %>%
    dplyr::pull() %>%
    purrr::map_dbl(~ cosine_dissimilarity(as.matrix(rbind(
      avg_cluster_embedding, .x
    ))))

  # Return the name of the embedding with smallest dissimilarity
  rownames(cluster_embeddings)[which.min(dissimilarities)]
}) %>%
  tibble::enframe(name = "Cluster", value = "UniqueID")

paver_annotation_res <- nearestpathways |>
  rename(cluster = Cluster,
         paver_cluster_label = UniqueID)

all_annotation_result <- enrichplot_annotation_res |>
  left_join(apear_annotation_res, by = "cluster") |>
  left_join(paver_annotation_res, by = "cluster")

save(all_annotation_result, file = "3_data_analysis/02_control_data/05_benchmarking/comparison_result/all_annotation_result.rda")

# Get embedding for name ====
setwd("3_data_analysis/02_control_data/05_benchmarking")
load("feifan-combine_scores.Rdata")
load("comparison_result/all_annotation_result.rda")

library(tidyr)

long_data <- all_annotation_result %>%
  pivot_longer(
    cols = c(enrich_plot_cluster_label, apear_cluster_label, paver_cluster_label),
    names_to = "method",
    values_to = "label"
  )

library(mapa)

long_data <- long_data |> mutate(name_emb = NA)

for (i in 1:nrow(long_data)) {
  label <- long_data$label[i]
  name_emb <- mapa::get_embedding(chunk = label,
                                  api_key = api_key,
                                  model_name = "text-embedding-3-small",
                                  api_provider = "openai")
  list_emb <- list(name_emb)
  long_data$name_emb[i] <- list_emb
}

human_expert <- unique(combined_scores$expert_name[combined_scores$expert_name != "LLM_Analysis"])

mapa_and_h_expert_result <- combined_scores |>
  select(module, expert_name, module_name, name_emb) |>
  mutate(cluster = as.numeric(sub("Functional_module_", "", module))) |>
  mutate(method = if_else(expert_name == "LLM_Analysis",
                          "mapa_cluster_label",
                          expert_name)) |>
  select(cluster, method, module_name, name_emb) |>
  rename(label = module_name)

all_result_long_data <- rbind(mapa_and_h_expert_result, long_data)
filtered_all_long_data <- all_result_long_data |> filter(!(cluster %in% c(1,2,3,6)))

save(filtered_all_long_data, file = "comparison_result/filtered_all_long_data.rda")

# Calculate cosine sim ====
# Function to calculate cosine similarity between two vectors
load("comparison_result/filtered_all_long_data.rda")
cosine_similarity <- function(vec1, vec2) {
  dot_product <- sum(vec1 * vec2)
  norm_vec1 <- sqrt(sum(vec1^2))
  norm_vec2 <- sqrt(sum(vec2^2))

  if (norm_vec1 == 0 | norm_vec2 == 0) {
    return(0)
  }

  return(dot_product / (norm_vec1 * norm_vec2))
}

# Separate human expert data and tool data
human_data <- filtered_all_long_data %>%
  filter(method %in% human_expert) %>%
  select(cluster, method, name_emb) %>%
  rename(expert_name = method, expert_embedding = name_emb, expert_cluster = cluster)

tool_data <- filtered_all_long_data %>%
  filter(!method %in% human_expert) %>%
  select(cluster, method, name_emb) %>%
  rename(tool_name = method, tool_embedding = name_emb, tool_cluster = cluster)

# Calculate cosine similarity for each tool-expert pair within the same cluster
similarity_results <- tool_data %>%
  crossing(human_data) %>%
  filter(tool_cluster == expert_cluster) %>%  # Only compare within same cluster
  mutate(
    cosine_sim = map2_dbl(tool_embedding, expert_embedding, ~ cosine_similarity(.x, .y))
  )

save(similarity_results, file = "comparison_result/similarity_results.rda")

# Calculate cosine similarity between combined summaries ====
setwd("3_data_analysis/02_control_data/05_benchmarking")
load("comparison_result/all_annotation_result.rda")
load("go_pathway_info.rda")
load("kegg_pathway_info.rda")
load("reactome_pathway_info.rda")
## Generate combined summaries
### aPEAR and PAVER select the representative pathway
### -> combined summary = representative pathway name + description
all_pathway_info <- c(go_pathway_info,
                      kegg_pathway_info,
                      reactome_pathway_info)

combined_methods <- c("mapa_cluster_label", human_expert)
all_result_long_data <- all_result_long_data |> mutate(combined_label = NA,
                                                       combined_emb = NA)

combined_scores <- combined_scores |>
  dplyr::mutate(cluster = as.numeric(sub("Functional_module_", "", module)))

for (i in 1:nrow(all_result_long_data)) {
  method_i <- all_result_long_data$method[i]
  cluster_i <- all_result_long_data$cluster[i]
  if (method_i %in% combined_methods) {
    if (method_i == "mapa_cluster_label") {
      index_i <- which(combined_scores$expert_name == "LLM_Analysis" & combined_scores$cluster == cluster_i)
    } else {
      index_i <- which(combined_scores$expert_name == method_i & combined_scores$cluster == cluster_i)
    }
    combined_label_i <- paste0(combined_scores$module_name[index_i], ". ", combined_scores$module_description[index_i])

    all_result_long_data$combined_label[i] <- combined_label_i
    all_result_long_data$combined_emb[i] <- combined_scores$combined_emb[index_i]
  } else if (method_i %in% c("apear_cluster_label", "paver_cluster_label")) {
    term_name_i <- all_result_long_data$label[i]
    for (m in seq_along(all_pathway_info)) {
      if (all_pathway_info[[m]]$term_name == term_name_i) {
        combined_label_i <- paste0(all_pathway_info[[m]]$term_name, ". ", all_pathway_info[[m]]$term_definition)
        combined_emb_i <- list(mapa::get_embedding(
          chunk = combined_label_i,
          api_key = api_key,
          model_name = "text-embedding-3-small",
          api_provider = "openai"
        ))
      }
    }

    all_result_long_data$combined_label[i] <- combined_label_i
    all_result_long_data$combined_emb[i] <- combined_emb_i
  } else if (method_i == "enrich_plot_cluster_label") {
    all_result_long_data$combined_label[i] <- all_result_long_data$label[i]
    all_result_long_data$combined_emb[i] <- all_result_long_data$name_emb[i]
  }
}

save(all_result_long_data, file = "comparison_result/all_result_long_data.rda")

filtered_all_result_long_data <- all_result_long_data |> filter(!(cluster %in% c(1,2,3,6)))
save(filtered_all_result_long_data, file = "comparison_result/filtered_all_result_long_data.rda")

## Calculate cosine sim between combined labels
# Separate human expert data and tool data
human_data <- filtered_all_result_long_data %>%
  filter(method %in% human_expert) %>%
  select(cluster, method, combined_emb) %>%
  rename(expert_name = method, expert_embedding = combined_emb, expert_cluster = cluster)

tool_data <- filtered_all_result_long_data %>%
  filter(!method %in% human_expert) %>%
  select(cluster, method, combined_emb) %>%
  rename(tool_name = method, tool_embedding = combined_emb, tool_cluster = cluster)

# Calculate cosine similarity for each tool-expert pair within the same cluster
combined_similarity_results <- tool_data %>%
  crossing(human_data) %>%
  filter(tool_cluster == expert_cluster) %>%  # Only compare within same cluster
  mutate(
    cosine_sim = map2_dbl(tool_embedding, expert_embedding, ~ cosine_similarity(.x, .y))
  )

save(combined_similarity_results, file = "comparison_result/combined_similarity_results.rda")

tool_names <- c("mapa_cluster_label", "apear_cluster_label", "enrich_plot_cluster_label", "paver_cluster_label")

annotation_result <-
  filtered_all_result_long_data |>
  dplyr::filter(method %in% tool_names) |>
  dplyr::select(cluster, method, combined_label) |>
  dplyr::rename(annotation = combined_label,
                functional_module = cluster) |>
  dplyr::mutate(method = sub("_cluster_label", "", method)) |>
  dplyr::arrange(functional_module)

sim_hm_vs_tools <-
  combined_similarity_results |>
  dplyr::select(tool_cluster, tool_name, expert_name, cosine_sim) |>
  dplyr::rename(functional_module = tool_cluster) |>
  dplyr::mutate(method = sub("_cluster_label", "", tool_name)) |>
  dplyr::select(-tool_name) |>
  left_join(hm_name_id, by = "expert_name")

sim_hm_vs_tools <- sim_hm_vs_tools |> dplyr::select(functional_module, expert_id, method, cosine_sim)

library(rio)
export(all_clustering_res, file = "clustering_result.xlsx")
export(annotation_result, file = "annotation_result.xlsx")
export(sim_hm_vs_tools, file = "similarity_hm_tool.xlsx")
