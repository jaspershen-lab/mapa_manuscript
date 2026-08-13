## Benchmarking on the informative GO set benchmark (k = 20)
## Step 9: LLM/annotation benchmarking of MAPA, aPEAR, PAVER and enrichplot.
## All four tools annotate the GROUND-TRUTH modules (not their own clustering
## results). The ground-truth annotation of each module is
## "[anchor_name]. [anchor_definition]" of its anchor GO term.
## Annotation generation follows
## 1_code/2_control_data/05_benchmarking/04_compare_annotation.R:
##   - MAPA:       llm_interpret_module() (gpt-4o-mini-2024-07-18 +
##                 text-embedding-3-small); combined label = name + summary
##   - enrichplot: over-represented word tags (get_wordcloud, nWords = 4)
##   - aPEAR:      pagerank representative pathway on the biotext embedding
##                 similarity matrix; combined label = rep. term name + definition
##   - PAVER:      term nearest to the cluster centroid embedding;
##                 combined label = rep. term name + definition
## Similarity to the ground truth: cosine similarity of
## text-embedding-3-small embeddings of the combined labels.

library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(org.Hs.eg.db)

api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "") stop("OPENAI_API_KEY not found in ~/.Renviron")

res_dir <- "3_data_analysis/11_informative_go_benchmarking"

mapa_input <- readr::read_csv("2_data/informative_go_set/k_20/mapa_input.csv",
                              show_col_types = FALSE)
modules_dt <- readr::read_csv("2_data/informative_go_set/k_20/modules.csv",
                              show_col_types = FALSE)
ground_truth_dt <- readr::read_csv("2_data/informative_go_set/k_20/ground_truth.csv",
                                   show_col_types = FALSE)

## Ground-truth annotation: "[anchor_name]. [anchor_definition]"
gt_annotation <-
  modules_dt |>
  dplyr::transmute(
    module_id = module_id,
    gt_label = paste0(anchor_name, ". ", anchor_definition)
  )
stopifnot(nrow(gt_annotation) == 57)

## Numeric ground-truth labels (same coding as the other scripts)
ground_truth_dt <-
  ground_truth_dt |>
  dplyr::mutate(ground_truth_label = as.numeric(factor(module_id)))

# 1. MAPA: LLM interpretation of the ground-truth modules ====
## Build a functional_module object whose modules are exactly the ground-truth
## modules: replace the similarity matrix by the ground-truth block matrix
## (1 within a module, 0 between modules); louvain then recovers each module
## as one connected clique, and mapa itself assembles the module result.
load("3_data_analysis/10_sensitivity_res_GO/bioembed_sim_res.rda")
load(file.path(res_dir, "pathway_genes.rda"))

## Real gene annotation so that the LLM prompts contain real gene symbols
go_result <- bioembed_sim_res$enriched_pathway@enrichment_go_result@result
go_result$geneID <- sapply(pathway_genes[go_result$ID], paste, collapse = "/")
bioembed_sim_res$enriched_pathway@enrichment_go_result@result <- go_result

all_accessions <- unique(unlist(pathway_genes))
uniprot_map <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db,
                        keys = all_accessions,
                        keytype = "UNIPROT",
                        columns = c("ENSEMBL", "ENTREZID", "SYMBOL"))
) |>
  dplyr::distinct(UNIPROT, .keep_all = TRUE)
uniprot_map <- uniprot_map[match(all_accessions, uniprot_map$UNIPROT), ]

bioembed_sim_res$enriched_pathway@variable_info <-
  data.frame(
    ensembl = ifelse(is.na(uniprot_map$ENSEMBL), all_accessions, uniprot_map$ENSEMBL),
    entrezid = ifelse(is.na(uniprot_map$ENTREZID), all_accessions, uniprot_map$ENTREZID),
    uniprot = all_accessions,
    symbol = ifelse(is.na(uniprot_map$SYMBOL), all_accessions, uniprot_map$SYMBOL),
    variable_id = paste0("gene_", seq_along(all_accessions))
  )

bioembed_sim_res$enriched_pathway@enrichment_kegg_result <- NULL
bioembed_sim_res$enriched_pathway@enrichment_reactome_result <- NULL

## Ground-truth block similarity matrix
ids <- rownames(bioembed_sim_res$sim_matrix)
gt_label_vec <- ground_truth_dt$ground_truth_label[match(ids, ground_truth_dt$go_id)]
stopifnot(!anyNA(gt_label_vec))
gt_block_matrix <- outer(gt_label_vec, gt_label_vec,
                         function(a, b) as.numeric(a == b))
rownames(gt_block_matrix) <- colnames(gt_block_matrix) <- ids
bioembed_sim_res$sim_matrix <- gt_block_matrix

set.seed(23)
fm_gt <- get_functional_modules(
  object = bioembed_sim_res,
  sim.cutoff = 0.5,
  cluster_method = "louvain"
)

## Verify the modules match the ground truth exactly (ARI = 1)
mapa_module_map <-
  fm_gt@merged_module$result_with_module |>
  dplyr::select(node, module) |>
  dplyr::left_join(ground_truth_dt[, c("go_id", "module_id")],
                   by = c("node" = "go_id"))
stopifnot(
  mclust::adjustedRandIndex(mapa_module_map$module, mapa_module_map$module_id) == 1
)
module2gt <- mapa_module_map |>
  dplyr::distinct(module, module_id)
stopifnot(nrow(module2gt) == 57)
save(fm_gt, file = file.path(res_dir, "fm_gt.rda"))

## LLM interpretation (same settings as the case study / control data:
## gpt-4o-mini-2024-07-18 + text-embedding-3-small)
dir.create(file.path(res_dir, "embedding_output"),
           showWarnings = FALSE, recursive = TRUE)
llm_interpreted_gt <- llm_interpret_module(
  object = fm_gt,
  module_content_number_cutoff = 1,
  llm_model = "gpt-4o-mini-2024-07-18",
  embedding_model = "text-embedding-3-small",
  api_key = api_key,
  embedding_output_dir = file.path(res_dir, "embedding_output"),
  orgdb = org.Hs.eg.db,
  output_prompt = FALSE,
  api_provider = "openai"
)

save(llm_interpreted_gt, file = file.path(res_dir, "llm_interpreted_gt.rda"))

mapa_annotation_res <-
  purrr::map_dfr(names(llm_interpreted_gt@llm_module_interpretation), function(m) {
    gn <- llm_interpreted_gt@llm_module_interpretation[[m]]$generated_name
    data.frame(module = m,
               mapa_name = ifelse(is.null(gn$module_name), NA_character_, gn$module_name),
               mapa_summary = ifelse(is.null(gn$summary), NA_character_, gn$summary))
  }) |>
  dplyr::left_join(module2gt, by = "module") |>
  dplyr::mutate(mapa_combined_label = paste0(mapa_name, ". ", mapa_summary)) |>
  dplyr::select(module_id, mapa_cluster_label = mapa_name, mapa_combined_label)
stopifnot(nrow(mapa_annotation_res) == 57, !anyNA(mapa_annotation_res$mapa_combined_label))

# 2. enrichplot: over-represented word tags ====
pathway_info_with_cluster_res <-
  mapa_input |>
  dplyr::transmute(id = go_id, Description = name) |>
  dplyr::left_join(ground_truth_dt[, c("go_id", "module_id", "ground_truth_label")],
                   by = c("id" = "go_id"))

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
  module_id_i <- unique(pathway_info_with_cluster_res$module_id[
    pathway_info_with_cluster_res$ground_truth_label == i])
  enrichplot_annotation_res <- rbind(enrichplot_annotation_res,
                                     data.frame(module_id = module_id_i,
                                                enrich_plot_cluster_label = cluster_label_i))
}
## enrichplot only produces a word-tag label -> combined label = label
enrichplot_annotation_res <-
  enrichplot_annotation_res |>
  dplyr::mutate(enrich_plot_combined_label = enrich_plot_cluster_label)

# 3. aPEAR: pagerank representative pathway ====
load("3_data_analysis/10_sensitivity_res_GO/embedding_sim_matrix.rda")

library(foreach)
library(tibble)
library(data.table)

named_clusters <- ground_truth_dt$ground_truth_label
names(named_clusters) <- ground_truth_dt$go_id
named_clusters <- named_clusters[rownames(embedding_sim_matrix)]

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

apear_rep_terms <- clusterNamesPagerank(sim = embedding_sim_matrix,
                                        clusters = named_clusters)

apear_annotation_res <-
  data.frame(go_id = names(apear_rep_terms),
             rep_id = unname(apear_rep_terms)) |>
  dplyr::left_join(ground_truth_dt[, c("go_id", "module_id")], by = "go_id") |>
  dplyr::distinct(module_id, rep_id) |>
  dplyr::left_join(mapa_input[, c("go_id", "name", "definition")],
                   by = c("rep_id" = "go_id")) |>
  dplyr::transmute(
    module_id = module_id,
    apear_cluster_label = name,
    apear_combined_label = paste0(name, ". ", definition)
  )
stopifnot(nrow(apear_annotation_res) == 57)

# 4. PAVER: term nearest to the cluster centroid embedding ====
load(file.path(res_dir, "PAVER_02_embedding_matrix.rda"))

paver_embedding_df <- embedding_matrix |> as.data.frame()
paver_embedding_df <- paver_embedding_df |>
  dplyr::mutate(UniqueID = rownames(paver_embedding_df)) |>
  dplyr::select(UniqueID, dplyr::everything())

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

  D_sim[D_sim < 0] <- 0
  D_sim[is.na(D_sim)] <- 0

  D_sim
}

nearestpathways <- purrr::map_chr(levels(clustering$Cluster), function(i) {
  cluster_embeddings <- clustering %>%
    dplyr::filter(.data$Cluster == i) %>%
    dplyr::select(-.data$Cluster) %>%
    tibble::column_to_rownames("UniqueID")

  avg_cluster_embedding <- avg_cluster_embeddings %>%
    dplyr::filter(.data$Cluster == i) %>%
    dplyr::select(-.data$Cluster) %>%
    as.numeric()

  dissimilarities <- cluster_embeddings %>%
    dplyr::group_nest(dplyr::row_number()) %>%
    dplyr::pull() %>%
    purrr::map_dbl(~ cosine_dissimilarity(as.matrix(rbind(
      avg_cluster_embedding, .x
    ))))

  rownames(cluster_embeddings)[which.min(dissimilarities)]
})

paver_annotation_res <-
  data.frame(ground_truth_label = as.numeric(levels(clustering$Cluster)),
             rep_id = nearestpathways) |>
  dplyr::left_join(ground_truth_dt |>
                     dplyr::distinct(module_id, ground_truth_label),
                   by = "ground_truth_label") |>
  dplyr::left_join(mapa_input[, c("go_id", "name", "definition")],
                   by = c("rep_id" = "go_id")) |>
  dplyr::transmute(
    module_id = module_id,
    paver_cluster_label = name,
    paver_combined_label = paste0(name, ". ", definition)
  )
stopifnot(nrow(paver_annotation_res) == 57)

# 5. Combine all annotations ====
all_annotation_result <-
  gt_annotation |>
  dplyr::left_join(mapa_annotation_res, by = "module_id") |>
  dplyr::left_join(apear_annotation_res, by = "module_id") |>
  dplyr::left_join(paver_annotation_res, by = "module_id") |>
  dplyr::left_join(enrichplot_annotation_res, by = "module_id")
stopifnot(nrow(all_annotation_result) == 57, !anyNA(all_annotation_result))

save(all_annotation_result, file = file.path(res_dir, "all_annotation_result.rda"))

# 6. Embeddings of combined labels and ground truth ====
long_data <-
  all_annotation_result |>
  tidyr::pivot_longer(
    cols = c(mapa_combined_label, apear_combined_label,
             paver_combined_label, enrich_plot_combined_label),
    names_to = "method",
    values_to = "combined_label"
  ) |>
  dplyr::mutate(method = sub("_combined_label", "_cluster_label", method)) |>
  dplyr::select(module_id, gt_label, method, combined_label)

long_data$combined_emb <- vector("list", nrow(long_data))
for (i in 1:nrow(long_data)) {
  long_data$combined_emb[[i]] <- mapa::get_embedding(
    chunk = long_data$combined_label[i],
    api_key = api_key,
    model_name = "text-embedding-3-small",
    api_provider = "openai"
  )
}

gt_annotation$gt_emb <- vector("list", nrow(gt_annotation))
for (i in 1:nrow(gt_annotation)) {
  gt_annotation$gt_emb[[i]] <- mapa::get_embedding(
    chunk = gt_annotation$gt_label[i],
    api_key = api_key,
    model_name = "text-embedding-3-small",
    api_provider = "openai"
  )
}

save(long_data, file = file.path(res_dir, "annotation_long_data.rda"))
save(gt_annotation, file = file.path(res_dir, "gt_annotation.rda"))

# 7. Cosine similarity between tool annotations and the ground truth ====
cosine_similarity <- function(vec1, vec2) {
  dot_product <- sum(vec1 * vec2)
  norm_vec1 <- sqrt(sum(vec1^2))
  norm_vec2 <- sqrt(sum(vec2^2))

  if (norm_vec1 == 0 | norm_vec2 == 0) {
    return(0)
  }

  return(dot_product / (norm_vec1 * norm_vec2))
}

annotation_similarity_results <-
  long_data |>
  dplyr::left_join(gt_annotation[, c("module_id", "gt_emb")], by = "module_id") |>
  dplyr::mutate(
    cosine_sim = purrr::map2_dbl(combined_emb, gt_emb, ~ cosine_similarity(.x, .y))
  ) |>
  dplyr::select(module_id, method, combined_label, cosine_sim)

save(annotation_similarity_results,
     file = file.path(res_dir, "annotation_similarity_results.rda"))

## Export tables
rio::export(all_annotation_result |> dplyr::select(-dplyr::ends_with("_emb")),
            file = file.path(res_dir, "annotation_result.xlsx"))
rio::export(annotation_similarity_results |> dplyr::select(-combined_label),
            file = file.path(res_dir, "similarity_gt_tool.xlsx"))

## Summary
annotation_similarity_results |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    median_cosine_sim = median(cosine_sim),
    mean_cosine_sim = mean(cosine_sim),
    n = dplyr::n()
  ) |>
  dplyr::arrange(dplyr::desc(median_cosine_sim)) |>
  print()
