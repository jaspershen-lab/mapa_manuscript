library(dplyr)
library(purrr)
library(httr2)
library(org.Hs.eg.db)
library(reactome.db)
library(AnnotationDbi)
library(rvest)
library(rentrez)
library(KEGGREST)
library(rbioapi)
library(httr)
library(httr2)
library(curl)
library(jsonlite)
library(rtiktoken)
library(bayesbio)
library(GOSemSim)
library(tidygraph)
library(igraph)

control_pathways <- control_data

# I. Similarity ====

# 1. Get annotated gene information for pathways ====
## 22 GO terms
go_info <-
  control_pathways %>%
  dplyr::filter(database == "GO") %>%
  dplyr::pull(id) %>%
  get_go_info(include_gene_name = FALSE)

## Test gene name influence ====
test_id <- c(top_dif$from, top_dif$to) %>% unique()
test_go_info <- get_go_info(go_ids = test_id, include_gene_name = TRUE)
only_gname <- combine_gene_name(info = test_go_info, include_gene_name = TRUE)
combine_gene_name <- function(info, include_gene_name = FALSE) {
  info %>%
    purrr::map(
      function(x) {
        if (include_gene_name) {
          text_info <- sprintf("Annotated gene names: %s", x$annotated_genename)
        } else {
          text_info <- sprintf("Pathway name: %s\nDefinition: %s", x$term_name, x$term_definition)
        }

        text <- list(
          "id" = x$id,
          "text_info" = text_info
        )

        return(text)
      }
    )
}

embedding_matrix_only_gname <- get_embedding_matrix(text = only_gname, include_gene_name = TRUE, api_provider = "openai", text_embedding_model = "text-embedding-3-small", api_key = api_key)
test_sim_matrix <- calculate_cosine_sim(m = embedding_matrix_only_gname)
test_semantic_sim_df <-
  as.data.frame.table(test_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

test_compare <- left_join(test_semantic_sim_df, go_combined_sim_df, by = c("from", "to"))
## Test end =====

## 14 KEGG pathways
kegg_info <-
  control_pathways %>%
  dplyr::filter(database == "KEGG") %>%
  dplyr::pull(id) %>%
  get_kegg_info(include_gene_name = FALSE)

## 8 Ractome pathways
reactome_info <-
  control_pathways %>%
  dplyr::filter(database == "Reactome") %>%
  dplyr::pull(id) %>%
  get_reactome_info(include_gene_name = FALSE)

all_text_info <- c(go_info, kegg_info, reactome_info)

all_combined_info <- combine_info(info = all_text_info, include_gene_name = FALSE)

# 2. Get semantic similarity using embedding model ====
embedding_matrix <- get_embedding_matrix(text = all_combined_info, include_gene_name = FALSE, api_provider = "openai", text_embedding_model = "text-embedding-3-small", api_key = api_key)

sim_matrix <- calculate_cosine_sim(m = embedding_matrix)

semantic_sim_df <-
  as.data.frame.table(sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

# 3. Get similarity based on gene overlap ====
## Get gene id for each pathway
go_id <-
  control_pathways %>%
  dplyr::filter(database == "GO") %>%
  dplyr::pull(id)

go2egs <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go_id, columns = "ENTREZID", keytype = "GOALL"))

go_gene_lists <- list()

go_pathways <-
  control_pathways %>%
  dplyr::filter(database == "GO")

for (i in 1:nrow(go_pathways)) {
  pathway_id <- go_pathways$id[i]

  # Add this gene set to the list, using the pathway_id as the name
  go_gene_lists[[pathway_id]] <-
    go2egs %>%
    dplyr::filter(GOALL == pathway_id) %>%
    pull(ENTREZID) %>%
    unique()
}

##
reactome_id <-
  control_pathways %>%
  dplyr::filter(database == "Reactome") %>%
  dplyr::pull(id)

reactome2egs <- AnnotationDbi::select(reactome.db, keys = reactome_id, columns = "ENTREZID", keytype = "PATHID")

reactome_gene_lists <- list()

reactome_pathways <-
  control_pathways %>%
  dplyr::filter(database == "Reactome")

for (i in 1:nrow(reactome_pathways)) {
  pathway_id <- reactome_pathways$id[i]

  # Add this gene set to the list, using the pathway_id as the name
  reactome_gene_lists[[pathway_id]] <-
    reactome2egs %>%
    dplyr::filter(PATHID == pathway_id) %>%
    pull(ENTREZID) %>%
    unique()
}

##
kegg_id <-
  control_pathways %>%
  dplyr::filter(database == "KEGG") %>%
  dplyr::pull(id)

chunk_size <- 10
chunks <- split(kegg_id, ceiling(seq_along(kegg_id) / chunk_size))
kegg_info <- list()
for (i in 1:length(chunks)) {
  sub_kegg_info <-
    KEGGREST::keggGet(dbentries = chunks[[i]]) %>%
    purrr::map(
      function(x) {
        # Initialize kegg2genename as empty
        kegg2genename <- NA

        # Only get annotated Entrez ID if include_gene_name is TRUE
        if ("GENE" %in% names(x)) {
          kegg2genename <-
            x$GENE[seq(1, length(x$GENE), 2)]
        }

        all_info <- list(
          "id" = unname(x$ENTRY),
          "genes" = kegg2genename
        )

        return(all_info)
      }
    )
  kegg_info <- c(kegg_info, sub_kegg_info)
}

kegg_gene_lists <- list()

for (i in 1:length(kegg_info)) {
  pathway_id <- kegg_info[[i]]$id

  kegg_gene_lists[[pathway_id]] <- kegg_info[[i]]$genes
}

all_gene_list <- c(go_gene_lists, kegg_gene_lists, reactome_gene_lists)

## Calculate similarity
jaccard_sim_matrix <- term_similarity_internal(gl = all_gene_list,
                                       measure.method = "jaccard")
jaccard_sim_df <-
  as.data.frame.table(jaccard_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

op_sim_matrix <- term_similarity_internal(gl = all_gene_list,
                                          measure.method = "overlap")
op_sim_df <-
  as.data.frame.table(op_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

kappa_sim_matrix <- term_similarity_internal(gl = all_gene_list,
                                             measure.method = "kappa")

kappa_sim_df <-
  as.data.frame.table(kappa_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

dice_sim_matrix <- term_similarity_internal(gl = all_gene_list,
                                            measure.method = "dice")
dice_sim_df <-
  as.data.frame.table(dice_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

# 4. Calculate GO term similarity based on GO hierarchical structure ====
## BP
go_bp_id <- control_pathways %>% dplyr::filter(database == "GO" & go_ontology == "BP") %>% pull(id)
bp_semgodata <- godata(OrgDb = org.Hs.eg.db, ont = "BP")
bp_sim_matrix <- termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata, method = "Wang")
bp_sim_df <-
  as.data.frame.table(bp_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

## MF
go_mf_id <- control_pathways %>% dplyr::filter(database == "GO" & go_ontology == "MF") %>% pull(id)
mf_semgodata <- godata(OrgDb = org.Hs.eg.db, ont = "MF")
mf_sim_matrix <- termSim(t1 = go_mf_id, t2 = go_mf_id, semData = mf_semgodata, method = "Wang")
mf_sim_df <-
  as.data.frame.table(mf_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

## CC
go_cc_id <- control_pathways %>% dplyr::filter(database == "GO" & go_ontology == "CC") %>% pull(id)
cc_semgodata <- godata(OrgDb = org.Hs.eg.db, ont = "CC")
cc_sim_matrix <- termSim(t1 = go_cc_id, t2 = go_cc_id, semData = cc_semgodata, method = "Wang")
cc_sim_df <-
  as.data.frame.table(cc_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

wang_sim_df <- rbind(bp_sim_df, mf_sim_df, cc_sim_df)


# 5. Compare results ====

## 5.1 Compare sim_semantic with other metrics ====
df_lists <- list(semantic_sim_df, jaccard_sim_df, op_sim_df, kappa_sim_df, dice_sim_df)
suffixes <- c("_semantic", "_jaccard", "_op", "_kappa", "_dice")

combined_sim_df <- combine_similarity_dataframes_tidyverse(
  dataframes = df_lists,
  suffixes = suffixes
)

save(combined_sim_df, file = "control_all_combined_sim_df.RData")

combine_similarity_dataframes_tidyverse <- function(dataframes, suffixes = NULL) {
  # Generate suffixes if not provided
  if (is.null(suffixes)) {
    suffixes <- paste0("_df", seq_along(dataframes))
  }

  # Rename the sim columns
  renamed_dfs <- map2(dataframes, suffixes, function(df, suffix) {
    df %>% rename_with(~paste0("sim", suffix), .cols = "sim")
  })

  # Reduce to merge all dataframes
  result <- purrr::reduce(renamed_dfs, function(x, y) {
    full_join(x, y, by = c("from", "to"))
  })

  return(result)
}

eva_sim_re <- create_single_panel_correlation(combined_sim_df)

pdf("control_similarity_semantic_vs_others.pdf", width=8, height=8)
p <- eva_sim_re$plot
print(p)
dev.off()

create_single_panel_correlation <- function(df) {
  # Required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is needed for this function to work. Please install it.")
  }

  # Calculate Pearson correlations
  cor_jaccard <- cor(df$sim_semantic, df$sim_jaccard, use = "complete.obs", method = "pearson")
  cor_op <- cor(df$sim_semantic, df$sim_op, use = "complete.obs", method = "pearson")
  cor_kappa <- cor(df$sim_semantic, df$sim_kappa, use = "complete.obs", method = "pearson")
  cor_dice <- cor(df$sim_semantic, df$sim_dice, use = "complete.obs", method = "pearson")

  # Create a data frame for correlation results
  cor_results <- data.frame(
    metric = c("Jaccard", "Overlap", "Kappa", "Dice"),
    correlation = c(cor_jaccard, cor_op, cor_kappa, cor_dice)
  )

  # Print the correlation results
  print(cor_results)

  # Reshape data for plotting
  plot_data <- data.frame(
    sim_semantic = rep(df$sim_semantic, 4),
    sim_value = c(df$sim_jaccard, df$sim_op, df$sim_kappa, df$sim_dice),
    metric = factor(rep(c("Jaccard", "Overlap", "Kappa", "Dice"), each = nrow(df)))
  )

  # Define colors for each metric
  metric_colors <- c("Jaccard" = "#e9c46b", "Overlap" = "#e66f51", "Kappa" = "#264653", "Dice" = "#299e8c")

  # Create the plot with all regression lines in one panel
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sim_semantic, y = sim_value, color = metric)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    ggplot2::scale_color_manual(values = metric_colors) +
    ggplot2::labs(
      title = "Correlation between Semantic Similarity and Other Metrics",
      x = "Semantic Similarity",
      y = "Similarity Value",
      color = "Similarity Metric"
    ) +
    ggplot2::theme_minimal() +
    # Add annotation with correlation values
    ggplot2::annotate(
      "text",
      x = rep(min(df$sim_semantic, na.rm = TRUE) + 0.05, 4),
      y = seq(
        from = max(c(df$sim_jaccard, df$sim_op, df$sim_kappa, df$sim_dice), na.rm = TRUE) - 0.05,
        to = max(c(df$sim_jaccard, df$sim_op, df$sim_kappa, df$sim_dice), na.rm = TRUE) - 0.35,
        length.out = 4
      ),
      label = sprintf(
        "%s: r = %.3f",
        c("Jaccard", "Overlap", "Kappa", "Dice"),
        c(cor_jaccard, cor_op, cor_kappa, cor_dice)
      ),
      color = metric_colors,
      hjust = 0,
      fontface = "bold"
    )

  # Return the correlation results and plot
  return(list(correlations = cor_results, plot = p))
}

## 5.2 Compare sim_jaccard with other metrics ====
eva_sim_re_jaccard <- create_single_panel_correlation_jaccard(combined_sim_df)
eva_sim_re_jaccard$plot
create_single_panel_correlation_jaccard <- function(df) {
  # Required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is needed for this function to work. Please install it.")
  }

  # Calculate Pearson correlations
  cor_semantic <- cor(df$sim_semantic, df$sim_jaccard, use = "complete.obs", method = "pearson")
  cor_op <- cor(df$sim_jaccard, df$sim_op, use = "complete.obs", method = "pearson")
  cor_kappa <- cor(df$sim_jaccard, df$sim_kappa, use = "complete.obs", method = "pearson")
  cor_dice <- cor(df$sim_jaccard, df$sim_dice, use = "complete.obs", method = "pearson")

  # Create a data frame for correlation results
  cor_results <- data.frame(
    metric = c("Semantic", "Overlap", "Kappa", "Dice"),
    correlation = c(cor_semantic, cor_op, cor_kappa, cor_dice)
  )

  # Print the correlation results
  print(cor_results)

  # Reshape data for plotting
  plot_data <- data.frame(
    sim_jaccard = rep(df$sim_jaccard, 4),
    sim_value = c(df$sim_semantic, df$sim_op, df$sim_kappa, df$sim_dice),
    metric = factor(rep(c("Semantic", "Overlap", "Kappa", "Dice"), each = nrow(df)))
  )

  # Define colors for each metric
  metric_colors <- c("Semantic" = "#e9c46b", "Overlap" = "#e66f51", "Kappa" = "#264653", "Dice" = "#299e8c")

  # Create the plot with all regression lines in one panel
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sim_jaccard, y = sim_value, color = metric)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    ggplot2::scale_color_manual(values = metric_colors) +
    ggplot2::labs(
      title = "Correlation between Jaccard Similarity and Other Metrics",
      x = "Jaccard Similarity",
      y = "Similarity Value",
      color = "Similarity Metric"
    ) +
    ggplot2::theme_minimal() +
    # Add annotation with correlation values
    ggplot2::annotate(
      "text",
      x = rep(min(df$sim_jaccard, na.rm = TRUE) + 0.05, 4),
      y = seq(
        from = max(c(df$sim_semantic, df$sim_op, df$sim_kappa, df$sim_dice), na.rm = TRUE) - 0.05,
        to = max(c(df$sim_semantic, df$sim_op, df$sim_kappa, df$sim_dice), na.rm = TRUE) - 0.35,
        length.out = 4
      ),
      label = sprintf(
        "%s: r = %.3f",
        c("Semantic", "Overlap", "Kappa", "Dice"),
        c(cor_semantic, cor_op, cor_kappa, cor_dice)
      ),
      color = metric_colors,
      hjust = 0,
      fontface = "bold"
    )

  # Return the correlation results and plot
  return(list(correlations = cor_results, plot = p))
}

## 5.3 Compare sim_wang with other metrics ====
df_lists <- list(wang_sim_df, semantic_sim_df, jaccard_sim_df, op_sim_df, kappa_sim_df, dice_sim_df)
suffixes <- c("_wang", "_semantic", "_jaccard", "_op", "_kappa", "_dice")

go_combined_sim_df <- combine_go_similarity_dataframes_tidyverse(
  dataframes = df_lists,
  suffixes = suffixes
)
save(go_combined_sim_df, file = "control_go_combined_sim_df.RData")
combine_go_similarity_dataframes_tidyverse <-
  function(dataframes, suffixes = NULL) {
    # Generate suffixes if not provided
    if (is.null(suffixes)) {
      suffixes <- paste0("_df", seq_along(dataframes))
    }

    # Rename the sim columns
    renamed_dfs <- map2(dataframes, suffixes, function(df, suffix) {
      df %>% rename_with(~paste0("sim", suffix), .cols = "sim")
    })

    # Reduce to merge all dataframes
    result <- purrr::reduce(renamed_dfs, function(x, y) {
      left_join(x, y, by = c("from", "to"))
    })

    return(result)
  }

eva_sim_re_go <- create_single_panel_correlation_go(go_combined_sim_df)
eva_sim_re_go$plot

pdf("control_similarity_wang_vs_others.pdf", width=8, height=8)
p <- eva_sim_re_go$plot
print(p)
dev.off()

create_single_panel_correlation_go <- function(df) {
  # Required libraries
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is needed for this function to work. Please install it.")
  }

  # Calculate Pearson correlations
  cor_semantic <- cor(df$sim_wang, df$sim_semantic, use = "complete.obs", method = "pearson")
  cor_op <- cor(df$sim_wang, df$sim_op, use = "complete.obs", method = "pearson")
  cor_kappa <- cor(df$sim_wang, df$sim_kappa, use = "complete.obs", method = "pearson")
  cor_dice <- cor(df$sim_wang, df$sim_dice, use = "complete.obs", method = "pearson")
  cor_jaccard <- cor(df$sim_wang, df$sim_jaccard, use = "complete.obs", method = "pearson")

  # Create a data frame for correlation results
  cor_results <- data.frame(
    metric = c("Semantic", "Overlap", "Kappa", "Dice", "Jaccard"),
    correlation = c(cor_semantic, cor_op, cor_kappa, cor_dice, cor_jaccard)
  )

  # Print the correlation results
  print(cor_results)

  # Reshape data for plotting
  plot_data <- data.frame(
    sim_wang = rep(df$sim_wang, 5),
    sim_value = c(df$sim_semantic, df$sim_op, df$sim_kappa, df$sim_dice, df$sim_jaccard),
    metric = factor(rep(c("Semantic", "Overlap", "Kappa", "Dice", "Jaccard"), each = nrow(df)))
  )

  # Define colors for each metric
  metric_colors <- c("Semantic" = "#e9c46b", "Overlap" = "#e66f51", "Kappa" = "#264653", "Dice" = "#299e8c", "Jaccard" = "#f2a361")

  # Create the plot with all regression lines in one panel
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sim_wang, y = sim_value, color = metric)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
    ggplot2::scale_color_manual(values = metric_colors) +
    ggplot2::labs(
      title = "Correlation between Wang Similarity and Other Metrics",
      x = "Wang Similarity",
      y = "Similarity Value",
      color = "Similarity Metric"
    ) +
    ggplot2::theme_minimal() +
    # Add annotation with correlation values
    ggplot2::annotate(
      "text",
      x = rep(min(df$sim_wang, na.rm = TRUE) + 0.05, 5),
      y = seq(
        from = max(c(df$sim_semantic, df$sim_op, df$sim_kappa, df$sim_dice, df$sim_jaccard), na.rm = TRUE) - 0.05,
        to = max(c(df$sim_semantic, df$sim_op, df$sim_kappa, df$sim_dice, df$sim_jaccard), na.rm = TRUE) - 0.35,
        length.out = 5
      ),
      label = sprintf(
        "%s: r = %.3f",
        c("Semantic", "Overlap", "Kappa", "Dice", "Jaccard"),
        c(cor_semantic, cor_op, cor_kappa, cor_dice, cor_jaccard)
      ),
      color = metric_colors,
      hjust = 0,
      fontface = "bold"
    )

  # Return the correlation results and plot
  return(list(correlations = cor_results, plot = p))
}

## Check the result ====
library(VennDiagram)
library(grid)

A <- all_gene_list[["hsa03030"]]
B <- all_gene_list[["R-HSA-68952"]]
C <- all_gene_list[["R-HSA-69306"]]

# Calculate overlaps
area1 <- length(A)
area2 <- length(B)
area3 <- length(C)
n12   <- length(intersect(A, B))
n23   <- length(intersect(B, C))
n13   <- length(intersect(A, C))
n123  <- length(Reduce(intersect, list(A, B, C)))

# Draw the triple Venn diagram
venn_plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12   = n12,
  n23   = n23,
  n13   = n13,
  n123  = n123,
  category = c("hsa03030 \nDNA replication", "R-HSA-68952 \nDNA replication initiation", "R-HSA-69306 \nDNA replication"),

  # Specify fill colors and transparency
  fill  = c("#d15356", "#ebb869", "#4f94d7"),
  alpha = c(0.5, 0.5, 0.5),

  # Specify outline (border) colors, line style, and line width
  col   = c("#d15356", "#ebb869", "#4f94d7"),
  lty   = "solid",
  lwd   = 2
)

# Render the diagram
grid.newpage()
grid.draw(venn_plot)

##
diff <- go_combined_sim_df %>% filter(sim_wang < 0.5 | sim_semantic < 0.5) %>% filter(!(sim_wang < 0.5 & sim_semantic < 0.5))
top_dif <- diff %>% mutate(diff = sim_semantic - sim_wang) %>% arrange(desc(diff)) %>% head(5)
top_dif

control_pathways <- control_sim_clustering
id_fm_4 <- control_pathways %>% filter(expected_module =="Functional_module_12") %>% pull(id)

combined_sim_df %>% filter(from %in% id_fm_4 & to %in% id_fm_4)

go_combined_sim_df %>% filter(from == "GO:0044183" | to == "GO:0044183")

id_single <- control_pathways %>% filter(expected_count == 1) %>% pull(id)
single_sim <- combined_sim_df %>% filter(from %in% id_single | to %in% id_single)

b4_single_sim <- b4_all_sim %>% filter(from %in% id_single | to %in% id_single)

single_sim_b4_aftr <- left_join(single_sim, b4_single_sim[,c(1:3)], by = c("from", "to"))

single_long_sim <- single_sim %>%
  pivot_longer(
    cols = c(sim_semantic, sim_jaccard, sim_op, sim_kappa, sim_dice),
    names_to = "sim_metrics",
    values_to = "sim"
  )

single_long_sim_b4_aftr <- single_sim_b4_aftr %>%
  pivot_longer(
    cols = c(sim_semantic, sim_jaccard, sim_op, sim_kappa, sim_dice, b4_mod_sim_semantic),
    names_to = "sim_metrics",
    values_to = "sim"
  )

single_long_sim %>%
  ggplot(aes(x = sim_metrics, y = sim, colour = sim_metrics)) +
  geom_boxplot() +
  geom_point() +
  theme_bw()


pdf("control_similarity_boxplot_single_sim.pdf", width=10, height=8)

single_long_sim_b4_aftr$sim_metrics <- factor(single_long_sim_b4_aftr$sim_metrics, levels = c("sim_jaccard", "sim_kappa", "", "sim_dice", "sim_op", "sim_semantic", "b4_mod_sim_semantic"))

p <-
  single_long_sim_b4_aftr %>%
  ggplot(aes(x = sim_metrics, y = sim, colour = sim_metrics)) +
  geom_boxplot() +
  geom_point() +
  scale_y_continuous(limits = c(0, 0.7)) +
  theme_bw()
print(p)
dev.off()


single_sim_b4_aftr %>%
  filter(sim_semantic >= 0.5) %>%
  arrange(desc(sim_semantic)) %>%
  head(5)


# II. Clustering ====

## 1. Create tidygraph object ====
### Collect edge data
control_pathways <- control_data

semantic_sim_df <-
  as.data.frame.table(sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

### For graph based clustering
sim.cutoff <- 0.5
edge_data <-
  semantic_sim_df %>%
  dplyr::filter(sim > sim.cutoff)

### For clustering (hcluster, binary cut)
edge_data <-
  semantic_sim_df

node_data <-
  control_pathways %>%
  dplyr::filter(id %in% edge_data$from | id %in% edge_data$to) %>%
  dplyr::rename(node = id) %>%
  dplyr::select(node, expected_module, expected_count, database, name)

#### Expected modules
#exp_fm <- levels(factor(node_data$expected_module))

#### Create tidygraph object
graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE,
                       node_key = "node") %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())

## 2. Clustering ====

### 2.1 edge_betweenness clustering ====
subnetwork <-
  suppressWarnings(igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(igraph::edge_attr(graph_data, "sim"))))
#### Assign functional module label for pathways
cluster <-
  paste("Functional_module", as.character(igraph::membership(subnetwork)), sep = "_")

### 2.2 hierarchical clustering ====
cosine_dist <- 1 - sim_matrix
#### Convert distance matrix to a 'dist' object
cosine_dist_obj <- as.dist(cosine_dist)
#### Perform hierarchical clustering
hc <- hclust(cosine_dist_obj, method = "complete")
# #### Plot the hierarchical tree
# plot(hc, main = "Hierarchical Clustering Dendrogram", xlab = "pathways")
# rect.hclust(hc, k = 12, border = "red")  # Highlight 3 clusters
#### Cut the tree to assign clusters (number of clusters(k) or height(h))
clusters <- cutree(hc, k = 9)
# clusters <- cutree(hc, h = 0.7)

#### Assign functional module label for pathways
cluster <-
  paste("Functional_module", as.character(clusters[node_data$node]), sep = "_")

### 2.3 Markov chain clustering ====
res <- MCL::mcl(sim_matrix,
           addLoops = FALSE,
           max.iter = 500,
           expansion = 2,
           inflation = 2.5,
           allow1 = FALSE)
res$Cluster

### 2.4 Spectral clustering ====
clusters <- Spectrum::Spectrum(sim_matrix, maxk = 100, showres = FALSE, silent = TRUE, clusteralg = 'km')
clusters <- clusters$assignments

### 2.5 Binary cut clustering ====
clusters <- simplifyEnrichment::binary_cut(mat = sim_matrix, cutoff = 0.7)
levels(as.factor(clusters))
simplifyEnrichment::plot_binary_cut(mat = sim_matrix, cutoff = 0.7)

cluster <-
  paste("Functional_module", as.character(clusters), sep = "_")
cluster_label <- data.frame(node = rownames(sim_matrix), module = cluster)

### Update graph data with clustering result
graph_data <-
  graph_data %>%
  igraph::upgrade_graph() %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::mutate(module = cluster)

g_node_data <-
  graph_data %>%
  tidygraph::activate(what = "nodes") %>%
  tidygraph::as_tibble()

#### For clustering (hcluster)
graph_data <- graph_data %>%
  activate(edges) %>%
  filter(
    # Get the module attribute for the from node
    .N()$module[from] ==
      # Get the module attribute for the to node
      .N()$module[to]
)

#### For binary cut clustering
graph_data <-
  graph_data %>%
  igraph::upgrade_graph() %>%
  tidygraph::activate(what = "nodes") %>%
  dplyr::left_join(cluster_label)

graph_data <- graph_data %>%
  activate(edges) %>%
  filter(
    # Get the module attribute for the from node
    .N()$module[from] ==
      # Get the module attribute for the to node
      .N()$module[to]
  )

## 3. Sankeyplot to show the relation between the expected functional module and actual module by clustering ====
# p_list = c("ggalluvial", "ggplot2")
# for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
#   library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
#
# library(devtools)
# if(!requireNamespace("amplicon", quietly = TRUE))
#   install_github("davidsjoberg/ggsankey")

suppressWarnings(suppressMessages(library(ggalluvial)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggsankey)))

# Data format transformation
data <- g_node_data %>% dplyr::select(expected_module, module)
data$expected_module <- paste0("Expected_", data$expected_module)
data$expected_module <- tolower(data$expected_module)
data$expected_module <- stringr::str_to_title(data$expected_module)

df <- to_lodes_form(data[, 1:ncol(data)], axes = 1:ncol(data), id = "value")

colors <- colorRampPalette(c('#0ca9ce', '#78cfe5', '#c6ecf1', '#ff6f81', '#ff9c8f', '#ffc2c0','#d386bf',
                             '#cdb1d2', '#fae6f0', '#eb6fa6', '#ff88b5', '#00b1a5',"#ffa68f","#ffca75","#97bc83","#acd295",
                             "#00ada1","#009f93","#ace2da","#448c99","#00b3bc","#b8d8c9","#db888e","#e397a4","#ead0c7",
                             "#8f9898","#bfcfcb"))(12)


pdf("control_cluster_binary_cut.pdf", width=10, height=8)
p1 <- ggplot(df, aes(x = x, stratum = stratum, alluvium = value, fill = stratum, label = stratum)) +
  geom_flow(width = 0.15, curve_type = "cubic", alpha = 0.7, color = 'gray80', size = 0.2) +
  geom_stratum(width = 0.15, color = "gray90") +
  geom_text(stat = 'stratum', size = 3, color = 'gray30') +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(5, 5, 5, 5))
print(p1)

dev.off()

## 4. Create network ====
plot1 <-
  graph_data %>%
  ggraph::ggraph(layout = 'fr',
                 circular = FALSE) +
  ggforce::geom_mark_ellipse(
    aes(x = x,
        y = y,
        group = expected_module,
        color = expected_module),
    alpha = 1,
    expand = unit(5, "mm"),
    linewidth = 1,
    label.fontsize = 9,
    con.cap = 0,
    fill = NA,
    con.type = "straight",
    show.legend = TRUE
  ) +
  ggraph::geom_edge_link(
    aes(width = sim),
    color = "black",
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    aes(fill = database),
    size = 6,
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = c(GO = "#eeca40", KEGG = "#fd7541", Reactome = "#23b9c7")) +
  guides(
    color = guide_legend(order = 1, ncol = 1, title = "Expected_functional_module", override.aes = list(linewidth = 0.5)),
    fill = guide_legend(order = 2, ncol = 1, title = "Database", override.aes = list(linewidth = 0)),
    edge_width = guide_legend(order = 3, title = "Similarity", override.aes = list(linewidth = 0))) +
  ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
  scale_size_continuous(range = c(3, 10)) +
  labs(fill = "Database") +
  ggraph::theme_graph() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "left",
    legend.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggraph::geom_node_text(
    aes(x = x,
        y = y,
        label = stringr::str_wrap(paste0(module, ":", name), width = 30)),
    check_overlap = TRUE,
    size = 3,
    repel = TRUE)

library(Cairo)
CairoPDF("control_cluster_binary_cut_network.pdf", width=18, height=18)
p1 <- plot1
print(p1)
dev.off()

## Evaluation of clustering ====
library(clValid)
node_data <-
  graph_data %>%
  tidygraph::activate(what = "nodes") %>%
  tidygraph::as_tibble()
cluster_label <-
  node_data %>%
  dplyr::select(node, module) %>%
  dplyr::mutate(module_id = as.numeric(sub(pattern = "Functional_module_", replacement = "", x = module)))

# Give single module a cluster label
module_node <- node_data$node %>% unique()
single_node <- setdiff(rownames(sim_matrix), module_node)

max_id <- max(cluster_label$module_id)  # Current maximum module_id
new_ids <- seq(max_id + 1, max_id + length(single_node))  # Increment for each new node

new_cluster_label <- data.frame(
  node = single_node,
  module = paste0("Functional_module_", seq_along(single_node)),  # Placeholder module names
  module_id = new_ids
)

cluster_label_with_single <- rbind(cluster_label, new_cluster_label)

## reorder label

cosine_dist <- 1 - openai_semantic_sim_matrix
cosine_dist_obj <- as.dist(cosine_dist)

dunn_index <- clValid::dunn(distance = cosine_dist_obj, clusters = clusters)
silhouette_index <- mean(cluster::silhouette(clusters, cosine_dist_obj)[, "sil_width"])
