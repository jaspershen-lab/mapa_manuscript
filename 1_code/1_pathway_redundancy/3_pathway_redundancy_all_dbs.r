library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)

#dir.create("3_data_analysis/1_pathway_redundancy/all_dbs")
# dir.create("2_data/pathway_redundancy/")
# dir.create("2_data/pathway_redundancy/within_db_sim")
# dir.create("2_data/pathway_redundancy/cross_db")

# Within database similarity ====
## Get GO terms and similarity =====
# "Wang" similarity
load("3_data_analysis/1_pathway_redundancy/2_go/sem_matrix_bp.rda")
load("3_data_analysis/1_pathway_redundancy/2_go/sem_matrix_mf.rda")
load("3_data_analysis/1_pathway_redundancy/2_go/sem_matrix_cc.rda")

## Get KEGG pathways and sim =====
library(KEGGREST)
library(pbapply)

### Get KEGG pathway id =====
organism <- "hsa"

# Get all pathway IDs for the organism
pathway_ID <- KEGGREST::keggList(database = "pathway", organism = organism) |> names()

grepl("hsa[0-9]{5}", pathway_ID) |> sum()

### Get KEGG pathway similarity =====
kegg_sim_res <- term_similarity_KEGG(term_id = pathway_ID,
                                     measure.method = "overlap")
save(kegg_sim_res, file = "2_data/pathway_redundancy/kegg_sim_res.rda")

term_similarity_KEGG <- function(term_id,
                                 measure.method = c("jaccard", "dice", "overlap", "kappa")) {

  measure.method <- match.arg(measure.method)

  species <- gsub("^([a-zA-Z]+)(\\d+$)", "\\1", term_id[1])

  kegg_data <- tryCatch(
    expr = {
      #A function that downloads KEGG data for a given species, KEGG type, and key type
      getFromNamespace("prepare_KEGG", "clusterProfiler")(species, "KEGG", "kegg")
    },
    error = function(e) {
      getFromNamespace("get_data_from_KEGG_db", "clusterProfiler")(species)
    })

  gl <- kegg_data$PATHID2EXTID[term_id]

  term_similarity_internal(gl = gl,
                           measure.method = measure.method)
}

## Get Reactome ====
### Get Reactome pathway =====
## Download pathway information: https://download.reactome.org/94/ReactomePathways.txt
reactome_all_pathways <- readr::read_delim(file = "2_data/pathway_redundancy/ReactomePathways.txt",
                                           delim = "\t",
                                           col_names = c("id", "term_name", "organism"))

human_reactome_pathway_id <- reactome_all_pathways |> dplyr::filter(organism == "Homo sapiens")
### Calculate Reactome pathway similarity =====
reactome_sim_res <- term_similarity_Reactome(term_id = human_reactome_pathway_id$id,
                                             measure.method = "overlap")

save(reactome_sim_res, file = "2_data/pathway_redundancy/reactome_sim_res.rda")

term_similarity_Reactome <- function(term_id,
                                     measure.method = c("jaccard", "dice", "overlap", "kappa")) {

  measure.method <- match.arg(measure.method)

  all <- as.list(reactome.db::reactomePATHID2EXTID)
  available_terms <- term_id[term_id %in% names(all)]
  gl <- all[available_terms]

  term_similarity_internal(gl = gl,
                           measure.method = measure.method)
}

term_similarity_internal <-
  function(gl,
           measure.method = c("jaccard", "dice", "overlap", "kappa"),
           remove_negative = TRUE) {

    measure.method <- match.arg(measure.method)

    all <- unique(unlist(gl))
    gl <- lapply(gl, function(x) as.numeric(factor(x, levels = all)))
    n <- length(gl)

    pathway_gene_m <- matrix(0, ncol = length(all), nrow = n)
    for(i in seq_len(n)) {
      pathway_gene_m[i, gl[[i]]] = 1
    }
    pathway_gene_m <- as(pathway_gene_m, "sparseMatrix")

    if(measure.method == "kappa") {
      mat <- kappa_dist(pathway_gene_m, remove_negative = remove_negative)
    } else if(measure.method == "overlap") {
      mat <- overlap_dist(pathway_gene_m)
    } else {
      mat <- proxyC::simil(pathway_gene_m, method = measure.method)
    }

    sim_matrix <-  as.matrix(mat)
    diag(sim_matrix) <- 1
    rownames(sim_matrix) = colnames(sim_matrix) = names(gl)

    return(sim_matrix)
  }

kappa_dist <- function(m, remove_negative = TRUE) {
  tab <- ncol(m)
  po <- proxyC::simil(m, method = "simple matching")
  m_yes <- Matrix::rowSums(m)
  m_no <- abs(Matrix::rowSums(m - 1))
  pe <- (outer(m_yes, m_yes, FUN = "*") + outer(m_no, m_no, FUN = "*"))/tab^2
  k <- (po - pe)/(1 - pe)
  if(remove_negative) k[k < 0] <- 0
  return(k)
}

overlap_dist <- function(m) {
  n = Matrix::rowSums(m)
  proxyC::simil(m, method = "dice")*outer(n, n, FUN = "+")/2/outer(n, n, pmin)
}

# Cross database similarity ====
## Get GO term gene list ====
go_data <- as.data.frame(org.Hs.egGO)

# Filter only "BP", "MF", or "CC"
bp_data <-
  go_data |>
  dplyr::filter(Ontology == "BP")

mf_data <-
  go_data |>
  dplyr::filter(Ontology == "MF")

cc_data <-
  go_data |>
  dplyr::filter(Ontology == "CC")


# Helper to build list: GO ID -> list of Entrez gene IDs
make_go_list <- function(df) {
  split(df$gene_id, df$go_id)
}

go2gene_bp <- make_go_list(bp_data)
go2gene_mf <- make_go_list(mf_data)
go2gene_cc <- make_go_list(cc_data)

go2gene_bp <- go2gene_bp[sapply(go2gene_bp, length) >= 10]
go2gene_mf <- go2gene_mf[sapply(go2gene_mf, length) >= 10]
go2gene_cc <- go2gene_cc[sapply(go2gene_cc, length) >= 10]

## Get KEGG pathways gene list ====
species <- "hsa"
kegg_data <- tryCatch(
  expr = {
    #A function that downloads KEGG data for a given species, KEGG type, and key type
    getFromNamespace("prepare_KEGG", "clusterProfiler")(species, "KEGG", "kegg")
  },
  error = function(e) {
    getFromNamespace("get_data_from_KEGG_db", "clusterProfiler")(species)
  })

kegg_gl <- kegg_data$PATHID2EXTID

## Get Reactome pathways gene list ====
reactome_gl <- as.list(reactome.db::reactomePATHID2EXTID)
reactome_all_pathways <- readr::read_delim(file = "2_data/pathway_redundancy/ReactomePathways.txt",
                                           delim = "\t",
                                           col_names = c("id", "term_name", "organism"))

human_reactome_pathway_id <- reactome_all_pathways |> dplyr::filter(organism == "Homo sapiens")
human_ava_ids <- human_reactome_pathway_id$id[human_reactome_pathway_id$id %in% names(reactome_gl)]
reactome_gl <- reactome_gl[human_ava_ids]

## Calculate BP-KEGG ====
gobp_kegg_gl <- c(kegg_gl, go2gene_bp)
gobp_kegg <- term_similarity_internal(gl = gobp_kegg_gl,
                                      measure.method = "overlap")
save(gobp_kegg, file = "2_data/pathway_redundancy/gobp_kegg.rda")

## Calculate BP-Reactome ====
gobp_reactome_gl <- c(reactome_gl, go2gene_bp)
gobp_reactome <- term_similarity_internal(gl = gobp_reactome_gl,
                                          measure.method = "overlap")

save(gobp_reactome, file = "2_data/pathway_redundancy/gobp_reactome.rda")

## Calculate KEGG-Reactome ====
kegg_reactome_gl <- c(kegg_gl, reactome_gl)
kegg_reactome <- term_similarity_internal(gl = kegg_reactome_gl,
                                          measure.method = "overlap")
save(kegg_reactome, file = "2_data/pathway_redundancy/kegg_reactome.rda")

# Plot data ====
load("3_data_analysis/1_pathway_redundancy/2_go/sem_matrix_bp.rda")
load("3_data_analysis/1_pathway_redundancy/2_go/sem_matrix_mf.rda")
load("3_data_analysis/1_pathway_redundancy/2_go/sem_matrix_cc.rda")

within_bp_bp_sim <- tibble(
  term1 = rownames(sem_matrix_bp)[row(sem_matrix_bp)],
  term2 = colnames(sem_matrix_bp)[col(sem_matrix_bp)],
  sim = as.vector(sem_matrix_bp)
)

within_bp_bp_sim <- within_bp_bp_sim[as.character(within_bp_bp_sim$term1) < as.character(within_bp_bp_sim$term2), ]

within_mf_mf_sim <- tibble(
  term1 = rownames(sem_matrix_mf)[row(sem_matrix_mf)],
  term2 = colnames(sem_matrix_mf)[col(sem_matrix_mf)],
  sim = as.vector(sem_matrix_mf)
)

within_mf_mf_sim <- within_mf_mf_sim[as.character(within_mf_mf_sim$term1) < as.character(within_mf_mf_sim$term2), ]

within_cc_cc_sim <- tibble(
  term1 = rownames(sem_matrix_cc)[row(sem_matrix_cc)],
  term2 = colnames(sem_matrix_cc)[col(sem_matrix_cc)],
  sim = as.vector(sem_matrix_cc)
)

within_cc_cc_sim <- within_cc_cc_sim[as.character(within_cc_cc_sim$term1) < as.character(within_cc_cc_sim$term2), ]

load("2_data/pathway_redundancy/kegg_sim_res.rda")
within_kegg_sim <- tibble(
  term1 = rownames(kegg_sim_res)[row(kegg_sim_res)],
  term2 = colnames(kegg_sim_res)[col(kegg_sim_res)],
  sim = as.vector(kegg_sim_res)
)

within_kegg_sim <- within_kegg_sim[as.character(within_kegg_sim$term1) < as.character(within_kegg_sim$term2), ]

load("2_data/pathway_redundancy/reactome_sim_res.rda")
within_reactome_sim <- tibble(
  term1 = rownames(reactome_sim_res)[row(reactome_sim_res)],
  term2 = colnames(reactome_sim_res)[col(reactome_sim_res)],
  sim = as.vector(reactome_sim_res)
)

within_reactome_sim <- within_reactome_sim[as.character(within_reactome_sim$term1) < as.character(within_reactome_sim$term2), ]

save(within_bp_bp_sim, file = "2_data/pathway_redundancy/within_db_sim/within_bp_bp_sim.rda")
save(within_mf_mf_sim, file = "2_data/pathway_redundancy/within_db_sim/within_mf_mf_sim.rda")
save(within_cc_cc_sim, file = "2_data/pathway_redundancy/within_db_sim/within_cc_cc_sim.rda")
save(within_kegg_sim, file = "2_data/pathway_redundancy/within_db_sim/within_kegg_sim.rda")
save(within_reactome_sim, file = "2_data/pathway_redundancy/within_db_sim/within_reactome_sim.rda")

## Cross database sim ====
load("2_data/pathway_redundancy/gobp_kegg.rda")
load("2_data/pathway_redundancy/gobp_reactome.rda")
load("2_data/pathway_redundancy/kegg_reactome.rda")

cross_gobp_kegg <-tibble(
  term1 = rownames(gobp_kegg)[row(gobp_kegg)],
  term2 = colnames(gobp_kegg)[col(gobp_kegg)],
  sim = as.vector(gobp_kegg)
)
cross_gobp_kegg <- cross_gobp_kegg[as.character(cross_gobp_kegg$term1) < as.character(cross_gobp_kegg$term2), ]
cross_gobp_kegg <- cross_gobp_kegg |>
  dplyr::filter(substr(term1, 1, 3) != substr(term2, 1, 3))

cross_gobp_reactome <-tibble(
  term1 = rownames(gobp_reactome)[row(gobp_reactome)],
  term2 = colnames(gobp_reactome)[col(gobp_reactome)],
  sim = as.vector(gobp_reactome)
)
cross_gobp_reactome <- cross_gobp_reactome[as.character(cross_gobp_reactome$term1) < as.character(cross_gobp_reactome$term2), ]
cross_gobp_reactome <- cross_gobp_reactome |>
  dplyr::filter(substr(term1, 1, 3) != substr(term2, 1, 3))

cross_kegg_reactome <-tibble(
  term1 = rownames(kegg_reactome)[row(kegg_reactome)],
  term2 = colnames(kegg_reactome)[col(kegg_reactome)],
  sim = as.vector(kegg_reactome)
)
cross_kegg_reactome <- cross_kegg_reactome[as.character(cross_kegg_reactome$term1) < as.character(cross_kegg_reactome$term2), ]
cross_kegg_reactome <- cross_kegg_reactome |>
  dplyr::filter(substr(term1, 1, 3) != substr(term2, 1, 3))

save(cross_gobp_kegg, file = "2_data/pathway_redundancy/cross_db/cross_gobp_kegg.rda")
save(cross_gobp_reactome, file = "2_data/pathway_redundancy/cross_db/cross_gobp_reactome.rda")
save(cross_kegg_reactome, file = "2_data/pathway_redundancy/cross_db/cross_kegg_reactome.rda")

## Put all together
load("2_data/pathway_redundancy/within_db_sim/within_bp_bp_sim.rda")
load("2_data/pathway_redundancy/within_db_sim/within_mf_mf_sim.rda")
load("2_data/pathway_redundancy/within_db_sim/within_cc_cc_sim.rda")
load("2_data/pathway_redundancy/within_db_sim/within_kegg_sim.rda")
load("2_data/pathway_redundancy/within_db_sim/within_reactome_sim.rda")
load("2_data/pathway_redundancy/cross_db/cross_gobp_kegg.rda")
load("2_data/pathway_redundancy/cross_db/cross_gobp_reactome.rda")
load("2_data/pathway_redundancy/cross_db/cross_kegg_reactome.rda")

within_bp_bp_sim <- within_bp_bp_sim |> dplyr::mutate(group = "within_database", pair_type = "BP_BP")
within_mf_mf_sim <- within_mf_mf_sim |> dplyr::mutate(group = "within_database", pair_type = "MF_MF")
within_cc_cc_sim <- within_cc_cc_sim |> dplyr::mutate(group = "within_database", pair_type = "CC_CC")
within_kegg_sim <- within_kegg_sim |> dplyr::mutate(group = "within_database", pair_type = "KEGG_KEGG")
within_reactome_sim <- within_reactome_sim |> dplyr::mutate(group = "within_database", pair_type = "Reactome_Reactome")

cross_gobp_kegg <- cross_gobp_kegg |> dplyr::mutate(group = "cross_database", pair_type = "BP_KEGG")
cross_gobp_reactome <- cross_gobp_reactome |> dplyr::mutate(group = "cross_database", pair_type = "BP_Reactome")
cross_kegg_reactome <- cross_kegg_reactome |> dplyr::mutate(group = "cross_database", pair_type = "KEGG_Reactome")

all_pairs <- rbind(within_bp_bp_sim,
                   within_mf_mf_sim,
                   within_cc_cc_sim,
                   within_kegg_sim,
                   within_reactome_sim,
                   cross_gobp_kegg,
                   cross_gobp_reactome,
                   cross_kegg_reactome)

save(all_pairs, file = "2_data/pathway_redundancy/all_pairs.rda")
