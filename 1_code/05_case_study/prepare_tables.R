library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

# Step 1: variables conversion

# Step 2: Pathway enrichment analysis
# Input: Pathway enrichment results from mapa
# Transcriptome
load("3_data_analysis/05_case_study/01_monkey/brain_sn_rna_seq/brain_up_enrich_pathway_res.rda")
T_input <- enrich_pathway_res
T_input@variable_info
# Metabolome
load("3_data_analysis/05_case_study/01_monkey/plasma_metabolomics/functional_module_res_u_and_d_met.rda")
M_input <- functional_module_res
# Proteome
load("3_data_analysis/05_case_study/01_monkey/plasma_proteomics/02_proteomics_up_enrich_pathway_res.rda")
P_input <- enrich_pathway_res

## Get pathway layer information
get_TP_enrich_res <- function(enrich_res) {
  if (!is.null(enrich_res@enrichment_go_result)) {
    go_res <- enrich_res@enrichment_go_result@result
    go_res <- go_res |>
      dplyr::select(ID, Description, p_adjust, geneID) |>
      dplyr::filter(p_adjust < 0.05)
  } else {
    go_res <- tibble::tibble(ID = NA_character_, Description = NA_character_,
                             p_adjust = NA, geneID = NA_character_)
  }
  if (!is.null(enrich_res@enrichment_kegg_result)) {
    kegg_res <- enrich_res@enrichment_kegg_result@result
    kegg_res <- kegg_res |>
      dplyr::select(ID, Description, p_adjust, geneID) |>
      dplyr::filter(p_adjust < 0.05)
  } else {
    kegg_res <- tibble::tibble(ID = NA_character_, Description = NA_character_,
                               p_adjust = NA, geneID = NA_character_)
  }
  if (!is.null(enrich_res@enrichment_reactome_result)) {
    reactome_res <- enrich_res@enrichment_reactome_result@result
    reactome_res <- reactome_res |>
      dplyr::select(ID, Description, p_adjust, geneID) |>
      dplyr::filter(p_adjust < 0.05)
  } else {
    reactome_res <- tibble::tibble(ID = NA_character_, Description = NA_character_,
                                   p_adjust = NA, geneID = NA_character_)
  }

  all_enrich_res <- dplyr::bind_rows(go_res, kegg_res, reactome_res) |>
    dplyr::filter(!is.na(ID)) |>
    dplyr::rename(pathway_id = ID, pathway_name = Description, mapped_id = geneID)

  all_enrich_res
}
T_enrich_res <- get_TP_enrich_res(enrich_res = T_input)
P_enrich_res <- get_TP_enrich_res(enrich_res = P_input)

enrich_res <- M_input
get_M_enrich_res <- function(enrich_res) {
  if (!is.null(enrich_res@enrichment_hmdb_result)) {
    hmdb_res <- enrich_res@enrichment_hmdb_result@result
    hmdb_res <- hmdb_res |>
      dplyr::select(pathway_id, pathway_name, p_adjust, mapped_id) |>
      dplyr::filter(p_adjust < 0.05)
  } else {
    hmdb_res <- tibble::tibble(pathway_id = NA_character_, pathway_name = NA_character_,
                               p_adjust = NA, mapped_id = NA_character_)
  }

  if (!is.null(enrich_res@enrichment_metkegg_result)) {
    metkegg_res <- enrich_res@enrichment_metkegg_result@result
    metkegg_res <- metkegg_res |>
      dplyr::select(pathway_id, pathway_name, p_adjust, mapped_id) |>
      dplyr::filter(p_adjust < 0.05)
  } else {
    metkegg_res <- tibble::tibble(pathway_id = NA_character_, pathway_name = NA_character_,
                                  p_adjust = NA, mapped_id = NA_character_)
  }

  all_enrich_res <- dplyr::bind_rows(hmdb_res, metkegg_res) |>
    dplyr::filter(!is.na(pathway_id))

  all_enrich_res
}
M_enrich_res <- get_M_enrich_res(enrich_res = M_input)

all_enrich_res <- dplyr::bind_rows(T_enrich_res, P_enrich_res, M_enrich_res)

pathway_molecule_pairs <- all_enrich_res |>
  tidyr::separate_rows(mapped_id, sep = "\\s*[/;]\\s*") |>
  dplyr::filter(!is.na(mapped_id), mapped_id != "") |>
  dplyr::distinct(pathway_id, mapped_id, .keep_all = TRUE)

# Step 3: Prepare table for multi-omics integration ====
# Step 3.1: Get multi-omics DE table =====
T_DE <- T_input@variable_info
P_DE <- P_input@variable_info
M_DE <- M_input@variable_info

# Step 4:Get edges ====
## Step 4.1: Get Gene-Gene (TF-target) edges ====
data("dorothea_hs", package = "dorothea")  # human TF-target prior knowledge

T_deg_symbol <- T_DE$symbol

edges_in_list <- dorothea_hs |>
  dplyr::filter(
    tf %in% T_deg_symbol,
    target %in% T_deg_symbol,
    confidence %in% c("A")
  ) |>
  dplyr::distinct(tf, target, mor, confidence) |>
  dplyr::arrange(confidence, tf, target)

edges_in_list

edges_in_list$confidence |> table()

## Step 4.2: Get Protein-Protein (PPI) edges ====
taxid <- 9541

string_db <- STRINGdb$new(
  version = "12.0",
  species = taxid,
  score_threshold = 900,
  input_directory = ""
)

P_deg_symbol <- P_DE[, "symbol", drop = FALSE]

# Map to STRING IDs
mapped <- string_db$map(P_deg_symbol, "symbol", removeUnmappedRows = TRUE)

# Interactions *among* the proteins in your list :contentReference[oaicite:3]{index=3}
ppi <- string_db$get_interactions(mapped$STRING_id)

# https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
# https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
# https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz

## Step 4.3: Get Protein-Metabolite edges ====
### From Reactome database ====
chebl2rn_url <- "https://download.reactome.org/95/ChEBI2Reactome_PE_Reactions.txt"

# https://download.reactome.org/96/ChEBI2Reactome_PE_Reactions.txt
# dir.create("2_data/case_study/01_monkey/reactome_db")
dest <- file.path("2_data/case_study/01_monkey/reactome_db", "ChEBI2Reactome_PE_Reactions.txt")
system2("curl", c("-L", "--retry", "3", "--fail", "-o", shQuote(dest), shQuote(chebl2rn_url)))

uniprot2rn_url <- "https://download.reactome.org/95/UniProt2Reactome_PE_Reactions.txt"
# dir.create("2_data/case_study/01_monkey/reactome_db")
dest <- file.path("2_data/case_study/01_monkey/reactome_db", "UniProt2Reactome_PE_Reactions.txt")
system2("curl", c("-L", "--retry", "3", "--fail", "-o", shQuote(dest), shQuote(uniprot2rn_url)))
# https://download.reactome.org/96/UniProt2Reactome_PE_Reactions.txt

prot2rn_url <- "https://download.reactome.org/95/ProteinRoleReaction.txt"
# dir.create("2_data/case_study/01_monkey/reactome_db")
dest <- file.path("2_data/case_study/01_monkey/reactome_db", "ProteinRoleReaction.txt")
system2("curl", c("-L", "--retry", "3", "--fail", "-o", shQuote(dest), shQuote(prot2rn_url)))
# https://download.reactome.org/96/ProteinRoleReaction.txt

sp <- "HSA"
rxn_pat <- paste0("^R-", sp, "-")

# protein role in reaction
prot_role <- readr::read_delim(
  "2_data/case_study/01_monkey/reactome_db/ProteinRoleReaction.txt",
  col_names = c("protein_id", "role", "reaction_id"),
  show_col_types = FALSE
) |>
  dplyr::filter(str_detect(reaction_id, rxn_pat))

enz_rxn <- prot_role |>
  dplyr::filter(role == "catalystActivity") |>
  dplyr::distinct(protein_id, reaction_id) |>
  dplyr::rename(uniprotkb = protein_id)

# taxid <- 9541 # Macaca fasciculari
P_deg_symbol <- P_DE$symbol

ensembl <- biomaRt::useMart("ensembl")
datasets <- biomaRt::listDatasets(ensembl)
head(datasets)
dataset_name <- datasets$dataset[grepl("hsapiens_", datasets$dataset)]
mart <- biomaRt::useMart("ensembl", dataset = dataset_name)

symbol_attr <- if ("hgnc_symbol" %in% biomaRt::listAttributes(mart)$name) "hgnc_symbol" else "external_gene_name"

attrs_sym2uni <- c(symbol_attr, "uniprotswissprot", "uniprotsptrembl")

sym2uni_raw <- biomaRt::getBM(
  attributes = attrs_sym2uni,
  filters    = symbol_attr,
  values     = P_deg_symbol,
  mart       = mart
)

sym2uni <- sym2uni_raw |>
  dplyr::mutate(
    symbol   = .data[[symbol_attr]],
    uniprotkb = dplyr::coalesce(uniprotswissprot, uniprotsptrembl)
  ) |>
  dplyr::filter(!is.na(uniprotkb) & uniprotkb != "") |>
  dplyr::select(symbol, uniprotkb) |>
  dplyr::distinct()

enz_rxn <- sym2uni |>
  dplyr::left_join(enz_rxn, by = "uniprotkb") |>
  dplyr::filter(!is.na(reaction_id))

# ChEBI2Reactome_Reactions
chebi2rn <- readr::read_delim(
  "2_data/case_study/01_monkey/reactome_db/ChEBI2Reactome_PE_Reactions.txt",
  col_names = c("ChEBI", "reactomeid", "reactome_name", "reaction_id",
                "URL", "reaction_name", "evidence_Code", "species"),
  show_col_types = FALSE) |>
  dplyr::select(ChEBI, reactomeid, reactome_name, reaction_id, reaction_name) |>
  dplyr::filter(stringr::str_detect(reaction_id, rxn_pat)) |>
  dplyr::distinct(ChEBI, reaction_id, .keep_all = TRUE)

chunk_vec <- function(x, chunk_size = 100) {
  split(x, ceiling(seq_along(x) / chunk_size))
}

kegg_to_chebi <- function(kegg_ids,
                          chunk_size = 100,
                          sleep_sec = 0.1,
                          verbose = TRUE) {
  kegg_unique <- unique(kegg_ids)
  chunks <- chunk_vec(kegg_unique, chunk_size = chunk_size)

  all_maps <- list()
  for (i in seq_along(chunks)) {
    q <- paste0("compound:", chunks[[i]])
    if (verbose) message(sprintf("Query chunk %d/%d (n=%d)", i, length(chunks), length(q)))

    res <- tryCatch(
      KEGGREST::keggConv("chebi", q, querySize = chunk_size),
      error = function(e) {
        warning(sprintf("Chunk %d failed: %s", i, conditionMessage(e)))
        character(0)
      }
    )

    if (length(res)) all_maps[[length(all_maps) + 1]] <- res
    if (sleep_sec > 0) Sys.sleep(sleep_sec)
  }

  m <- if (length(all_maps)) unlist(all_maps, use.names = TRUE) else character(0)

  out <- data.frame(
    keggid  = sub("^cpd:", "", names(m)),
    ChEBI = as.numeric(sub("^chebi:", "", unname(m))),
    stringsAsFactors = FALSE
  )

  missing <- setdiff(kegg_unique, out$keggid)
  if (length(missing)) {
    out <- rbind(out, data.frame(keggid = missing, ChEBI = NA_integer_))
  }
  out <- out[order(out$keggid), , drop = FALSE]

  out
}

M_de_chebi <- kegg_to_chebi(kegg_ids = M_DE$keggid)

chebi2rn <- M_de_chebi |>
  dplyr::left_join(chebi2rn, by = "ChEBI") |>
  dplyr::filter(!is.na(reaction_id))

reactome_edges_P2M <- dplyr::inner_join(enz_rxn, chebi2rn,
                                        by = "reaction_id",
                                        relationship = "many-to-many") |>
  dplyr::distinct(uniprotkb, ChEBI, reaction_id, .keep_all = TRUE) |>
  dplyr::select(reaction_id, reaction_name, symbol, uniprotkb, keggid, ChEBI, reactomeid, reactome_name)

### From KEGG database ====
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
chunk_vec <- function(x, size = 10) split(x, ceiling(seq_along(x) / size))

get_kegg_enzyme_compound_reactions <- function(kegg_ids, symbols, organism = "hsa", sleep_sec = 0.2) {
  # 1) Map SYMBOL -> Entrez -> KEGG gene IDs
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(symbols),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ) |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::distinct(SYMBOL, ENTREZID) |>
    dplyr::rename(symbol = SYMBOL, entrezid = ENTREZID) |>
    dplyr::mutate(kegg_gene = paste0(organism, ":", entrezid))

  if (nrow(gene_map) == 0) {
    return(NA)
  }

  # 2) Reactions linked to compounds
  cpd_kegg <- paste0("cpd:", unique(kegg_ids))
  Sys.sleep(sleep_sec)
  cpd2rn <- KEGGREST::keggLink("reaction", cpd_kegg)

  cpd2rn_tbl <- tibble::tibble(
    keggid = sub("^cpd:", "", names(cpd2rn)),
    reaction_id = sub("^rn:", "", unname(cpd2rn))
  ) |> dplyr::distinct()

  # 3) Reactions linked to genes (enzymes)
  Sys.sleep(sleep_sec)
  # gene -> EC
  g2ec <- KEGGREST::keggLink("enzyme", unique(gene_map$kegg_gene))

  # EC -> reaction
  ec_ids <- unname(g2ec)
  ec2rn  <- if (length(ec_ids)) KEGGREST::keggLink("reaction", ec_ids) else character(0)

  # tidy mapping (gene -> reaction via ec)
  g2ec_tbl <- tibble::tibble(kegg_gene = names(g2ec), ec = unname(g2ec))
  ec2rn_tbl <- tibble::tibble(ec = names(ec2rn), reaction_id = sub("^rn:", "", unname(ec2rn)))

  gene_to_reaction <- dplyr::left_join(g2ec_tbl, ec2rn_tbl, by = "ec", relationship = "many-to-many") |>
    dplyr::filter(!is.na(reaction_id))

  gene2rn_tbl <- gene_to_reaction |>
    dplyr::left_join(gene_map, by = "kegg_gene")

  # Intersect on reaction => enzyme–compound via same reaction
  edges <- dplyr::inner_join(cpd2rn_tbl, gene2rn_tbl, by = "reaction_id") |>
    dplyr::distinct(symbol, entrezid, ec, kegg_gene, keggid, reaction_id) |>
    dplyr::select(reaction_id, symbol, entrezid, ec, kegg_gene, keggid)

  if (nrow(edges) == 0) return(edges)

  # 5) add equation for the reaction
  rn_ids <- unique(edges$reaction_id)
  rxn_list <- purrr::map(chunk_vec(paste0("rn:", rn_ids), size = 10), function(chunk) {
    Sys.sleep(sleep_sec)
    KEGGREST::keggGet(chunk)
  }) |>  purrr::flatten()

  rxn_tbl <- purrr::map_dfr(rxn_list, function(x) {
    tibble::tibble(
      reaction_id  = unname(x$ENTRY) %||% NA_character_,
      reaction_name = x$NAME %||% NA_character_,
      definition = x$DEFINITION %||% NA_character_,
      equation  = x$EQUATION %||% NA_character_
    )
  }) |> dplyr::distinct(reaction_id, .keep_all = TRUE)

  edges <- edges |> left_join(rxn_tbl, by = "reaction_id")
}

kegg_edges_P2M <- get_kegg_enzyme_compound_reactions(kegg_ids = M_DE$keggid,
                                                     symbols = P_DE$symbol)

### Combine KEGG and Reactome results =====
all_cols <- c(
  "source",
  "reaction_id", "reaction_info",
  "symbol", "entrezid", "uniprotkb", "ec", "kegg_gene",
  "keggid", "ChEBI", "reactomeid", "reactome_name"
)

kegg_edges_P2M_norm <- kegg_edges_P2M |>
  tibble::as_tibble() |>
  dplyr::mutate(
    source = "KEGG",
    # Reactome-only columns -> NA
    uniprotkb    = NA_character_,
    ChEBI        = NA_character_,
    reactomeid   = NA_character_,
    reactome_name= NA_character_,
    reaction_info = paste(reaction_name, definition, equation, sep = "\n")
  ) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
  dplyr::select(dplyr::any_of(all_cols))

reactome_edges_P2M_norm <- reactome_edges_P2M |>
  tibble::as_tibble() |>
  dplyr::mutate(
    source = "Reactome",
    # KEGG-only columns -> NA
    entrezid   = NA_character_,
    ec         = NA_character_,
    kegg_gene  = NA_character_
  ) |>
  dplyr::rename(reaction_info = reaction_name) |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
  dplyr::select(dplyr::any_of(all_cols))

edges_P2M_merged <- dplyr::bind_rows(kegg_edges_P2M_norm, reactome_edges_P2M_norm)



