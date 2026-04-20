library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

setwd("2_data/case_study/02_monkey_multi_omics/")

plot_similarity_network <-
  function(object,
           module_ids,
           node_size = 8,
           llm_text = TRUE,
           text = TRUE,
           text_all = FALSE) {

    sim_method <- object@process_info$merge_pathways@function_name

    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }


    # Determine analysis type
    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
      query_type <- object@process_info$enrich_pathway@parameter$query_type
    } else {
      analysis_type <- "do_gsea"
      query_type <- "gene"
    }

    graph_data <-
      object@merged_module$graph_data|>
      tidygraph::activate(what = "nodes")

    result_with_module <-
      object@merged_module$result_with_module

    graph_data <-
      graph_data |>
      tidygraph::activate(what = "nodes") |>
      dplyr::filter(module %in% module_ids) |>
      # dplyr::mutate(label = if (query_type == "gene") Description else pathway_name)
      dplyr::mutate(label = node)

    if (igraph::gorder(graph_data) == 0) {
      warning("No functional modules have degree > ", degree_cutoff)
      return(ggplot() +
               geom_blank())
    }

    ###plot to show the clusters of GO terms
    cluster_label_module <-
      tryCatch(
        expr = {
          df <- igraph::as_data_frame(graph_data, what = "vertices")

          if (analysis_type == "enrich_pathway") {
            if ("Count" %in% colnames(df)) {
              df |>
                dplyr::group_by(module) |>
                dplyr::arrange(p_adjust, desc(Count), .by_group = TRUE) |>
                dplyr::slice_head(n = 1) |>
                dplyr::pull(label)
            } else {
              df |>
                dplyr::group_by(module) |>
                dplyr::arrange(p_adjust, desc(mapped_number), .by_group = TRUE) |>
                dplyr::slice_head(n = 1) |>
                dplyr::pull(label)
            }
          } else {
            df |>
              dplyr::group_by(module) |>
              dplyr::arrange(desc(abs(NES)), desc(Count), .by_group = TRUE) |>
              dplyr::slice_head(n = 1) |>
              dplyr::pull(label)
          }
        },
        error = function(e) {

        },
        warning = function(w) {
          # Handle warning here if needed
        }
      )

    if (is.null(cluster_label_module)) {
      cluster_label_module <- ""
    }

    cluster_label_all <-
      graph_data |>
      tidygraph::activate(what = "nodes") |>
      dplyr::pull(label)

    lay <- ggraph::create_layout(graph_data, layout = "fr")

    # colors <- colorRampPalette(c('#0ca9ce', '#78cfe5', '#c6ecf1', '#ff6f81', '#ff9c8f', '#ffc2c0','#d386bf',
    #                              '#cdb1d2', '#fae6f0', '#eb6fa6', '#ff88b5', '#00b1a5',"#ffa68f","#ffca75","#97bc83","#acd295",
    #                              "#00ada1","#009f93","#ace2da","#448c99","#00b3bc","#b8d8c9","#db888e","#e397a4","#ead0c7",
    #                              "#8f9898","#bfcfcb"))(length(unique(lay$module)))

    plot <-
      ggraph::ggraph(lay) +
      ggraph::geom_edge_link(
        # aes(width = sim),
        width = 0.5,
        color = "grey",
        alpha = 1,
        show.legend = FALSE
      ) +
      ggraph::geom_node_point(
        size = node_size,
        shape = 22,
        fill = "#bcd59b",
        alpha = 1,
        show.legend = FALSE
      ) +
      # guides(fill = guide_legend(ncol = 1)) +
      # scale_fill_manual(values = colors) +
      # ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
      # scale_size_continuous(range = c(1, 7)) +
      # labs(size = if(analysis_type == "enrich_pathway") "-log10(FDR adjusted P-values)" else "abs(NES)") +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )

    if (text_all) {
      plot <-
        plot +
        ggraph::geom_node_text(aes(x = x,
                                   y = y,
                                   label = label),
                               size = 3,
                               repel = TRUE)
    } else if (llm_text && level == "functional_module") {
      node_tbl  <- as_tibble(lay)
      centroids <- node_tbl |>
        group_by(module) |>
        summarise(cx = mean(x), cy = mean(y), .groups = "drop")
      module_name_df <- object@merged_module$functional_module_result |> dplyr::select(module, llm_module_name)
      centroids <- centroids |> dplyr::left_join(module_name_df, by = "module")

      plot <-
        plot +
        ggraph::geom_node_text(
          data = centroids,
          aes(x = cx,
              y = cy,
              label = stringr::str_wrap(llm_module_name, 30)),
          check_overlap = TRUE,
          size = 3,
          repel = TRUE
        )
    } else if (text) {
      plot <-
        plot +
        ggraph::geom_node_text(aes(
          x = x,
          y = y,
          label = ifelse(label %in% cluster_label_module, label, NA)
        ),
        size = 3,
        repel = TRUE)
    }

    plot
  }

# Spleen up ====
load("taxid_9606/spleen/up/llm_t_res.rda")

# Single module
## mrna
set.seed(14)
plot <-
plot_similarity_network(
  object = llm_t_res,
  node_size = 10,
  # module_id = c("Functional_module_2", "Functional_module_6", "Functional_module_10",
  #               "Functional_module_15", "Functional_module_16", "Functional_module_17"),
  module_ids = "Functional_module_2",
  llm_text = TRUE,
  text_all = TRUE
)

plot

library(Cairo)
CairoPDF("taxid_9606/spleen/up/spleen_up_t.pdf", width = 5, height = 4)
plot
dev.off()

## prot
load("taxid_9606/spleen/up/llm_p_res.rda")

set.seed(1)

plot <-
  plot_similarity_network(
    object = llm_p_res,
    node_size = 10,
    module_ids = c("Functional_module_6", "Functional_module_10", "Functional_module_21"),
    llm_text = TRUE,
    text_all = TRUE
  )

plot

library(Cairo)
CairoPDF("taxid_9606/spleen/up/spleen_up_p.pdf", width = 5, height = 4)
plot
dev.off()

## metabolite
load("taxid_9606/spleen/met_enriched_pathways.rda")

# aortic arch ====
## mrna
load("taxid_9606/aortic_arch/up/llm_t_res.rda")

set.seed(11)

plot <-
  plot_similarity_network(
    object = llm_t_res,
    node_size = 10,
    module_id = c("Functional_module_60"),
    llm_text = TRUE,
    text_all = TRUE
  )

plot

library(Cairo)
CairoPDF("taxid_9606/aortic_arch/up/aortic_arch_t.pdf", width = 5, height = 4)
plot
dev.off()

## prot
load("taxid_9606/aortic_arch/up/llm_p_res.rda")

set.seed(15)

plot <-
  plot_similarity_network(
    object = llm_p_res,
    node_size = 10,
    module_ids = c("Functional_module_68"),
    llm_text = TRUE,
    text_all = TRUE
  )

plot

library(Cairo)
CairoPDF("taxid_9606/aortic_arch/up/aortic_arch_up_p.pdf", width = 5, height = 4)
plot
dev.off()

## metabolite
load("taxid_9606/spleen/met_enriched_pathways.rda")



