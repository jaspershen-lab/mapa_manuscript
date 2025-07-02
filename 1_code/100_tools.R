library(mapa)
library(tidyverse)
library(ggplot2)
library(rio)
library(data.table)
library(AnnotationHub)
library(ggthemes)
library(paletteer)

database_color <-
  c(
    GO = "#eeca40",
    KEGG = "#fd7541",
    Reactome = "#23b9c7",
    SMPDB = "#7998ad"
  )

intra_inter_database_color <-
  c(
    intra_database = "#d8d113",
    inter_database = "#2373cd"
  )


same_different_module_color <-
  c(
    same_module = "#ff0000",
    different_module = "#65684c"
  )

# clustering algorithms color
class_colors <- c("Distance_based" = "#1BB6AFFF", "Graph_based" = "#FFAD0AFF")


real_random_module_color <-
  c(
    Real = "#e53ba4",
    Random = "#084d68"
  )
