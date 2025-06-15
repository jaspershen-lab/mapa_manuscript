library(mapa)
library(tidyverse)
library(ggplot2)
library(rio)
library(data.table)
library(AnnotationHub)

database_color <-
  c(
    GO = "#eeca40",
    KEGG = "#fd7541",
    Reactome = "#23b9c7"
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
