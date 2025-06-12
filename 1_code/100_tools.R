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
