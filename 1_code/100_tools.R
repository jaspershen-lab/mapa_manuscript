library(mapa)
library(tidyverse)
library(ggplot2)
library(rio)
library(data.table)
library(AnnotationHub)

database_color =
  c(
    GO = "#1F77B4FF",
    KEGG = "#FF7F0EFF",
    Reactome = "#2CA02CFF"
  )

