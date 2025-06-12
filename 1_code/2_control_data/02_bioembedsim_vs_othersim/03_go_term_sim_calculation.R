library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
# source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

library(GOSemSim)
library(org.Hs.eg.db)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

## BP
go_bp_id <- control_dt %>% dplyr::filter(database == "GO" & go_ontology == "BP") %>% pull(id)
bp_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "BP")
bp_sim_matrix <- termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata, method = "Wang")
bp_sim_df <-
  as.data.frame.table(bp_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

## MF
go_mf_id <- control_dt %>% dplyr::filter(database == "GO" & go_ontology == "MF") %>% pull(id)
mf_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "MF")
mf_sim_matrix <- termSim(t1 = go_mf_id, t2 = go_mf_id, semData = mf_semgodata, method = "Wang")
mf_sim_df <-
  as.data.frame.table(mf_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

## CC
go_cc_id <- control_dt %>% dplyr::filter(database == "GO" & go_ontology == "CC") %>% pull(id)
cc_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "CC")
cc_sim_matrix <- termSim(t1 = go_cc_id, t2 = go_cc_id, semData = cc_semgodata, method = "Wang")
cc_sim_df <-
  as.data.frame.table(cc_sim_matrix, responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

wang_sim_df <- rbind(bp_sim_df, mf_sim_df, cc_sim_df)
save(wang_sim_df, file = "3_data_analysis/02_control_data/02_bioembedsim_vs_othersim/wang_sim_df.rda")
