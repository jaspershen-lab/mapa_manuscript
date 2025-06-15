library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
# source('1_code/2_control_data/02_bioembedsim_vs_othersim/utils.R')

library(GOSemSim)
library(org.Hs.eg.db)

control_dt <- readxl::read_excel("2_data/control_data.xlsx", sheet = 1)

setwd("3_data_analysis/02_control_data/02_bioembedsim_vs_othersim")

## BP
go_bp_id <- control_dt %>% dplyr::filter(database == "GO" & go_ontology == "BP") %>% pull(id)
bp_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "BP")

##Wang
bp_sim_df_wang <-
  termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata, method = "Wang") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Resnik
bp_sim_df_resnik <-
  bp_sim_df_wang <-
  termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata, method = "Resnik") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Rel
bp_sim_df_rel <-
  bp_sim_df_wang <-
  termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata, method = "Rel") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Jiang
bp_sim_df_jiang <-
  termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata, method = "Jiang") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Lin
bp_sim_df_lin <-
  termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata, method = "Lin") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

# ##TCSS
# bp_semgodata_tcss <- godata(annoDb = org.Hs.eg.db, ont = "BP", processTCSS = TRUE)
# bp_sim_df_tcss <-
#   termSim(t1 = go_bp_id, t2 = go_bp_id, semData = bp_semgodata_tcss, method = "TCSS") %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)

## MF
go_mf_id <- control_dt %>% dplyr::filter(database == "GO" & go_ontology == "MF") %>% pull(id)
mf_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "MF")

##wang
mf_sim_df_wang <-
  termSim(t1 = go_bp_id, t2 = go_bp_id, semData = mf_semgodata, method = "Wang") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Resnik
mf_sim_df_resnik <-
  mf_sim_df_wang <-
  termSim(t1 = go_mf_id, t2 = go_mf_id, semData = mf_semgodata, method = "Resnik") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Rel
mf_sim_df_rel <-
  mf_sim_df_wang <-
  termSim(t1 = go_mf_id, t2 = go_mf_id, semData = mf_semgodata, method = "Rel") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Jiang
mf_sim_df_jiang <-
  termSim(t1 = go_mf_id, t2 = go_mf_id, semData = mf_semgodata, method = "Jiang") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Lin
mf_sim_df_lin <-
  termSim(t1 = go_mf_id, t2 = go_mf_id, semData = mf_semgodata, method = "Lin") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

# ##TCSS
# mf_semgodata_tcss <- godata(annoDb = org.Hs.eg.db, ont = "BP", processTCSS = TRUE)
# mf_sim_df_tcss <-
#   termSim(t1 = go_mf_id, t2 = go_mf_id, semData = mf_semgodata_tcss, method = "TCSS") %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)


## CC
go_cc_id <- control_dt %>% dplyr::filter(database == "GO" & go_ontology == "CC") %>% pull(id)
cc_semgodata <- godata(annoDb = org.Hs.eg.db, ont = "CC")

##wang
cc_sim_df_wang <-
  termSim(t1 = go_bp_id, t2 = go_bp_id, semData = cc_semgodata, method = "Wang") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Resnik
cc_sim_df_resnik <-
  cc_sim_df_wang <-
  termSim(t1 = go_cc_id, t2 = go_cc_id, semData = cc_semgodata, method = "Resnik") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Rel
cc_sim_df_rel <-
  cc_sim_df_wang <-
  termSim(t1 = go_cc_id, t2 = go_cc_id, semData = cc_semgodata, method = "Rel") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Jiang
cc_sim_df_jiang <-
  termSim(t1 = go_cc_id, t2 = go_cc_id, semData = cc_semgodata, method = "Jiang") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

##Lin
cc_sim_df_lin <-
  termSim(t1 = go_cc_id, t2 = go_cc_id, semData = cc_semgodata, method = "Lin") %>%
  as.data.frame.table(responseName = "sim") %>%
  dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
  dplyr::rename(from = Var1, to = Var2) %>%
  dplyr::mutate(across(c(from, to), as.character)) %>%
  dplyr::filter(from < to)

# ##TCSS
# cc_semgodata_tcss <- godata(annoDb = org.Hs.eg.db, ont = "BP", processTCSS = TRUE)
# cc_sim_df_tcss <-
#   termSim(t1 = go_cc_id, t2 = go_cc_id, semData = cc_semgodata_tcss, method = "TCSS") %>%
#   as.data.frame.table(responseName = "sim") %>%
#   dplyr::filter(Var1 != Var2) %>%                 # Remove self-edges
#   dplyr::rename(from = Var1, to = Var2) %>%
#   dplyr::mutate(across(c(from, to), as.character)) %>%
#   dplyr::filter(from < to)

wang_sim_df <- rbind(bp_sim_df_wang, mf_sim_df_wang, cc_sim_df_wang)
resnik_sim_df <- rbind(bp_sim_df_resnik, mf_sim_df_resnik, cc_sim_df_resnik)
rel_sim_df <- rbind(bp_sim_df_rel, mf_sim_df_rel, cc_sim_df_rel)
jiang_sim_df <- rbind(bp_sim_df_jiang, mf_sim_df_jiang, cc_sim_df_jiang)
lin_sim_df <- rbind(bp_sim_df_lin, mf_sim_df_lin, cc_sim_df_lin)

save(wang_sim_df, file = "wang_sim_df.rda")
save(resnik_sim_df, file = "resnik_sim_df.rda")
save(rel_sim_df, file = "rel_sim_df.rda")
save(jiang_sim_df, file = "jiang_sim_df.rda")
save(lin_sim_df, file = "lin_sim_df.rda")
