setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

dt <- import("2_data/case_study/mmc2.xlsx",
             skip = 1,
             header = TRUE,
             sheet = 2)
dt <- dt |> rename(ensembl = GeneID)

# cluster U
dt_u <- dt |> filter(Cluster == "Cluster U")
dt_u <- dt_u |> distinct(ensembl, .keep_all = TRUE)

# cluster D
dt_d <- dt |> filter(Cluster == "Cluster D")
dt_d <- dt_d |> distinct(ensembl, .keep_all = TRUE)

# cluster UD
dt_ud <- dt |> filter(Cluster == "Cluster UD")
dt_ud <- dt_ud |> distinct(ensembl, .keep_all = TRUE)

# cluster DU
dt_du <- dt |> filter(Cluster == "Cluster DU")
dt_du <- dt_du |> distinct(ensembl, .keep_all = TRUE)

save(dt_u, file = "2_data/case_study/bulk_rna_seq/dt_u.rda")
save(dt_d, file = "2_data/case_study/bulk_rna_seq/dt_d.rda")
save(dt_ud, file = "2_data/case_study/bulk_rna_seq/dt_ud.rda")
save(dt_du, file = "2_data/case_study/bulk_rna_seq/dt_du.rda")
