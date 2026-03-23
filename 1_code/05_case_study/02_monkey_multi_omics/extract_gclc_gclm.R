library(r4projects)
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

# mrna
file_path <- "2_data/case_study/02_monkey_multi_omics/raw_data/41592_2025_2912_MOESM4_ESM.xlsx"

sheet_names <- excel_sheets(file_path)

# Loop through each sheet, filter for GCLM or GCLC, and add tissue_name column
result <- lapply(sheet_names, function(sheet) {
  df <- read_excel(file_path, sheet = sheet)
  df %>%
    # filter(ID %in% c("GCLM", "GCLC")) %>%
    filter(ID %in% c("SRXN1", "PRDX2")) %>%
    mutate(tissue_name = sheet)
}) %>%
  bind_rows()

T_res_sig <- result |> dplyr::filter(Pvalue < 0.05)

# prot
file_path <- "2_data/case_study/02_monkey_multi_omics/raw_data/41592_2025_2912_MOESM5_ESM.xlsx"

sheet_names <- excel_sheets(file_path)

# Loop through each sheet, filter for GCLM or GCLC, and add tissue_name column
result <- lapply(sheet_names, function(sheet) {
  df <- read_excel(file_path, sheet = sheet)
  df %>%
    # filter(ID %in% c("GCLM", "GCLC")) %>%
    # filter(ID %in% c("SRXN1", "PRDX2")) %>%
    filter(ID %in% c("LAMTOR1", "LAMTOR2", "RRAGA", "MIOS")) %>%
    mutate(tissue_name = sheet)
}) %>%
  bind_rows()

P_res_sig <- result |> dplyr::filter(Pvalue < 0.05)

all_res <- bind_rows(T_res_sig |> mutate(src = "T"), P_res_sig |> mutate(src = "P"))

rio::export(all_res, file = "2_data/case_study/02_monkey_multi_omics/taxid_9606/spleen/up/srxn1_prdx2_up_res.xlsx")
