library(rio)
library(dplyr)
library(stringr)

metadata <- import("2_data/metadata/CRA015001_snrna_seq.xlsx", sheet = 3)

r1_link <- metadata %>%
  dplyr::pull(`DownLoad Read file1`) %>%
  stringr::str_replace("ftp://download\\.big\\.ac\\.cn/gsa/",
              "ftp://download.big.ac.cn/gsa2/")


r2_link <- metadata %>%
  dplyr::pull(`DownLoad  Read file2`) %>%
  stringr::str_replace("ftp://download\\.big\\.ac\\.cn/gsa/",
              "ftp://download.big.ac.cn/gsa2/")


all_link <- c(r1_link, r2_link)

writeLines(all_link, "2_data/metadata/links_CRA015128_multi-tissue_metformin.txt")

# Extract links for Proteomics data downloading links from xml file
library(xml2)
library(dplyr)

all_data <- xml2::read_xml("2_data/metadata/PX_IPX0008294000.xml")

kids <- xml_children(all_data)
sub_kids <- xml_children(kids[12])

http_links <- sub_kids %>%
  purrr::map(function(x) {
    ssub_kids <- xml_children(x)
    xml_attr(ssub_kids, "value")
  }) %>% unlist()

# convert a vector into a list of elements in txt
writeLines(http_links, "2_data/metadata/proteomics_dl_links.txt")


