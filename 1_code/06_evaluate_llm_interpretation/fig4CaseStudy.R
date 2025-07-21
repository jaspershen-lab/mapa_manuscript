library(dplyr)
library(readxl)
library(org.Hs.eg.db)
library(reactome.db)
library(AnnotationDbi)
library(KEGGREST)
library(mapa)
load("/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa doc/fig4CaseStudy/ora_enriched_functional_module.rda")
control_data <- read_excel("/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa doc/fig4CaseStudy/control_data.xlsx", sheet = 1)
module_counts <- table(control_data$expected_module)
modules_to_keep <- names(module_counts[module_counts > 1])
control_data_filtered <- control_data[control_data$expected_module %in% modules_to_keep, ]

map_pathway_to_genes <- function(pid) {
  ## 1) GO 通路
  if (grepl("^GO", pid, perl = TRUE)) {
    go_genes <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = pid,
      columns  = "ENSEMBL", # 将 SYMBOL 改为 ENSEMBL
      keytype  = "GO"
    )
    return(paste(unique(go_genes$ENSEMBL), collapse = "/"))
  }

  ## 2) Reactome 通路（ID 形如 R-HSA-69306）
  if (grepl("^R-", pid, perl = TRUE)) {
    eid <- unlist(as.list(reactomePATHID2EXTID[pid]))
    gene_df <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = eid,
      columns  = "ENSEMBL", # 将 SYMBOL 改为 ENSEMBL
      keytype  = "ENTREZID"
    )
    return(paste(unique(gene_df$ENSEMBL), collapse = "/"))
  }

  ## 3) KEGG 通路（ID 形如 hsa00190）
  if (grepl("^hsa", pid, perl = TRUE)) {
    eid_raw <- keggLink("hsa", paste0("path:", pid))
    eid     <- sub("hsa:", "", unname(eid_raw))
    gene_df <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = eid,
      columns  = "ENSEMBL", # 将 SYMBOL 改为 ENSEMBL
      keytype  = "ENTREZID"
    )
    return(paste(unique(gene_df$ENSEMBL), collapse = "/"))
  }

  ## 4) 其它情况：返回 NA
  return(NA_character_)
}

# ── 2. 批量应用到数据框 ────────────────────────────────────────
control_data_filtered <- control_data_filtered %>%
  mutate(
    gene_symbols = vapply(id, map_pathway_to_genes, FUN.VALUE = character(1))
  )


randomize_modules <- function(df, module_col_name = "expected_module") {
  # 检查指定的列是否存在
  if (!module_col_name %in% names(df)) {
    stop(paste("数据框中找不到列:", module_col_name))
  }

  # 1. 统计每个模块的原始数量
  module_counts <- df %>%
    count(!!sym(module_col_name), name = "count")

  # 2. 创建新的带有“Random”前缀的模块名称
  new_modules <- paste0("Random ", module_counts[[module_col_name]])

  # 3. 根据原始数量重新构建一个包含所有新模块的向量
  # 这个向量将用于随机抽样
  all_new_modules_for_sampling <- rep(new_modules, times = module_counts$count)

  # 4. 随机分配这些新模块
  # sample函数默认是不放回抽样，这里我们要保证总数不变
  df[[module_col_name]] <- sample(all_new_modules_for_sampling)

  return(df)
}

control_data_randomized = randomize_modules(control_data_filtered,)
combined_data <- rbind(control_data_filtered, control_data_randomized)

formated_result <- aggregate(cbind(id, name) ~ expected_module,
                     data = combined_data,
                     FUN = function(x) paste(x, collapse = ";"))
formated_result$geneID <- aggregate(gene_symbols ~ expected_module,
                                    data = combined_data,
                                    FUN = function(x) paste(x, collapse = "/"))$gene_symbols

module_content_counts <- combined_data %>%
  group_by(expected_module) %>%
  summarise(module_content_number = n())


formated_result <- left_join(formated_result, module_content_counts, by = "expected_module")
colnames(formated_result)[1:3] = c("module","pathway_id","Description")

#加上一列$module_content_number
openai_key = "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"
enriched_functional_module@merged_module[["functional_module_result"]] = formated_result
# processed_module = llm_interpret_module(enriched_functional_module, api_key=openai_key,embedding_output_dir = "./embeddings")
object <- enriched_functional_module
api_key <- openai_key
embedding_output_dir <- "./embeddings" 

llm_model = "gpt-4o-mini-2024-07-18"
embedding_model = "text-embedding-3-small"
phenotype = NULL
local_corpus_dir = NULL
chunk_size = 5
years = 5
retmax = 10
similarity_filter_num = 20
GPT_filter_num = 5
orgdb = org.Hs.eg.db
output_prompt = TRUE

# -- 流程开始 --

# 步骤 1: 初始检查和数据提取
# ---------------------------------
message("Step 1: Running initial checks and extracting data...")
if (!is(object, "functional_module")) {
  stop("Error: 'object' should be a functional_module class object.")
}
if (all(names(object@process_info) != "merge_modules")) {
  stop("Error: Please use the merge_modules() function to process the object first.")
}

# 从 object 中提取功能模块结果
functional_module_result <- object@merged_module$functional_module_result
message(" -> Intermediate variable 'functional_module_result' is created.")

# 步骤 2: 设置本地语料库参数
# ---------------------------------
message("Step 2: Setting up local corpus parameters...")
local_corpus <- FALSE
save_dir_local_corpus_embed <- NULL # 初始化变量

if (!is.null(local_corpus_dir)) {
  local_corpus <- TRUE
  save_dir_local_corpus_embed <- "local"
  message(" -> Local corpus will be used.")
} else {
  message(" -> Local corpus not provided, skipping.")
}
clear_output_dir <- function(output_dir = NULL) {
  # Check if directory exists
  if (file.exists(output_dir)) {
    # Get all files and subdirectories in the directory
    files <- list.files(output_dir, full.names = TRUE)
    
    if (length(files) != 0) {
      # Recursively delete files and subdirectories
      unlink(files, recursive = TRUE)
      cat("Successfully cleared the directory:", output_dir, "\n")
    }
    
    # # Delete the empty directory
    # unlink(output_dir, recursive = TRUE)
  } else {
    stop("Directory does not exist:", output_dir, "\n")
  }
}
# 步骤 3: 清理输出目录并处理本地语料库 (如果提供)
# ---------------------------------
message("Step 3: Preparing output directory and processing local corpus...")
clear_output_dir(output_dir = embedding_output_dir)

if (!is.null(local_corpus_dir)) {
  message(" -> Generating embeddings for the local corpus...")
  embedding_local_corpus(embedding_model = embedding_model, 
                         api_key = api_key, 
                         local_corpus_dir = local_corpus_dir, 
                         embedding_output_dir = embedding_output_dir, 
                         save_dir_local_corpus_embed = save_dir_local_corpus_embed)
}

# 步骤 4: 预处理模块数据
# ---------------------------------
message("Step 4: Preprocessing module data...")
source('/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa/R/16_llm_module_pathway_infor.R')
source('/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa/R/16_llm_interpret_module.R')
source('/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa/R/16_llm_module_get_pathway_name_desc_pmid.R')
source('/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa/R/16_llm_module_embedding_database.R')
source('/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa/R/16_llm_module_online_retrieval.R')
source('/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa/R/16_llm_module_output_generation.R')
source('/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa/R/16_llm_module_RAG_strategy.R')
source('/Users/cgxjdzz/Desktop/NTU phd/KEGG RAG mapa/mapa/R/16_llm_module_utils.R')
processed_data <- preprocess_module(df = functional_module_result, orgdb = orgdb)
message(" -> Intermediate variable 'processed_data' is created.")

# 步骤 5: 从 PubMed 搜索相关文献
# ---------------------------------
message("Step 5: Searching PubMed for relevant literature...")
pubmed_result <- pubmed_search(processed_data = processed_data, 
                               chunk_size = chunk_size, 
                               years = years, 
                               retmax = 2)
message(" -> Intermediate variable 'pubmed_result' is created.")

# 步骤 6: 合并和整理 PubMed ID
# ---------------------------------
message("Step 6: Consolidating PubMed IDs...")
for (module in names(pubmed_result)) {
  pmids <- pubmed_result[[module]][["PathwayReferencePMID"]]
  if (!is.null(pmids)) {
    pmids <- pmids[pmids != "" & !is.na(pmids)]
    if (length(pmids) > 0) {
      if (!is.null(pubmed_result[[module]][["PubmedIDs"]])) {
        pubmed_result[[module]][["PubmedIDs"]] <- unique(c(pubmed_result[[module]][["PubmedIDs"]], pmids))
      } else {
        pubmed_result[[module]][["PubmedIDs"]] <- pmids
      }
    }
  }
}
message(" -> 'pubmed_result' has been updated with consolidated PMIDs.")

# 步骤 7: 为 PubMed 搜索结果生成嵌入向量
# ---------------------------------
message("Step 7: Generating embeddings for PubMed search results...")
library(pbmcapply)
library(data.table)
embedding_pubmed_search(pubmed_result = pubmed_result, 
                        embedding_model = embedding_model, 
                        api_key = api_key, 
                        embedding_output_dir = embedding_output_dir)
  



library(httr2)
# 步骤 8: 检索和筛选相关文献
# ---------------------------------
message("Step 8: Retrieving and filtering relevant papers...")
related_paper <- retrieve_strategy(pubmed_result = pubmed_result, 
                                   model = llm_model, 
                                   embedding_model = embedding_model, 
                                   api_key = api_key, 
                                   similarity_filter_num = similarity_filter_num, 
                                   GPT_filter_num = GPT_filter_num, 
                                   local_corpus = local_corpus, 
                                   embedding_output_dir = embedding_output_dir, 
                                   save_dir_local_corpus_embed = save_dir_local_corpus_embed)

message(" -> Intermediate variable 'related_paper' is created.")

# 步骤 9: 组合文献结果
# ---------------------------------
message("Step 9: Combining paper results...")
paper_result <- Map(function(x, y) {
  return(list(related_paper = x, pubmed_result = y))
}, related_paper, pubmed_result)
message(" -> Intermediate variable 'paper_result' is created.")

# 步骤 10: 生成模块的功能性名称
# ---------------------------------
message("Step 10: Generating functional names for modules using LLM...")
final_result <- module_name_generation(paper_result = paper_result, 
                                       phenotype = phenotype, 
                                       model = llm_model, 
                                       api_key = api_key, 
                                       output_prompt = output_prompt)
message(" -> Intermediate variable 'final_result' is created.")

# # 假设你的数据存储在 final_result 中
# module_names <- names(final_result)
# true_modules <- module_names[!grepl("^Random ", module_names)]  # 找到真实模块名
# 
# # 初始化空的data.frame用于收集结果
# library(dplyr)
# all_scores <- data.frame()
# 
# # 遍历每个真实模块
# for (mod in true_modules) {
#   random_mod <- paste0("Random ", mod)
#   
#   # 提取真实与随机模块的 score
#   score_true <- final_result[[mod]]$generated_name$confidence_score
#   score_rand <- final_result[[random_mod]]$generated_name$confidence_score
#   
#   # 转换为data.frame并添加标签
#   df_true <- data.frame(score = score_true, group = "True", module = mod)
#   df_rand <- data.frame(score = score_rand, group = "Random", module = mod)
#   
#   # 合并
#   all_scores <- bind_rows(all_scores, df_true, df_rand)
# }
# # # 拆出两个向量
# # true_scores <- as.numeric(as.character(true_scores))
# # random_scores <- as.numeric(as.character(random_scores))
# # 
# # # 做 t 检验
# # t.test(true_scores, random_scores)
# all_scores$score <- as.numeric(all_scores$score)
# 
# library(ggplot2)
# library(ggpubr)
# 
# ggplot(all_scores, aes(x = group, y = score, fill = group)) +
#   geom_boxplot( width = 0.4, alpha = 0.5) +  # 箱型图
#   geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.6) +  # 点图
#   stat_compare_means(method = "t.test", label = "p.format") +  # 显著性
#   theme_minimal(base_size = 14) +
#   scale_fill_manual(values = c("True" = "#1b9e77", "Random" = "#d95f02")) +
#   labs(title = "Confidence Score: True vs Random Modules",
#        x = "", y = "Confidence Score")
# 
module_names <- names(final_result)
true_modules <- module_names[!grepl("^Random ", module_names)]  # 找到真实模块名

# 初始化空的data.frame用于收集结果
library(dplyr)
library(ggplot2)
library(ggpubr)

all_scores <- data.frame()

# 遍历每个真实模块
for (mod in true_modules) {
  random_mod <- paste0("Random ", mod)
  
  # 提取真实与随机模块的 score
  score_true <- final_result[[mod]]$generated_name$confidence_score
  score_rand <- final_result[[random_mod]]$generated_name$confidence_score
  
  # 转换为data.frame并添加标签
  df_true <- data.frame(score = score_true, group = "True", module = mod)
  df_rand <- data.frame(score = score_rand, group = "Random", module = mod)
  
  # 合并
  all_scores <- bind_rows(all_scores, df_true, df_rand)
}

# 确保score是数值型
all_scores$score <- as.numeric(all_scores$score)

# 创建一个固定的jitter位置，这样连线和点用的是同一个jitter
set.seed(123)  # 设置随机种子保证结果可重复
jitter_width <- 0.1

# 确保group的顺序一致
all_scores$group <- factor(all_scores$group, levels = c("Random", "True"))

# 为每个点计算jitter后的坐标（横向和纵向都加jitter）
all_scores$x_jittered <- ifelse(all_scores$group == "Random", 
                                1 + runif(nrow(all_scores), -jitter_width, jitter_width),
                                2 + runif(nrow(all_scores), -jitter_width, jitter_width))

# 添加纵向jitter（但幅度要小，避免影响数据解读）
y_jitter_width <- 0.005  # 很小的纵向抖动
all_scores$y_jittered <- all_scores$score + runif(nrow(all_scores), -y_jitter_width, y_jitter_width)

# 绘图
p = ggplot(all_scores, aes(x = group, y = score)) +
  geom_boxplot(aes(fill = group), width = 0.4, alpha = 0.5) +  # 箱型图
  
  # 添加连接线 - 使用jitter后的坐标
  geom_line(aes(x = x_jittered, y = y_jittered, group = module), 
            color = "gray50", 
            alpha = 0.7, 
            size = 0.8) +
  
  # 添加点图 - 使用相同的jitter坐标
  geom_point(aes(x = x_jittered, y = y_jittered, color = group), 
             size = 3, 
             alpha = 0.8) +
  
  stat_compare_means(method = "t.test", label = "p.format") +  # 显著性
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Random" = "#d95f02", "True" = "#1b9e77")) +
  scale_color_manual(values = c("Random" = "#d95f02", "True" = "#1b9e77")) +
  labs(title = "Confidence Score: True vs Random Modules",
       x = "", y = "Confidence Score")
p
ggsave("Jun30_confidence_score_comparison.pdf", plot = p, 
       width = 10, height = 8, units = "in", dpi = 300, device = "pdf")
save.image("Jun30_fig4CaseStudy_workspace.RData")

