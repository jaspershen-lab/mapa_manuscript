library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("2_data/Jun30_final_result.RData")
load("3_data_analysis/4_random_module_dataset/formated_result.Rdata")

names(final_result)

dir.create("3_data_analysis/4_random_module_dataset/confidence_score_comparison")
setwd("3_data_analysis/4_random_module_dataset/confidence_score_comparison")

temp_data <-
  final_result %>%
  purrr::map(function(x) {
    x$generated_name$confidence_score
  }) %>%
  unlist() %>%
  as.numeric()

all_scores <-
  data.frame(
    module = names(final_result),
    module_name = names(final_result),
    score = temp_data
  ) %>%
  dplyr::mutate(group = case_when(stringr::str_detect(module, "Random") ~ "Random", TRUE ~ "Real")) %>%
  dplyr::mutate(module = stringr::str_replace(module, "Random ", "")) %>%
  dplyr::mutate(module = stringr::str_replace(module, "Functional_module_", "")) %>%
  dplyr::mutate(module =  factor(module, levels = stringr::str_sort(unique(module), numeric = TRUE)))

module_size <-
  data.frame(
    module_name = formated_result$module,
    size = stringr::str_split(formated_result$pathway_id, ";") %>%
      lapply(length) %>%
      unlist()
  )

all_scores <-
  all_scores %>%
  dplyr::left_join(module_size, by = "module_name") %>%
  dplyr::mutate(module = factor(module, levels = stringr::str_sort(unique(module), numeric = TRUE)))


# 确保score是数值型
all_scores$score <- as.numeric(all_scores$score)

# 创建一个固定的jitter位置，这样连线和点用的是同一个jitter
set.seed(123)  # 设置随机种子保证结果可重复
jitter_width <- 0.1

# 确保group的顺序一致
all_scores$group <- factor(all_scores$group, levels = c("Random", "Real"))

# 为每个点计算jitter后的坐标（横向和纵向都加jitter）
all_scores$x_jittered <- ifelse(
  all_scores$group == "Random",
  1 + runif(nrow(all_scores), -jitter_width, jitter_width),
  2 + runif(nrow(all_scores), -jitter_width, jitter_width)
)

# 添加纵向jitter（但幅度要小，避免影响数据解读）
y_jitter_width <- 0.005  # 很小的纵向抖动
all_scores$y_jittered <- all_scores$score + runif(nrow(all_scores), -y_jitter_width, y_jitter_width)

# 绘图
library(ggpubr)
library(ggside)
plot  <-
  ggplot(all_scores, aes(x = group, y = score)) +
  geom_boxplot(
    aes(fill = group),
    width = 0.4,
    alpha = 0.5,
    outlier.shape = NA
  ) +  # 箱型图
  # 添加连接线 - 使用jitter后的坐标
  geom_line(
    aes(x = x_jittered, y = y_jittered, group = module),
    color = "black",
    alpha = 1
  ) +
  # 添加点图 - 使用相同的jitter坐标
  geom_point(
    aes(
      x = x_jittered,
      y = y_jittered,
      fill = group,
      size = size
    ),
    shape = 21,
    alpha = 0.8
  ) +
  geom_text(
    aes(
      x = x_jittered,
      y = y_jittered,
      label = module_name
    ),
    size = 3,
    vjust = -1.5,
    hjust = 0.5,
    color = "black"
  ) +
  stat_compare_means(method = "t.test", label = "p.signif") +  # 显著性
  theme_bw() +
  scale_fill_manual(values = real_random_module_color) +
  scale_color_manual(values = real_random_module_color) +
  labs(x = "", y = "Confidence Score") +
  scale_size_continuous(range = c(3, 8)) +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))

plot

ggsave(
  plot = plot,
  filename = "confidence_score_comparison.pdf",
  height = 6,
  width = 8
)
