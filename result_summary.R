##此脚本用于分别对d1-d5群体中每个文件计算cor和bias的均值和标准差
library(tidyverse)
data_dir <- "/public/home/liujf/workspace/shikp/TempWork/20250612_piggtex/runs/d5/"
files <- list.files(
  path = data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)
# 汇总结果
result_summary <- tibble(file_path = files) %>%
  mutate(
    file_name = basename(file_path),
    file_name_noext = tools::file_path_sans_ext(file_name),
    file_name_noext_clean = str_remove(file_name_noext, "_ee$"),
    # 假设文件名格式固定为 trait_tissue.csv
    trait = str_extract(file_name_noext_clean, "^[^_]+"),
    tissue = str_extract(file_name_noext_clean, "(?<=_).+$")
  ) %>%
  rowwise() %>%
  mutate(
    dat = list(read_csv(file_path, show_col_types = FALSE)),
    cor_mean = mean(dat$cor, na.rm = TRUE),
    cor_sd = sd(dat$cor, na.rm = TRUE),
    bias_mean = mean(dat$bias, na.rm = TRUE),
    bias_sd = sd(dat$bias, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(trait, tissue, cor_mean, cor_sd, bias_mean, bias_sd, file_path)
write_csv(result_summary, file = file.path(data_dir, "summary_result.csv"))
