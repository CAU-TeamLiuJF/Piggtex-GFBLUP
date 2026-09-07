setwd("/public/home/liujf/workspace/shikp/TempWork/20250612_piggtex")
library(tidyverse)
library(data.table)
library(rrBLUP)
library(lme4)
library(lme4qtl)
library(cli)
# receive parameter
arg = as.numeric(commandArgs(t = T)[1])

if(arg == 1 | arg == 2) {
  snp_map <- read_delim("data/CAU50K_v11.1genome.map", col_names = FALSE, show_col_types = FALSE)
} else if(arg == 3) {
  snp_map <- read_delim("data/American_Duroc_pigs_genotypes.bim", col_names = FALSE, show_col_types = FALSE)
} else if(arg == 4) {
  snp_map <- read_delim("data/Canadian_Duroc_pigs_genotypes.bim", col_names = FALSE, show_col_types = FALSE)
} else if(arg == 5) {
  snp_map <- read_delim("data/中芯一号v1plus_to_11.1.map", col_names = FALSE, show_col_types = FALSE)
} else{
  stop("Error parameter !!!")
}

#cat("\nArg : ", arg, " | data : dataset", arg, "\n")
cli_h2("\nArg :  {arg} | data : dataset {arg} \n")

snp_map <- snp_map %>%
  filter(X1 %in% 1:18) %>%
  mutate(X1 = as.character(X1))


fl = fs::dir_ls("data/PigGTEx_v0.significant_lncQTL/")

st <- function(.x) {
  t1 <- read_delim(.x, show_col_types = F) %>%
    select(variant_id) %>%
    distinct(variant_id) %>%
    separate(
      col = variant_id,
      into = c("chr", "pos", "ref", "alt"),
      sep = "_",
      remove = FALSE
    )
  
  #精准匹配
  exact_match <- t1 %>%
    mutate(pos = as.integer(pos)) %>%
    inner_join(snp_map, by = c("chr" = "X1", "pos" = "X4")) %>%
    distinct(X2, .keep_all = TRUE) %>%
    mutate(match_type = "exact")%>%
    select(X2, match_type)

  #范围匹配
  t1_dt <- as.data.table(t1)
  snp_dt <- as.data.table(snp_map)
  t1_dt[, `:=`(chr = as.character(chr),pos = as.integer(pos),pos_start = as.integer(pos) - 50000L,pos_end = as.integer(pos) + 50000L)]
  snp_dt[, `:=`(chr = as.character(X1),snp_pos = as.integer(X4))]
  snp_dt[, snp_start := snp_pos]
  snp_dt[, snp_end   := snp_pos]
  setkey(snp_dt, chr, snp_start, snp_end)
  setkey(t1_dt, chr, pos_start, pos_end)
  range_match <- foverlaps(snp_dt,t1_dt,by.x = c("chr", "snp_start", "snp_end"),by.y = c("chr", "pos_start", "pos_end"),type = "within",nomatch = 0)
  range_match <- range_match[!X2 %in% exact_match$X2][!duplicated(X2)][, match_type := "range_50kb"]
  range_match <- range_match[,c("X2", "match_type")]
  #依据精准匹配和范围匹配结果，对snp_map中所有位点进行标记exact、range_50kb、no_match
  final_result <- snp_map %>%
    rename(chr = X1, snp = X2, pos = X4) %>%
    left_join(
      bind_rows(exact_match, range_match) %>%
        select(snp = X2, match_type),
      by = "snp"
    ) %>%
    mutate(
      match_type = coalesce(match_type, "no_match")
    )
  
  #统计不同类型SNP的数量，并建立不同类型SNP的列表，存入对应SNPs
  final_result %>%
    count(match_type) %>%
    tidyr::pivot_wider(names_from = match_type,values_from = n,values_fill = 0) %>%
    { df <- .; miss <- setdiff(c("exact","no_match","range_50kb"), names(df)); if(length(miss)>0) df[miss] <- 0; df }%>%
    mutate(tissue = stringr::str_match(.x,"PigGTEx_v0\\.significant_lncQTL/([^/]+)\\.cis_qtl")[,2]) %>%
    mutate(
      exact_snp = list(final_result %>% filter(match_type == "exact") %>% pull(snp)),
      range50kb_snp = list(final_result %>% filter(match_type == "range_50kb") %>% pull(snp)),
      er_snp = list(final_result %>% filter(match_type %in% c("exact","range_50kb")) %>% pull(snp))) %>%
    select(tissue, exact, no_match, range_50kb, exact_snp, range50kb_snp, er_snp) %>%
    return()
}

statistic <- map_dfr(fl, st, .progress = T)
saveRDS(statistic, paste0("runs/significant_lncQTL_d", arg, ".rds"))
