setwd("/public/home/liujf/workspace/shikp/TempWork/20250612_piggtex")
library(tidyverse)
library(rrBLUP)
library(lme4)
library(lme4qtl)
library(cli)
library(sommer)

# receive parameter
arg = as.numeric(commandArgs(t = T)[1])
print(arg)
load(paste0("data/data5",".Rdata"))


# 去掉基因型文件geno列名的最后两个字符
colnames(geno) <- str_sub(colnames(geno), 1, -3)
colnames(geno)[1] <- "id"
geno <- geno %>% column_to_rownames("id")
geno[is.na(geno)] = matrix(
  colMeans(geno, na.rm = T),
  nr = nrow(geno),
  nc = ncol(geno),
  byrow = T
)[is.na(geno)]
geno <- geno %>% rownames_to_column("id")



## functions
# crva_lme4qtl <- function(
#   pheno,
#   pcol = 2,#默认值是第二列，后面ccff函数会覆盖
#   pid = 1,
#   G,
#   G2 = NULL,
#   fold = 5,
#   rep = 2,
#   pb = FALSE
# ) {
#   if (is.null(pheno)) {
#     return(NULL)
#   }
#   #G=G_sig;G2=G_res
#   # 筛选样本并准备数据
# 
#   ids <- intersect(colnames(G), pull(pheno, pid))
# 
#   phenos <- pheno %>%
#     dplyr::select(pid = pid, any_of(pcol)) %>%
#     filter(!is.na(!!sym(colnames(.)[2])), pid %in% ids) %>%
#     mutate(gid = pid, gid = as.character(gid), gid2 = gid) %>%
#     as.data.frame()
# 
#   G <- G[phenos$gid, phenos$gid]
#   if (!is.null(G2)) {
#     G2 <- G2[phenos$gid, phenos$gid]
#   }
# 
#   results <- data.frame()
# 
#   if (pb) {
#     cli::cli_progress_bar("iterator", total = fold * rep)
#   }
#   for (i in seq_len(rep)) {
#     print(paste("Starting rep", i))
#     for (j in seq_len(fold)) {
#       print(paste("  Fold", j, "of", fold, "- rep", i))
# 	#i=1;j=1
#       set.seed(i)
#       phenos$partition <- sample(
#         seq_len(fold),
#         size = nrow(phenos),
#         replace = TRUE,
#         prob = c(rep((1 / fold), times = fold))
#       )
#       phenos$yNA <- phenos[, pcol]
#       nas <- phenos$partition == j
#       phenos$yNA[nas] <- NA
# 
#       goo <- G[!nas, !nas]
#       gno <- G[nas, !nas]
#       too <- if(!is.null(G2)) G2[!nas, !nas] else NULL
#       tno <- if(!is.null(G2)) G2[nas, !nas] else NULL
# 
#       # 拟合模型
#       if(!is.null(G2)) {
# 		res_gt <- relmatLmer(
# 			yNA ~ (1 | gid) + (1 | gid2),
# 			data = phenos[!nas, ],
# 			relmat = list(gid = G, gid2 = G2)
# 		)
# 	  } else {
#         res_gt <- relmatLmer(
#           yNA ~ (1 | gid),
#           data = phenos[!nas, ],
#           relmat = list(gid = G)
#         )
#       }
# 
#       u1 <- as.matrix(t(res_gt@optinfo$relmat$relfac$gid) %*% as.matrix(ranef(res_gt)$gid))[rownames(goo), 1]
#       if(!is.null(G2)) {
# 		u2 <- as.matrix(t(res_gt@optinfo$relmat$relfac$gid2) %*% as.matrix(ranef(res_gt)$gid2))[rownames(too), 1]
# 	  } else {
#         u2 <- 0
#       }
#       # 预测
#       if(!is.null(G2)) {
# 		pred_lme4qtl_gt <- gno %*% MASS::ginv(goo) %*% u1 + tno %*% MASS::ginv(too) %*% u2
# 	  } else {
#         pred_lme4qtl_gt <- gno %*% MASS::ginv(goo) %*% u1
#       }
#       # 保存结果
#       results <- rbind(
#         results,
#         data.frame(
#           rep = i,
#           fold = j,
#           cor = cor(
#             phenos[nas, pcol],
#             pred_lme4qtl_gt,
#             use = 'pairwise.complete.obs'
#           ),
#           bias = lm(phenos[nas, pcol] ~ pred_lme4qtl_gt)$coefficients[2]
#         )
#       )
#       
#       # 打印每折完成情况
#       print(paste("    Completed fold", j, "of", fold, "- rep", i))
#       
#       if (pb) cli::cli_progress_update()
#     }
#   }
# 
#   if (pb) {
#     cli::cli_process_done()
#   }
#   return(results)
# }
#system.time(a <- crva(pheno, G = G_sig, G2 = G_res))

##——crva_sommer
crva_sommer <- function(
    pheno,
    pcol = 2,#默认值是第二列，后面ccff函数会覆盖
    pid = 1,
    G,
    G2 = NULL,
    fold = 5,
    rep = 2,
    pb = FALSE
) {
  if (is.null(pheno)) {
    return(NULL)
  }
  #G=G_sig;G2=G_res
  # 筛选样本并准备数据
  
  ids <- intersect(colnames(G), pull(pheno, pid))
  
  phenos <- pheno %>%
    dplyr::select(pid = pid, any_of(pcol)) %>%
    filter(!is.na(!!sym(colnames(.)[2])), pid %in% ids) %>%
    mutate(gid = pid, gid = as.character(gid), gid2 = gid) %>%
    as.data.frame()
  
  G <- G[phenos$gid, phenos$gid]
  if (!is.null(G2)) {
    G2 <- G2[phenos$gid, phenos$gid]
  }
  
  results <- data.frame()
  
  if (pb) {
    cli::cli_progress_bar("iterator", total = fold * rep)
  }
  for (i in seq_len(rep)) {
    print(paste("Starting rep", i))
    for (j in seq_len(fold)) {
      print(paste("  Fold", j, "of", fold, "- rep", i))
      #i=1;j=1
      set.seed(i)
      phenos$partition <- sample(
        seq_len(fold),
        size = nrow(phenos),
        replace = TRUE,
        prob = c(rep((1 / fold), times = fold))
      )
      phenos$yNA <- phenos[, pcol]
      nas <- phenos$partition == j
      phenos$yNA[nas] <- NA
      
      
      # 拟合模型
      df_test  <- phenos[nas, ]   # 测试集
      
      if (!is.null(G2)) {
        mod <- mmer(
          yNA ~ 1,
          random = ~ vs(gid, Gu = G) + vs(gid2, Gu = G2),
          data = phenos
        )
      } else {
        mod <- mmer(
          yNA ~ 1,
          random = ~ vs(gid, Gu = G),
          data = phenos
        )
      }
      u1 <- mod$U$`u:gid`$yNA[as.character(phenos$gid)]
      if (!is.null(G2)) {
        u2 <- mod$U$`u:gid2`$yNA[as.character(phenos$gid)]
        pred <- u1[nas] + u2[nas]
      } else {
        pred <- u1[nas]
      }
      pred <- as.numeric(pred)
      
      
      # 保存结果
      results <- rbind(
        results,
        data.frame(
          rep = i,
          fold = j,
          cor = cor(df_test[, pcol], pred, use = "pairwise.complete.obs"),
          bias = lm(df_test[, pcol] ~ pred)$coefficients[2]
        )
      )
      
      # 打印每折完成情况
      print(paste("    Completed fold", j, "of", fold, "- rep", i))
      
      if (pb) cli::cli_progress_update()
    }
  }
  
  if (pb) {
    cli::cli_process_done()
  }
  return(results)
}




ccff <- function(.x) {
  test <- statistic %>%
    filter(tissue == result$tissue[.x]) %>%
    pull(er_snp) %>%
    unlist()
  all_snps <- colnames(geno)[-1]  # 去掉 id
  if (length(test) == 0) {
    # erSNP为空，只用剩余 SNP
    geno_res <- geno %>%
      dplyr::select(any_of(c("id", all_snps))) %>%
      column_to_rownames("id")
    
    G_sig <- NULL
    G_res <- rrBLUP::A.mat(geno_res)
    diag(G_res) <- diag(G_res) + 0.001
  } else {
    # erSNP有值
    geno_sig <- geno %>%
      dplyr::select(any_of(c("id", test))) %>%
      column_to_rownames("id")
    
    geno_res <- geno %>%
      dplyr::select(c("id", any_of(all_snps[!all_snps %in% test]))) %>%
      column_to_rownames("id")
    
    G_sig <- rrBLUP::A.mat(geno_sig)
    G_res <- rrBLUP::A.mat(geno_res)
    diag(G_sig) <- diag(G_sig) + 0.001
    diag(G_res) <- diag(G_res) + 0.001
  }

  return(crva_sommer(
    pheno = pheno,
    pcol = result$trait[.x],
    G = G_res,
    G2 = G_sig,
    fold = 5,
    rep = 2,
    pb = T
  ))
}


statistic <- readRDS(paste0("runs/significant_lncQTL_d5", ".rds"))      

traits <- c("AGE", "BF", "TNB")    #data1的性状


result <- expand_grid(
  trait = traits,
  tissue = pull(statistic, tissue)
)


#tmp <- ccff(arg)

tmp <- walk(seq_len(nrow(result)), ~{
  if (.x == arg) {
  trait_i <- result$trait[.x]
  tissue_i <- result$tissue[.x]
  filepth2 <- paste0("runs/d5", "/", trait_i, "_", tissue_i, "_lncQTL.csv");
	write_csv(ccff(.x), filepth2);
	cli_alert_success(paste0("file saved to : ", filepth2))
  }
})

