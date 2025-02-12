options(warn=-1)
library(optparse)
library(AlphaSimR)
library(nadiv)
# 导入基础公共参数
source('common.R')

work_dir <- getwd()

# 获取命令行参数
option_parser <- OptionParser()
# 情景编号
option_parser <- add_option(option_parser,
                            c("-s", "--strategy"),
                            dest="strategy_name",
                            type="character",
                            help="strategy name",
                            default="baseline")

# 重复次数
option_parser <- add_option(option_parser,
                            c("-r", "--repeat"),
                            dest="repeat_num",
                            type="integer",
                            help="repeat num",
                            default=20)

# 输出目录
option_parser <- add_option(option_parser,
                            c("-o", "--output"),
                            dest="output_name",
                            type="character",
                            help="output dir",
                            default='default_output')
# 解析命令行参数
parsed_args <- parse_args(option_parser)
strategy_name <- parsed_args$strategy_name
repeat_num <- parsed_args$repeat_num
output_name <- parsed_args$output_name

output_dir <- file.path(work_dir, 'output', output_name)
tmp_dir <- file.path(work_dir, 'tmp', output_name)
init_data_file <- file.path(tmp_dir, 'init_data.RData')
load(file=init_data_file)

geno_dir <- file.path(output_dir, 'geno')
if (!dir.exists(geno_dir)){
  dir.create(geno_dir)
}

print(target_gene_pos)

# 获得全基因组信息
geno_map <- data.frame(getGenMap(founderPop, sex='A'))
geno_map_file <- paste0('geno_map_', strategy_name, '.txt')
write.table(geno_map, file.path(geno_dir, geno_map_file), quote=F, row.names=F, col.names=T)
# 获得 snp map
snp_map <- data.frame(getSnpMap(snpChip = 1, sex='A', simParam = SP))
snp_map_file <- paste0('snp_map_', strategy_name, '.txt')
write.table(snp_map, file.path(geno_dir, snp_map_file), quote=F, row.names=F, col.names=T)
# 获得 qtl map
snp_map <- data.frame(getQtlMap(trait = 1, sex='A', simParam = SP))
snp_map_file <- paste0('qtl_map_', strategy_name, '.txt')
write.table(snp_map, file.path(geno_dir, snp_map_file), quote=F, row.names=F, col.names=T)


for (repeat_id in 1:repeat_num){
  tmp_result_file <- paste0('result_', strategy_name, '_repeat_', repeat_id, '.RData')
  tmp_result_file <- file.path(tmp_dir, tmp_result_file)
  print(tmp_result_file)
  print(paste0('reading result of strategy [', strategy_name, '] repeat ', repeat_id, '...'))
  load(tmp_result_file)
  result <- data.frame(matrix(ncol = 17, nrow = 0))
  for (k in 1:maxGeneration){
    print(paste0('generation ', k))
    if (repeat_id > 1 && k < 6){
      next
    }

    # inb
    avg_inb <- 0
    if (k >= 6){
      print(paste0(format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), ' Calculating inb...'))
      pedi_total <- data_list[[k-2]]$pedi
      for (l in (k-2):k){
        pedi_total <- rbind(pedi_total, data_list[[l]]$pedi)
      }
      pedi_total <- pedi_total[, c(1, 2, 3)]
      pped <- prepPed(pedi_total)
      A <- as.matrix(makeA(pped))
      A_row <- nrow(A)
      inb_pop <- diag(A) - 1
      # 只取当前代的部分
      avg_inb <- mean(inb_pop[(A_row - n_rows):A_row])
    }

    # 目标基因含量
    data <- data_list[[k]]
    data_pop <- data$data_pop
    n_rows <- nrow(data_pop)
    avg_target_percent <- (mean(data_pop$gene_sum) / (2 * targetGeneNum))

    # 计算每个性状 gv 的均值和方差
    mean_gv_CZS <- mean(data$gv$CZS)
    var_gv_CZS <- var(data$gv$CZS)

    mean_gv_JZRL <- mean(data$gv$JZRL)
    var_gv_JZRL <- var(data$gv$JZRL)

    mean_gv_JZBBH <- mean(data$gv$JZBBH)
    var_gv_JZBBH <- var(data$gv$JZBBH)

    # 计算 rst
    for (i in 1:n_rows) {
      data$data_pop$rst[i] <- w_tnb*(data$gv$CZS[i]/var_gv_CZS) - w_age*(data$gv$JZRL[i]/var_gv_JZRL) - w_bf*(data$gv$JZBBH[i]/var_gv_JZBBH)
    }

    mean_rst <- mean(data$data_pop$rst)
    var_rst <- var(data$data_pop$rst)

    # 计算每个性状 pheno 的均值和方差
    mean_pheno_CZS <- mean(data$pheno$CZS)
    var_pheno_CZS <- var(data$pheno$CZS)

    mean_pheno_JZRL <- mean(data$pheno$JZRL)
    var_pheno_JZRL <- var(data$pheno$JZRL)

    mean_pheno_JZBBH <- mean(data$pheno$JZBBH)
    var_pheno_JZBBH <- var(data$pheno$JZBBH)

    # 不同阳性基因数的个体分布
    target_gene_str <- ''
    for (t in 1:targetGeneNum){
      target_gene_pop <- data_pop[data_pop$none_zero == t,]
      n_pop <- nrow(target_gene_pop)
      target_gene_str <- paste0(target_gene_str, n_pop, ',')
    }

    # 获取 snp 信息计算 MAF

    print(paste0(format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), ' Calculating maf...'))
    snp_geno <- data.frame(pullSnpGeno(data$pop, simParam = SP))

    col_sums <- lapply(snp_geno, sum)
    MAF <- data.frame(
      POS = names(col_sums),
      MAF = unlist(col_sums)
    )
    MAF$MAF <- (2 * n_rows - MAF$MAF) / (2 * n_rows)
    MAF$MAF <- ifelse(MAF$MAF > 0.5, 1 - MAF$MAF, MAF$MAF)
    maf_file <- paste0('maf_strategy_', strategy_name, '_repeat_', repeat_id, '_generation_', k, '.maf')
    write.table(MAF, file.path(geno_dir, maf_file), quote=F, row.names=F, col.names=T)
    rm(snp_geno)

    avg_MAF <- mean(MAF$MAF)
    result <- rbind(result, c(avg_target_percent, mean_gv_CZS, var_gv_CZS, mean_gv_JZRL, var_gv_JZRL, mean_gv_JZBBH,
                              var_gv_JZBBH, mean_rst,
                              mean_pheno_CZS, var_pheno_CZS, mean_pheno_JZRL, var_pheno_JZRL,
                              mean_pheno_JZBBH, var_pheno_JZBBH, avg_MAF, avg_inb, target_gene_str))
  }
  colnames(result) <- c("avg_target_percent", "mean_gv_CZS", "var_gv_CZS", "mean_gv_JZRL", "var_gv_JZRL", "mean_gv_JZBBH",
                        "var_gv_JZBBH", "mean_rst",
                        "mean_pheno_CZS", "var_pheno_CZS", "mean_pheno_JZRL", "var_pheno_JZRL",
                        "mean_pheno_JZBBH", "var_pheno_JZBBH", "avg_MAF", "avg_inb", "target_gene_str")
  result_file <- paste0('result_strategy_', strategy_name, '_repeat_', repeat_id, '.csv')
  write.table(result, file=file.path(output_dir, result_file), sep=',', quote=F, row.names=F, col.names=T)
}
