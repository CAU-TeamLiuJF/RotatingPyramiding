options(warn=-1)
library(optparse)
library(AlphaSimR)
# 导入基础公共参数
source('common.R')
# 表型权重
w_age <- 0.4
w_bf <- 0.3
w_tnb <- 0.3

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
first5_data_file <- file.path(tmp_dir, 'data_5.RData')
load(file=first5_data_file)

dis_dir <- file.path(output_dir, 'gene_dis')
if (!dir.exists(dis_dir)){
  dir.create(dis_dir)
}

print(target_gene_pos)
print(all_familys)
familyNum <- 5
for (repeat_id in 1:repeat_num){
  tmp_result_file <- paste0('result_', strategy_name, '_repeat_', repeat_id, '.RData')
  tmp_result_file <- file.path(tmp_dir, tmp_result_file)
  print(tmp_result_file)
  print(paste0('reading result of scenario [', strategy_name, '] repeat ', repeat_id, '...'))
  load(tmp_result_file)
  result <- data.frame(matrix(ncol = 2, nrow = 0))


  for (k in 6:maxGeneration){

    # 不同阳性基因数的个体分布
    target_gene_dis <- ''

    for (f in 1:familyNum){
      chosen_family <- all_familys[f]
      chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
      chosen_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
      for (t in 1:targetGeneNum){
        target_gene_pop <- chosen_pop[chosen_pop$none_zero == t,]
        n_pop <- nrow(target_gene_pop)
        target_gene_dis <- paste0(target_gene_dis, n_pop, ':')
      }
      target_gene_dis <- paste0(target_gene_dis, ',')
    }

    result <- rbind(result, c(paste0('G', k - 6), target_gene_dis))
  }
  col_names <- c('G')
  for (f in 1:familyNum){
    col_names <- c(col_names, paste0('f', f))
  }
  colnames(result) <- c('G', "target_gene_dis")
  result_file <- paste0('target_gene_dis_strategy_', strategy_name, '_repeat_', repeat_id, '.csv')
  write.table(result, file=file.path(dis_dir, result_file), sep=',', quote=F, row.names=F, col.names=T)
}
