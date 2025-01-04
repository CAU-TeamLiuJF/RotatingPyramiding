# 关闭 warning
options(warn=-1)
# 加载各种包
library(optparse)
library(AlphaSimR)
library(nadiv)
# 导入基础公共参数
source('common.R')
# 设置各种路径
work_dir <- getwd()
lib_dir <- file.path(work_dir, 'lib')
params_dir <- file.path(work_dir, 'params')
strategy_dir <- file.path(work_dir, 'strategy')

# 导入各种函数
source(file.path(lib_dir, 'sort.R'))
source(file.path(lib_dir, 'functions.R'))
source(file.path(lib_dir, 'cross.R'))

# 获取命令行参数
option_parser <- OptionParser()
# 使用的参数卡文件
option_parser <- add_option(option_parser,
                            c("-p", "--params"),
                            dest="params_file",
                            type="character",
                            help="params file",
                            default='example.txt')
# 重复次数
option_parser <- add_option(option_parser,
                            c("-r", "--repeat"),
                            dest="repeat_num",
                            type="integer",
                            help="repeat num",
                            default=20)

option_parser <- add_option(option_parser,
                            c("-s", "--name"),
                            dest="strategy_name",
                            type="character",
                            help="strategy dir",
                            default='baseline')

# 输出目录
option_parser <- add_option(option_parser,
                            c("-o", "--output"),
                            dest="output_name",
                            type="character",
                            help="output dir",
                            default='default_output')
# 解析命令行参数
parsed_args <- parse_args(option_parser)
params_file <- parsed_args$params_file
params_file_path <- file.path(params_dir, params_file)
repeat_num <- parsed_args$repeat_num
strategy_name <- parsed_args$strategy_name
output_name <- parsed_args$output_name
# 读取参数文件
print(params_file_path)
data <- read.table(params_file_path, header=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8")
# 将第一列作为变量名，第二列作为变量值
variable_names <- data[, 1]
variable_values <- data[, 2]
# 创建变量
for (i in seq_along(variable_names)) {
  assign(variable_names[i], variable_values[i])
}

# 一些文件的保存路径
tmp_dir <- file.path(work_dir, 'tmp')
if (!dir.exists(tmp_dir)){
  dir.create(tmp_dir)
}
tmp_dir <- file.path(tmp_dir, output_name)
if (!dir.exists(tmp_dir)){
  dir.create(tmp_dir)
}
# 一些文件的保存路径
output_dir <- file.path(work_dir, 'output')
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}
output_dir <- file.path(output_dir, output_name)
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}

# 初始化群体
init_data_file <- file.path(tmp_dir, 'init_data.RData')
if (file.exists(init_data_file)){
  print('已有初始群体, 直接读取...')
  load(init_data_file)
} else{
  source('init_pop.R')
}

for (repeat_idx in 1:repeat_num){

  tmp_result_file <- paste0('result_', strategy_name, '_repeat_', repeat_idx, '.RData')
  tmp_result_file <- file.path(tmp_dir, tmp_result_file)
  print(paste0('目标基因数 ', targetGeneNum, ' 情景 ', strategy_name, ' 第 ', repeat_idx, ' 次重复...'))
  # 断点续跑
  print(tmp_result_file)
  if (file.exists(tmp_result_file)){
    print("file exists, skip")
    next
  }

  # 前五代
  data_5_file <- file.path(tmp_dir, 'data_5.RData')
  print(data_5_file)
  if (file.exists(data_5_file)){
    print('已有前 5 代群体, 直接读取...')
    load(data_5_file)
  } else{
    source('prepare_pop.R')
  }

  # 从第 6 代开始各种情况
  scenario_dir <- file.path(strategy_dir, paste0(strategy_name, '.R'))
  print(scenario_dir)
  source(scenario_dir)

  save(data_list, file=tmp_result_file)
}



