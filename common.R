# 关闭 warning
options(warn=-1)

# 加载各种包
library(optparse)
library(AlphaSimR)
library(nadiv)

# 共用的各种参数
# 表型参数矩阵输入
corA <- matrix(c(1,0,0,0,1,-0.3,0,-0.3,1),nrow=3, ncol=3)
mean <- matrix(c(10, 180, 20), nrow=3, ncol=1)
var <- matrix(c(3, 80, 2), nrow=3, ncol=1)
h2 <- matrix(c(0.1, 0.3, 0.4), nrow=3, ncol=1)
names <- matrix(c("CZS","JZRL", "JZBBH"), nrow=3, ncol=1)
# 表型权重
w_age <- 0.4
w_bf <- 0.3
w_tnb <- 0.3
r_CZS <- 0.2
r_JZRL <- 0.4
r_JZBBH <- 0.5
