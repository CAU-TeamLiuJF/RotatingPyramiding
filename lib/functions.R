count_nonzero <- function(row) {
  sum(row != 0)
}

get_pop_data <- function (pop){
  pedi_pop <- getPed(pop)
  pheno_pop <- pheno(pop)
  pheno_pop <- as.data.frame(pheno_pop)
  id_pop <- pop@id
  sex_pop <- pop@sex
  id_pop <- as.data.frame(id_pop)
  sex_pop <- as.data.frame(sex_pop)
  data_pop <- data.frame(id_pop, sex_pop)
  data <- NA
  data$pop <- pop
  data$data_pop <- data_pop
  data$pedi <- pedi_pop
  data$pheno <- pheno_pop
  gv <- gv(data$pop)
  gv <- data.frame(gv)
  data$gv <- gv
  names(data$data_pop) <- c('id_pop', 'sex_pop')

  # 获取每个目标基因的基因型
  for (g in 1:targetGeneNum){
    gene_check <- pullMarkerGeno(data$pop, markers=target_gene_pos[[g]])
    new_column_name <- paste0('gene_', g)
    data$data_pop[[new_column_name]] <- gene_check
  }
  # 计算目标基因的阳性等位基因和
  data$data_pop$gene_sum <- rowSums(data$data_pop[, 3:(targetGeneNum+2)])
  data$data_pop$none_zero <- apply(data$data_pop[, 3:(targetGeneNum+2)], 1, count_nonzero)
  data <- get_rst(data)
  # 默认按照目标基因数量排序后按rst排序
  data$data_pop <- sort_pop_by_none_zero(sort_pop_by_gene_sum(sort_pop_by_rst(data$data_pop)))
  # print(head(data$data_pop))
  return(data)
}

get_rst <- function(data){
  n <- nrow(data$data_pop)
  data$data_pop$rst <- NA
  data$data_pop$rst_ebv <- NA
  data$data_pop$ebv_CZS <- NA
  data$data_pop$ebv_JZRL <- NA
  data$data_pop$ebv_JZBBH <- NA

  #计算每个性状育种值的均值和标准差
  mean_CZS <- mean(data$gv$CZS)
  sd_CZS <- sd(data$gv$CZS)

  mean_JZRL <- mean(data$gv$JZRL)
  sd_JZRL <- sd(data$gv$JZRL)

  mean_JZBBH <- mean(data$gv$JZBBH)
  sd_JZBBH <- sd(data$gv$JZBBH)
  e_norm <- rnorm(n)

  #根据评估准确性为每个性状模拟出估计育种值ebv
  for (i in 1:n) {
    data$data_pop$ebv_CZS[i] <- mean_CZS + r_CZS * (data$gv$CZS[i] - mean_CZS) + sqrt(1 - r_CZS * r_CZS) * e_norm[i] * sd_CZS
    data$data_pop$ebv_JZRL[i] <- mean_JZRL + r_JZRL * (data$gv$JZRL[i] - mean_JZRL) + sqrt(1 - r_JZRL * r_JZRL) * e_norm[i] * sd_JZRL
    data$data_pop$ebv_JZBBH[i] <- mean_JZBBH + r_JZBBH * (data$gv$JZBBH[i] - mean_JZBBH) + sqrt(1 - r_JZBBH * r_JZBBH) * e_norm[i] * sd_JZBBH
  }

  #计算模拟出的ebv的均值方差
  #计算模拟出的ebv的均值方差
  mean_CZS_ebv <- mean(data$data_pop$ebv_CZS)
  sd_CZS_ebv <- sd(data$data_pop$ebv_CZS)

  mean_JZRL_ebv <- mean(data$data_pop$ebv_JZRL)
  sd_JZRL_ebv <- sd(data$data_pop$ebv_JZRL)

  mean_JZBBH_ebv <- mean(data$data_pop$ebv_JZBBH)
  sd_JZBBH_ebv <- sd(data$data_pop$ebv_JZBBH)

  # print(cor(data$data_pop$ebv_CZS, data$gv$CZS))
  # print(cor(data$data_pop$ebv_JZRL, data$gv$JZRL))
  # print(cor(data$data_pop$ebv_JZBBH, data$gv$JZBBH))

  #根据ebv计算rst用于选种排序
  for (i in 1:n) {
      data$data_pop$rst[i] <- w_tnb*((data$data_pop$ebv_CZS[i]-mean_CZS_ebv)/sd_CZS_ebv)-w_age*((data$data_pop$ebv_JZRL[i]-mean_JZRL_ebv)/sd_JZRL_ebv)-w_bf*((data$data_pop$ebv_JZBBH[i]-mean_JZBBH_ebv)/sd_JZBBH_ebv);
  }
  return(data)
}

get_rel_matrix <- function(iteration, back=2){
  # 获取亲缘关系矩阵
  pedi_total <- data_list[[iteration-back]]$pedi
  for (l in (iteration-(back-1)):iteration){
    pedi_total <- rbind(pedi_total, data_list[[l]]$pedi)
  }
  pedi_total <- pedi_total[, c(1, 2, 3)]
  pped <- prepPed(pedi_total)
  A <- as.matrix(makeA(pped))
  A <- as.data.frame(A)
  return(A)
}

get_relation <- function(male, female, relMatrix){
  return(relMatrix[male, female])
}
