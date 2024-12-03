
sort_pop_by_rst <- function(data_pop){
  # 按 rst 排
  data_pop_sorted <- data_pop[order(data_pop$rst, decreasing=T),]
  return(data_pop_sorted)
}

sort_pop_by_gene_sum <- function (data_pop){
  # 按阳性等位基因总数
  data_pop_sorted <- data_pop[order(data_pop$gene_sum, decreasing=T),]
  return(data_pop_sorted)
}

sort_pop_by_none_zero <- function (data_pop){
  # 按阳性表型数量
  data_pop_sorted <- data_pop[order(data_pop$none_zero, decreasing=T),]
  return(data_pop_sorted)
}

sort_pop_by_inbreed <- function (femaleList, maleID, relMatrix){
  # 按照近交系数
  female_list_df <- data.frame(f_id=unlist(femaleList))
  female_list_df$rel <- NA
  n <- nrow(female_list_df)
  for (f in 1:n){
    female_list_df$rel[f] <- get_relation(maleID, female_list_df$f_id[f], relMatrix)
  }
  female_list_df <- female_list_df[order(female_list_df$rel, decreasing=F), ]$f_id
  return(female_list_df)
}

