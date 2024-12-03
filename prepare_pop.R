load(file=init_data_file)

data_list <- list()
data_list[[1]] <- data_0

# 产生前五代
# 每代随机选择 familyNum 个家系, 每个家系选择 rst 最高的 malePerFamily 头公猪
k <- 1
print(paste0('第', k, '代交配...'))
chosen_males <- data_list[[k]]$data_pop[data_list[[k]]$data_pop[2] == 'M', ]$id_pop[1:(familyNum * malePerFamily)]
female_list <- data_list[[k]]$data_pop[data_list[[k]]$data_pop[2] == 'F', 1]
cross_plan <- cross_plan_homo(chosen_males, female_list, data_list[[k]]$pedi)
data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)

for (k in 2:5){
  print(paste0('第', k, '代交配...'))
  chosen_males <- list()
  all_familys <- t(as.data.frame(table(data_list[[k]]$pedi$family))[1])
  chosen_familys <- sample(all_familys, size = familyNum)
  for (chosen_family in chosen_familys){
    sample_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
    sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% sample_ids, ]
    sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
    chosen_male <- sample_pop_male$id_pop[1:malePerFamily]
    chosen_males <- c(chosen_males, chosen_male)
  }
  if (length(chosen_males) != (familyNum * malePerFamily)){
    stop('公猪数量不正确!!!')
  }
  female_list <- data_list[[k]]$data_pop[data_list[[k]]$data_pop[2] == 'F', 1]
  cross_plan <- cross_plan_homo(chosen_males, female_list, data_list[[k]]$pedi)
  data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)
}
# 保存初始环境
print(all_familys)
print('保存前五代...')
save(data_list, SP, all_familys, file=data_5_file)
print('保存前五代成功')

