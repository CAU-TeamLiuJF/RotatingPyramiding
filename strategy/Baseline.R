
edit_male_per_family <- 5

for (k in 6:(maxGeneration - 1)){
  chosen_males <- list()
  ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)
  print(paste0('第', k, '代交配...'))
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

