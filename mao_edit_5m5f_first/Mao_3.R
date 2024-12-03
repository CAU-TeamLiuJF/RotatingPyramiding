
edit_male_per_family <- 5

# 第 6 代给每个家系编辑全部个体
k <- 6
ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)
for (f in 1:familyNum){
  chosen_family <- all_familys[f]
  chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
  sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
  sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
  chosen_male <- sample_pop_male$id_pop[1:malePerFamily]

  # 编辑公猪和母猪
  edit_ind <- sapply(chosen_pop_ids, function(x) as.numeric(x) - ind_offset)
  data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_ind,
                                   chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

  sample_pop_female <- sample_pop[sample_pop[2] == 'F', ]
  chosen_females <- sample_pop_female$id_pop[1:edit_male_per_family]
  edit_female_ind <- sapply(chosen_females, function(x) as.numeric(x) - ind_offset)
  data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                   chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)
}

cross_times <- ceiling(log(targetGeneNum, 2))
for (k in 6:(6 + cross_times - 2)){
  cross_interval <- (2^(k - 6))
  print(paste0('第', k, '代交配(基因整合)...'))
  all_female_list <- data_list[[k]]$data_pop[data_list[[k]]$data_pop[2] == 'F', ]
  ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)
  for (f in 1:familyNum){

    male_family <- all_familys[f]
    female_family_idx <- (f+cross_interval) %% familyNum
    if (female_family_idx == 0){
      female_family_idx <- 5
    }
    female_family <- all_familys[female_family_idx]
    male_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == male_family, 1]
    male_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% male_pop_ids, ]
    chosen_males <- male_pop[male_pop[2] == 'M', ]$id_pop[1:malePerFamily]

    female_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == female_family, 1]
    female_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% female_pop_ids, ]
    chosen_females <- female_pop[female_pop[2] == 'F', 1]

    chosen_cross_plan <- cross_plan_homo(chosen_males, chosen_females, data_list[[k]]$pedi)
    if (f == 1){
      cross_plan <- chosen_cross_plan
    } else{
      cross_plan <- rbind(cross_plan, chosen_cross_plan)
    }
  }
  if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
    stop('cross plan 数量不正确!!!')
  }
  data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)
}

for (k in (6 + cross_times - 1):15){
  ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)
  chosen_males <- list()
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

