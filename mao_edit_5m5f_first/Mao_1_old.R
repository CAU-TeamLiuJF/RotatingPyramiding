
edit_male_per_family <- 5

# 6-10 每代选择一个家系的 edit_male_per_family 头公猪导入一个新的目标基因，分五代导入五个基因到五个家系
# 11-15 正常按同质选配
for (k in 6:15){

  chosen_males <- list()
  ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)
  print(paste0('第', k, '代交配...'))

  for (f in 1:familyNum){
    # 6-9 代选择第 k-5 个家系的全部个体导入目标基因 第 6 代额外还要再编辑一个家系
    if ((k == 6 & (f == 1)) | (k == 6 & f == 2) | (k <= 9 & k > 6 & (f == (k - 4)))){
      chosen_family <- all_familys[f]
      chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
      sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
      sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
      chosen_male <- sample_pop_male$id_pop[1:malePerFamily]
      chosen_males <- c(chosen_males, chosen_male)

      # 编辑公猪和母猪各 edit_male_per_family 头
      edit_male_ind <- sapply(chosen_male[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
      data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_male_ind,
                                       chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

      sample_pop_female <- sample_pop[sample_pop[2] == 'F', ]
      chosen_females <- sample_pop_female$id_pop
      edit_female_ind <- sapply(chosen_male, function(x) as.numeric(x) - ind_offset)
      data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                       chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

    } else {
      chosen_family <- all_familys[f]
      chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
      sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
      sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
      chosen_male <- sample_pop_male$id_pop[1:malePerFamily]
      chosen_males <- c(chosen_males, chosen_male)
    }
    sample_pop_female <- sample_pop[sample_pop[2] == 'F', ]
    chosen_females <- sample_pop_female$id_pop
    chosen_cross_plan <- cross_plan_homo(chosen_male, chosen_females, data_list[[k]]$pedi)
    if (f == 1){
      cross_plan <- chosen_cross_plan
    } else{
      cross_plan <- rbind(cross_plan, chosen_cross_plan)
    }
  }

  if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
    print(nrow(cross_plan))
    stop('cross plan 数量不正确!!!')
  }
  data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)
}
