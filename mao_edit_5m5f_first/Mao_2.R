
edit_male_per_family <- 5

# 第 6 代给 4 个家系编辑 edit_male_per_family头公猪，然后与家系内母猪杂交
k <- 6
print(paste0('第', k, '代交配...'))
ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)
for (f in 1:5){
  chosen_family <- all_familys[f]
  chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
  sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
  sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
  chosen_male <- sample_pop_male$id_pop[1:malePerFamily]

  # 给 4 个家系编辑
  if (f <= 4){

    # 编辑公猪和母猪各 edit_male_per_family 头
    edit_male_ind <- sapply(chosen_male[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
    data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_male_ind,
                                     chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

    sample_pop_female <- sample_pop[sample_pop[2] == 'F', ]
    chosen_females <- sample_pop_female$id_pop[1:edit_male_per_family]
    edit_female_ind <- sapply(chosen_females, function(x) as.numeric(x) - ind_offset)
    data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                     chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)
  }

}

# 家系 1 与 家系 3 杂交, 家系 2 与 家系 4 杂交
for (f in 1:2){
  family_a <- all_familys[f]
  family_b <- all_familys[f+2]
  family_a_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == family_a, 1]
  family_b_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == family_b, 1]
  family_a_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% family_a_ids, ]
  family_b_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% family_b_ids, ]
  family_a_male <- family_a_pop[family_a_pop[2] == 'M', ]$id_pop[1:malePerFamily]
  family_b_male <- family_b_pop[family_b_pop[2] == 'M', ]$id_pop[1:malePerFamily]
  family_a_female <- family_a_pop[family_a_pop[2] == 'F', 1]
  family_b_female <- family_b_pop[family_b_pop[2] == 'F', 1]

  cross_plan_a <- cross_plan_homo(family_a_male, family_b_female, data_list[[k]]$pedi)
  cross_plan_b <- cross_plan_homo(family_b_male, family_a_female, data_list[[k]]$pedi)
  chosen_cross_plan <- rbind(cross_plan_a, cross_plan_b)
  if (f == 1){
    cross_plan <- chosen_cross_plan
  } else{
    cross_plan <- rbind(cross_plan, chosen_cross_plan)
  }
}
# 家系 5 不作任何操作
family_5 <- all_familys[5]
chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == family_5, 1]
chosen_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
chosen_males <- chosen_pop[chosen_pop[2] == 'M', ]$id_pop[1:malePerFamily]
chosen_females <- chosen_pop[chosen_pop[2] == 'F', 1]
chosen_cross_plan <- cross_plan_homo(chosen_males, chosen_females, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan, chosen_cross_plan)
if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
  print(nrow(cross_plan))
  stop('cross plan 数量不正确!!!')
}
data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)

# 家系 1 和 3 为 AaBb, 家系 2 和 4 为 CcDd, 家系 5 不做操作
# 家系 1,3 与 家系 2,4 交叉选配
k <- 7
ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)
print(paste0('第', k, '代交配...'))
family_1_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[1], 1]
family_2_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[2], 1]
family_3_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[3], 1]
family_4_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[4], 1]
family_1_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% family_1_ids, ]
family_2_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% family_2_ids, ]
family_3_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% family_3_ids, ]
family_4_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% family_4_ids, ]
male_1 <- family_1_pop[family_1_pop[2] == 'M', ]$id_pop[1:malePerFamily]
male_2 <- family_2_pop[family_2_pop[2] == 'M', ]$id_pop[1:malePerFamily]
male_3 <- family_3_pop[family_3_pop[2] == 'M', ]$id_pop[1:malePerFamily]
male_4 <- family_4_pop[family_4_pop[2] == 'M', ]$id_pop[1:malePerFamily]

male_1_3 <- list()
male_2_4 <- list()
for (i in 1:malePerFamily){
  male_1_3 <- c(male_1_3, male_1[i:i], male_3[i:i])
  male_2_4 <- c(male_2_4, male_2[i:i], male_4[i:i])
}
family_1_3_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[1] | data_list[[k]]$pedi$family == all_familys[3], 1]
family_2_4_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[2] | data_list[[k]]$pedi$family == all_familys[4], 1]
family_1_3_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% family_1_3_ids, ]
family_2_4_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% family_2_4_ids, ]
female_1_3 <- family_1_3_pop[family_1_3_pop[2] == 'F', 1]
female_2_4 <- family_2_4_pop[family_2_4_pop[2] == 'F', 1]
cross_plan_1_3 <- cross_plan_homo(male_1_3, female_2_4, data_list[[k]]$pedi)
cross_plan_2_4 <- cross_plan_homo(male_2_4, female_1_3, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan_1_3, cross_plan_2_4)

# 家系 5 不做任何操作
family_5 <- all_familys[5]
chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == family_5, 1]
chosen_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
chosen_males <- chosen_pop[chosen_pop[2] == 'M', ]$id_pop[1:malePerFamily]
chosen_females <- chosen_pop[chosen_pop[2] == 'F', 1]

chosen_cross_plan <- cross_plan_homo(chosen_males, chosen_females, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan, chosen_cross_plan)
if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
  print(nrow(cross_plan))
  stop('cross plan 数量不正确!!!')
}
data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)

# 家系 1,2,3,4 为 A_B_C_D_, 家系 5 导入为 EE 交叉选配
k <- 8
ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)
print(paste0('第', k, '代交配...'))
chosen_males_abcd <- list()
for (f in 1:4){
  chosen_family <- all_familys[f]
  chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
  chosen_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
  chosen_males <- chosen_pop[chosen_pop[2] == 'M', ]$id_pop[1:malePerFamily]
  chosen_males_abcd <- c(chosen_males_abcd, chosen_males)
}
chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family != all_familys[5], 1]
chosen_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
chosen_females_abcd <- chosen_pop[chosen_pop[2] == 'F', 1]

chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[5], 1]
chosen_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
chosen_males_e <- chosen_pop[chosen_pop[2] == 'M', ]$id_pop[1:malePerFamily]
chosen_females_e <- chosen_pop[chosen_pop[2] == 'F', 1]

# 第五个家系导入基因
# 编辑全部
edit_ind <- sapply(chosen_pop_ids, function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_ind,
                                 chr=target_gene_chr[[5]], segSites=target_gene_site[[5]], allele=1)

edit_female_ind <- sapply(chosen_females_e[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

cross_plan_abcd <- cross_plan_homo(chosen_males_abcd, chosen_females_e, data_list[[k]]$pedi)
cross_plan_e <- cross_plan_homo(chosen_males_e, chosen_females_abcd, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan_abcd, cross_plan_e)
if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
  print(nrow(cross_plan))
  stop('cross plan 数量不正确!!!')
}
data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)

for (k in 9:15){
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

