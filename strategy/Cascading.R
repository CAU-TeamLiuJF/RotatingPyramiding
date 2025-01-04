
edit_male_per_family <- 5

# 第 6 代 家系 1 和 2 导入基因
k <- 6
print(paste0('第', k, '代交配...'))
ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)

f <- 1
chosen_family_1 <- all_familys[f]
chosen_pop_ids_1 <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family_1, 1]
sample_pop_1 <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids_1, ]
sample_pop_male_1 <- sample_pop_1[sample_pop_1[2] == 'M', ]
sample_pop_female_1 <- sample_pop_1[sample_pop_1[2] == 'F', ]
chosen_male_1 <- sample_pop_male_1$id_pop[1:malePerFamily]
chosen_female_1 <- sample_pop_female_1$id_pop
# 编辑公猪和母猪各 edit_male_per_family 头
edit_male_ind <- sapply(chosen_male_1[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_male_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)
edit_female_ind <- sapply(chosen_female_1[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

f <- 2
chosen_family_2 <- all_familys[f]
chosen_pop_ids_2 <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family_2, 1]
sample_pop_2 <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids_2, ]
sample_pop_male_2 <- sample_pop_2[sample_pop_2[2] == 'M', ]
sample_pop_female_2 <- sample_pop_2[sample_pop_2[2] == 'F', ]
chosen_male_2 <- sample_pop_male_2$id_pop[1:malePerFamily]
chosen_female_2 <- sample_pop_female_2$id_pop
# 编辑公猪和母猪各 edit_male_per_family 头
edit_male_ind <- sapply(chosen_male_2[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_male_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)
edit_female_ind <- sapply(chosen_female_2[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

chosen_cross_plan_1 <- cross_plan_homo(chosen_male_1, chosen_female_2, data_list[[k]]$pedi)
chosen_cross_plan_2 <- cross_plan_homo(chosen_male_2, chosen_female_1, data_list[[k]]$pedi)
cross_plan <- rbind(chosen_cross_plan_1, chosen_cross_plan_2)

chosen_males_345 <- list()
for (f in 3:5){
  chosen_family <- all_familys[f]
  chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
  sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
  sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
  chosen_male <- sample_pop_male$id_pop[1:malePerFamily]
  chosen_males_345 <- c(chosen_males_345, chosen_male)
}
sample_345_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[3] |
                                        data_list[[k]]$pedi$family == all_familys[4] | data_list[[k]]$pedi$family == all_familys[5], 1]
sample_345_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% sample_345_ids, ]
female_345_ids <- sample_345_pop[sample_345_pop[2] == 'F', 1]
cross_plan_345 <- cross_plan_homo(chosen_males_345, female_345_ids, data_list[[k]]$pedi)

cross_plan <- rbind(cross_plan, cross_plan_345)
if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
  print(nrow(cross_plan))
  stop('cross plan 数量不正确!!!')
}
data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)

# 第 7 代 家系 3 导入基因
k <- 7
print(paste0('第', k, '代交配...'))
ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)

f <- 3
chosen_family_3 <- all_familys[f]
chosen_pop_ids_3 <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family_3, 1]
sample_pop_3 <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids_3, ]
sample_pop_male_3 <- sample_pop_3[sample_pop_3[2] == 'M', ]
sample_pop_female_3 <- sample_pop_3[sample_pop_3[2] == 'F', ]
chosen_male_3 <- sample_pop_male_3$id_pop[1:malePerFamily]
chosen_female_3 <- sample_pop_female_3$id_pop
# 编辑公猪和母猪各 edit_male_per_family 头
edit_male_ind <- sapply(chosen_male_3[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_male_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)
edit_female_ind <- sapply(chosen_female_3[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

chosen_males_12 <- list()
for (f in 1:2){
  chosen_family <- all_familys[f]
  chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
  sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
  sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
  chosen_male <- sample_pop_male$id_pop[1:malePerFamily]
  chosen_males_12 <- c(chosen_males_12, chosen_male)
}
sample_12_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[1] |
                                        data_list[[k]]$pedi$family == all_familys[2], 1]
sample_12_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% sample_12_ids, ]
female_12_ids <- sample_12_pop[sample_12_pop[2] == 'F', 1]
cross_plan_12_3 <- cross_plan_homo(chosen_males_12, chosen_female_3, data_list[[k]]$pedi)
cross_plan_3_12 <- cross_plan_homo(chosen_male_3, female_12_ids, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan_12_3, cross_plan_3_12)

chosen_males_45 <- list()
for (f in 4:5){
  chosen_family <- all_familys[f]
  chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
  sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
  sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
  chosen_male <- sample_pop_male$id_pop[1:malePerFamily]
  chosen_males_45 <- c(chosen_males_45, chosen_male)
}
sample_45_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[4] | data_list[[k]]$pedi$family == all_familys[5], 1]
sample_45_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% sample_45_ids, ]
female_45_ids <- sample_45_pop[sample_45_pop[2] == 'F', 1]
cross_plan_45 <- cross_plan_homo(chosen_males_45, female_45_ids, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan, cross_plan_45)

if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
  print(nrow(cross_plan))
  stop('cross plan 数量不正确!!!')
}
data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)

# 第 8 代 家系 4 导入基因
k <- 8
print(paste0('第', k, '代交配...'))
ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)

f <- 4
chosen_family_4 <- all_familys[f]
chosen_pop_ids_4 <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family_4, 1]
sample_pop_4 <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids_4, ]
sample_pop_male_4 <- sample_pop_4[sample_pop_4[2] == 'M', ]
sample_pop_female_4 <- sample_pop_4[sample_pop_4[2] == 'F', ]
chosen_male_4 <- sample_pop_male_4$id_pop[1:malePerFamily]
chosen_female_4 <- sample_pop_female_4$id_pop
# 编辑公猪和母猪各 edit_male_per_family 头
edit_male_ind <- sapply(chosen_male_4[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_male_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)
edit_female_ind <- sapply(chosen_female_4[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

chosen_males_123 <- list()
for (f in 1:3){
  chosen_family <- all_familys[f]
  chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
  sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
  sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
  chosen_male <- sample_pop_male$id_pop[1:malePerFamily]
  chosen_males_123 <- c(chosen_males_123, chosen_male)
}
sample_123_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == all_familys[1] |
                                        data_list[[k]]$pedi$family == all_familys[2] | data_list[[k]]$pedi$family == all_familys[3], 1]
sample_123_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% sample_123_ids, ]
female_123_ids <- sample_123_pop[sample_123_pop[2] == 'F', 1]
cross_plan_123_4 <- cross_plan_homo(chosen_males_123, chosen_female_4, data_list[[k]]$pedi)
cross_plan_4_123 <- cross_plan_homo(chosen_male_4, female_123_ids, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan_123_4, cross_plan_4_123)

f <- 5
chosen_family_5 <- all_familys[f]
chosen_pop_ids_5 <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family_5, 1]
sample_pop_5 <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids_5, ]
sample_pop_male_5 <- sample_pop_5[sample_pop_5[2] == 'M', ]
sample_pop_female_5 <- sample_pop_5[sample_pop_5[2] == 'F', ]
chosen_male_5 <- sample_pop_male_5$id_pop[1:malePerFamily]
chosen_female_5 <- sample_pop_female_5$id_pop
cross_plan_5 <- cross_plan_homo(chosen_male_5, chosen_female_5, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan, cross_plan_5)

if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
  print(nrow(cross_plan))
  stop('cross plan 数量不正确!!!')
}
data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)

# 第 9 代 家系 5 导入基因
k <- 9
print(paste0('第', k, '代交配...'))
ind_offset <- initPopNum + familyNum * malePerFamily * femalePerMale * nProgenyPerLitter * (k - 2)

f <- 5
chosen_family_5 <- all_familys[f]
chosen_pop_ids_5 <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family_5, 1]
sample_pop_5 <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids_5, ]
sample_pop_male_5 <- sample_pop_5[sample_pop_5[2] == 'M', ]
sample_pop_female_5 <- sample_pop_5[sample_pop_5[2] == 'F', ]
chosen_male_5 <- sample_pop_male_5$id_pop[1:malePerFamily]
chosen_female_5 <- sample_pop_female_5$id_pop
# 编辑公猪和母猪各 edit_male_per_family 头
edit_male_ind <- sapply(chosen_male_5[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_male_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)
edit_female_ind <- sapply(chosen_female_5[1:edit_male_per_family], function(x) as.numeric(x) - ind_offset)
data_list[[k]]$pop <- editGenome(data_list[[k]]$pop, ind=edit_female_ind,
                                 chr=target_gene_chr[[f]], segSites=target_gene_site[[f]], allele=1)

chosen_males_1234 <- list()
for (f in 1:4){
  chosen_family <- all_familys[f]
  chosen_pop_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family == chosen_family, 1]
  sample_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% chosen_pop_ids, ]
  sample_pop_male <- sample_pop[sample_pop[2] == 'M', ]
  chosen_male <- sample_pop_male$id_pop[1:malePerFamily]
  chosen_males_1234 <- c(chosen_males_1234, chosen_male)
}
sample_1234_ids <- data_list[[k]]$pedi[data_list[[k]]$pedi$family != all_familys[1], 1]
sample_1234_pop <- data_list[[k]]$data_pop[data_list[[k]]$data_pop$id_pop %in% sample_1234_ids, ]
female_1234_ids <- sample_1234_pop[sample_1234_pop[2] == 'F', 1]
cross_plan_1234_5 <- cross_plan_homo(chosen_males_1234, chosen_female_5, data_list[[k]]$pedi)
cross_plan_5_1234 <- cross_plan_homo(chosen_male_5, female_1234_ids, data_list[[k]]$pedi)
cross_plan <- rbind(cross_plan_1234_5, cross_plan_5_1234)

if (nrow(cross_plan) != familyNum * malePerFamily * femalePerMale){
  print(nrow(cross_plan))
  stop('cross plan 数量不正确!!!')
}
data_list[[k+1]] <- get_progency(data_list[[k]]$pop, cross_plan)

# 第 9 代开始普通杂交
for (k in 9:(maxGeneration - 1)){
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


