
# 同质选配
cross_plan_homo <- function(maleList, sortedFemalePop, pedi){
  crossPlan <- data.frame(id_M = numeric(), id_F = numeric(), family = numeric())
  colnames(crossPlan) <- c('id_M', 'id_F', 'family')
  female_chosen_list <- list()
  for (male_id in maleList){
    male_pedi <- pedi[pedi$id == male_id, ]
    male_father <- male_pedi$father
    male_mother <- male_pedi$mother
    male_family <- male_pedi$family

    female_for_male_num <- 0
    female_idx <- 1
    while (female_for_male_num < femalePerMale){
      female_id <- sortedFemalePop[female_idx]
      female_pedi <- pedi[pedi$id == female_id, ]
      female_father <- female_pedi$father
      female_mother <- female_pedi$mother
      # 母猪不能重复选
      if (female_id %in% female_chosen_list){
        female_idx <- female_idx + 1
        next
      }
      # 母猪与公猪不能同属于一窝
      # print("---------------------")
      # print(male_father)
      # print(female_idx)
      # print(head(sortedFemalePop))
      # print(sortedFemalePop[1])
      # print(female_id)
      # # print(female_pedi)
      # # print(female_father)
      # # print(male_mother)
      # # print(female_mother)
      # print("---------------------")
      if (male_father != female_father || male_mother != female_mother){
        new_plan <- data.frame(id_M = male_id, id_F = female_id, family = male_family)
        colnames(new_plan) <- c('id_M', 'id_F', 'family')
        # print(crossPlan)
        crossPlan <- rbind(crossPlan, new_plan)
        # print(new_plan)
        female_for_male_num <- female_for_male_num + 1
        female_chosen_list <- c(female_chosen_list, female_id)
      }
      female_idx <- female_idx + 1
    }
  }
  if (nrow(crossPlan) != (length(maleList) * femalePerMale)){
    stop('cross plan 数量不正确!!!')
  }
  return(crossPlan)
}

# 近交最小选配
cross_plan_inbreed <- function(maleList, sortedFemalePop, pedi, relMatrix){
  crossPlan <- data.frame(id_M = numeric(), id_F = numeric(), family = numeric())
  female_chosen_list <- list()
  for (male_id in maleList){
    male_pedi <- pedi[pedi$id == male_id, ]
    male_father <- male_pedi$father
    male_mother <- male_pedi$mother
    male_family <- male_pedi$family

    female_for_male_num <- 0
    female_idx <- 1
    female_sorted_by_inbreed <- sort_pop_by_inbreed(sortedFemalePop, male_id, relMatrix)
    while (female_for_male_num < femalePerMale){
      female_id <- female_sorted_by_inbreed[female_idx]
      female_pedi <- pedi[pedi$id == female_id, ]
      female_father <- female_pedi$father
      female_mother <- female_pedi$mother
      # 母猪不能重复选
      if (female_id %in% female_chosen_list){
        female_idx <- female_idx + 1
        next
      }
      # 母猪与公猪不能同属于一窝
      if (male_father != female_father || male_mother != female_mother){
        new_plan <- data.frame(id_M = male_id, id_F = female_id, family = male_family)
        crossPlan <- rbind(crossPlan, new_plan)
        female_for_male_num <- female_for_male_num + 1
        female_chosen_list <- c(female_chosen_list, female_id)
      }
      female_idx <- female_idx + 1
    }
  }
  if (nrow(crossPlan) != (length(maleList) * femalePerMale)){
    stop('cross plan 数量不正确!!!')
  }
  return(crossPlan)
}

get_progency <- function (pop, crossPlan){
  cross_plan_tmp <- matrix(c(crossPlan$id_F, crossPlan$id_M), nrow=nrow(crossPlan), ncol=2)
  progency <- makeCross(pop, cross_plan_tmp, nProgeny=nProgenyPerLitter, simParam=SP)
  progency <- setPheno(progency, h2=h2)
  data_progency <- get_pop_data(progency)
  data_progency$pedi$family <- NA
  n <- nrow(crossPlan)
  for (i in 1:n){
    for (j in ((i-1)*nProgenyPerLitter + 1):(i*nProgenyPerLitter)){
      data_progency$pedi$family[j] <- crossPlan$family[i]
    }
  }
  return(data_progency)
}
