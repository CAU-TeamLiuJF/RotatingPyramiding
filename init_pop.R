
print('初始化群体...')
founderPop <- runMacs(nInd=initPopNum, nChr=nChr, segSites=segSites)
SP <- SimParam$new(founderPop)
# 表型、基因型、性别
SP$addTraitA(nQtlPerChr=nQtlPerChr,mean=mean,var=var,corA=corA,name=names)
SP$addSnpChip(nSnpPerChr=nSnpPerChr)
SP$setSexes("yes_sys")
qtl_map <- getQtlMap(trait=1, simParam=SP)
all_qtl_pos <- t(as.data.frame(qtl_map$id))
# 生成 targetGeneNum 个不属于 qtl_map 内的目标基因, 每个基因在不同染色体上
target_gene_chr <- list()
target_gene_site <- list()
target_gene_pos <- list()
while (length(target_gene_chr) < targetGeneNum){
  chosen_chr <- sample(1:nChr, 1)
  chosen_site <- sample(1:segSites, 1)
  chosen_pos <- paste0(chosen_chr, '_', chosen_site)
  if (chosen_pos %in% all_qtl_pos || chosen_chr %in% target_gene_chr){
    next
  }
  target_gene_chr <- c(target_gene_chr, chosen_chr)
  target_gene_site <- c(target_gene_site, chosen_site)
  target_gene_pos <- c(target_gene_pos, chosen_pos)
}
print(unlist(target_gene_pos))
# 产生第 0 代
pop0 <- newPop(founderPop)
pop0 <- setPheno(pop0, h2=h2)
# 将群体的目标基因位点全部设置为“0”，可用pullMarkerGeno(pop, markers=c(target_gene_pos))进行基因型检查
for (k in 1:targetGeneNum){
  pop0 <- editGenome(pop0, ind=1:initPopNum, chr=target_gene_chr[[k]], segSites=target_gene_site[[k]], allele=0)
}
data_0 <- get_pop_data(pop0)
data_0$pedi$family <- data_0$pedi$id
data_0$pedi$mother <- data_0$pedi$id
data_0$pedi$father <- data_0$pedi$id
# 保存初始环境
save(founderPop, SP, pop0, data_0, qtl_map, targetGeneNum,
     target_gene_chr, target_gene_site, target_gene_pos, maxGeneration, file=init_data_file)
