# Differential Analysis

# dmrseq
cores <- 50
regions <- dmrseq(bs.filtered, testCovariate = "Condition", maxPerms = 10, minNumRegion = 5, cutoff = 0.05, BPPARAM= BiocParallel::MulticoreParam(workers = cores), adjustCovariate = c("Age","Sex","BMI","type1","type2","type3","type4","type5","PC1","PC2"))
save(regions,bs.filtered,file="Results/Groupwise/dmrseqResults.LCvH.rdata")
