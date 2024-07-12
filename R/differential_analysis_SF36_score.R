# SF.36 Score Analysis
bs.sub <- bs.filtered[,which(targets$Group=="LC")]
targets.sub <- targets[which(targets$Group=="LC"),]
pData(bs.sub)$SF.36_Score <- targets.sub$SF.36_Score
dss.formula <- formula(~ SF.36_Score + Age + Sex + BMI + type1 + type2 + type3 + type4 + type5 + PC1 + PC2)
design.df <- pData(bs.sub)
dml.fit <- DMLfit.multiFactor(bs.sub, design = design.df, smoothing = TRUE,  smoothing.span = 500, formula = dss.formula)
test.result <- DMLtest.multiFactor(dml.fit, coef = 2)
median(qchisq(1-test.result$pvals,1))/qchisq(0.5,1)
# 0.4675222
ss <- fdrtool(test.result$stat, statistic="normal", plot=F)
#median(qchisq(1-ss$pval,1))/qchisq(0.5,1)
# 1.038744
test.result$p.from.ss <- ss$pval
test.result$lfdr <- ss$lfdr
save(bs.sub, test.result, file = "Results/SF_Score/dssResults.SF36_Score.smooth500.rdata")
