# PASC whole-blood single-cell DNA methylation analysis (CpG-only)
# Samples S1-S6: controls S1/S5/S6, cases (PASC) S2/S3/S4.
# Builds the CG/CH bin matrices, runs QC, dimensionality reduction, clustering,
# doublet removal, cell-type annotation, and Case-vs-Control pseudobulk DMRs.

# load packages
library(amethyst)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(cluster)
library(limma)
library(readxl)
library(GenomicRanges)

# set directories
resultsDir <- "/media/data/andy/WGBS/Datasets/Albany/scMethylSeq/results"
sampleDir <- file.path(resultsDir, "samples")
outDir <- file.path(resultsDir, "clustering")
figDir <- file.path(outDir, "figures")
resDir <- file.path(outDir, "results")
ckptDir <- file.path(outDir, "checkpoints")
dir.create(figDir, recursive = TRUE, showWarnings = FALSE)
dir.create(resDir, recursive = TRUE, showWarnings = FALSE)
dir.create(ckptDir, recursive = TRUE, showWarnings = FALSE)

# set cases
controls <- c("S1", "S5", "S6")
cases    <- c("S2", "S3", "S4")
stdChrs  <- paste0("chr", c(1:22, "X"))
diagnosisOf <- function(s) ifelse(s %in% cases, "Case (PASC)", "Control")

set.seed(714)
options(future.globals.maxSize = 64 * 1024^3)


###############################
# Matrix generation + cell QC
###############################

obj <- createScaleObject(directory = sampleDir, genomeMatrices = list("CG.score", "CH"))

meta <- as.data.frame(obj@metadata)
qcVars <- c("cov", "cg_cov", "mcg_pct", "ch_cov", "mch_pct", "ch_high_reads_percent")

dir.create(file.path(figDir, "qc"), recursive = TRUE, showWarnings = FALSE)
pdf(file.path(figDir, "qc", "qc_distributions_prefilter.pdf"), width = 12, height = 9)
for (v in qcVars) {
  p <- ggplot(meta, aes(x = .data[[v]])) +
    geom_histogram(bins = 80, fill = "steelblue", color = "white") +
    facet_wrap(~sample, scales = "free_y") +
    theme_bw() + labs(title = paste("Distribution of", v), x = v, y = "Cells")
  if (v %in% c("cov", "cg_cov", "ch_cov")) p <- p + scale_x_log10()
  print(p)
}
print(ggplot(meta, aes(x = cov, y = mcg_pct, color = sample)) +
    geom_point(size = 0.3, alpha = 0.3) + scale_x_log10() + theme_bw() +
    labs(title = "Total coverage vs mCG%"))
print(ggplot(meta, aes(x = cov, y = ch_high_reads_percent, color = sample)) +
    geom_point(size = 0.3, alpha = 0.3) + scale_x_log10() + theme_bw() +
    labs(title = "Total coverage vs CH-high-reads%"))
dev.off()

# Percentile-based bounds, with an absolute coverage floor so cells have enough
absCovMin <- 50000
covLo <- max(quantile(meta$cov, 0.01, na.rm = TRUE), absCovMin)
covHi <- quantile(meta$cov, 0.995, na.rm = TRUE)
mcgLo <- quantile(meta$mcg_pct, 0.005, na.rm = TRUE)
mcgHi <- quantile(meta$mcg_pct, 0.995, na.rm = TRUE)
chHi <- quantile(meta$ch_high_reads_percent, 0.995, na.rm = TRUE)

keep <- meta$cov >= covLo & meta$cov <= covHi & meta$mcg_pct >= mcgLo & meta$mcg_pct <= mcgHi & meta$ch_high_reads_percent <= chHi
keepIds <- rownames(meta)[keep]

obj@metadata <- obj@metadata[keepIds, ]
obj@h5paths  <- obj@h5paths[obj@h5paths$barcode %in% keepIds, ]
for (m in names(obj@genomeMatrices)) {
    obj@genomeMatrices[[m]] <- obj@genomeMatrices[[m]][, intersect(colnames(obj@genomeMatrices[[m]]), keepIds)]
}

metaF <- as.data.frame(obj@metadata)
metaF$diagnosis <- diagnosisOf(metaF$sample)

pdf(file.path(figDir, "qc", "qc_distributions_postfilter.pdf"), width = 12, height = 9)
for (v in qcVars) {
  p <- ggplot(metaF, aes(x = .data[[v]])) +
    geom_histogram(bins = 80, fill = "darkgreen", color = "white") +
    facet_wrap(~sample, scales = "free_y") +
    theme_bw() + labs(title = paste("Post-filter:", v), x = v, y = "Cells")
  if (v %in% c("cov", "cg_cov", "ch_cov")) p <- p + scale_x_log10()
  print(p)
}
print(ggplot(metaF, aes(x = sample, fill = diagnosis)) +
    geom_bar() + theme_bw() +
    scale_fill_manual(values = c("Case (PASC)" = "#E05C5C", "Control" = "#4A90D9")) +
    labs(title = "Cells per sample after QC filtering", y = "Cells"))
dev.off()

write.csv(metaF, file.path(resDir, "cell_metadata_filtered.csv"))


###############################
# Bin filtering (CpG-only)
###############################

# CpG-only analysis, removing CH matrix just in case, ya know?
obj@genomeMatrices[["CH"]] <- NULL

mat <- obj@genomeMatrices[["CG.score"]]
binChrs <- sub("_[0-9]+_[0-9]+$", "", rownames(mat))
obj@genomeMatrices[["CG.score"]] <- mat[binChrs %in% stdChrs, ]

# Keep bins with coverage (nonzero CG.score) in at least half the cells.
cgThresh <- 0.50
mat <- as.matrix(obj@genomeMatrices[["CG.score"]])
nzFrac <- rowMeans(mat != 0)
rm(mat)

pdf(file.path(figDir, "bin_nonzero_fraction.pdf"), width = 9, height = 4)
print(ggplot(data.frame(frac = nzFrac), aes(x = frac)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    geom_vline(xintercept = cgThresh, linetype = "dashed", color = "red") +
    theme_bw() + labs(title = "CG.score: fraction of cells with nonzero value per bin", x = "Fraction nonzero", y = "Bins"))
dev.off()

obj@genomeMatrices[["CG.score"]] <- obj@genomeMatrices[["CG.score"]][nzFrac >= cgThresh, ]


###############################
# Dimensionality reduction
###############################

dimsCg <- max(unname(dimEstimate(obj, genomeMatrices = c("CG.score"), dims = c(50), threshold = 0.95)["CG.score"]), 10)

obj@reductions[["irlba"]] <- runIrlba(obj, genomeMatrices = c("CG.score"), dims = c(dimsCg), replaceNA = c(0))

# Check whether reduction dimensions track sequencing depth and regress it out if so.
irlbaDf <- obj@reductions[["irlba"]]
meta <- as.data.frame(obj@metadata)
common <- intersect(rownames(irlbaDf), rownames(meta))
covVec <- log(meta[common, "cov"])
corWithCov <- sapply(irlbaDf[common, ], function(x) cor(x, covVec, method = "spearman"))

pdf(file.path(figDir, "dimreduction_diagnostics.pdf"), width = 10, height = 6)
dfCor <- data.frame(dim = factor(names(corWithCov), levels = names(corWithCov)), cor = as.numeric(corWithCov))
print(ggplot(dfCor, aes(x = dim, y = cor)) +
    geom_col(fill = "steelblue") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", color = "red") +
    labs(title = "IRLBA dim correlation with log(cov)", y = "Spearman r"))
plotDf <- irlbaDf[common, 1:2]
colnames(plotDf) <- c("DIM1", "DIM2")
plotDf$sample <- meta[common, "sample"]
plotDf$cov    <- meta[common, "cov"]
print(ggplot(plotDf, aes(DIM1, DIM2, color = sample)) + geom_point(size = 0.3, alpha = 0.5) + theme_bw() + labs(title = "IRLBA DIM1 vs DIM2 by sample"))
print(ggplot(plotDf, aes(DIM1, DIM2, color = log10(cov))) + geom_point(size = 0.3, alpha = 0.5) + theme_bw() + scale_color_viridis_c() + labs(title = "IRLBA DIM1 vs DIM2 by log10(cov)"))
dev.off()

if (sum(abs(corWithCov) > 0.3) >= 3) {
  obj@reductions[["irlba_gam"]] <- regressCovBias(obj, reduction = "irlba", method = "gam")
  reductionFinal <- "irlba_gam"
} else {
  reductionFinal <- "irlba"
}


###############################
# Clustering + UMAP
###############################

red <- as.matrix(obj@reductions[["irlba"]])
silIdx  <- sample(seq_len(nrow(red)), min(5000, nrow(red)))
silDist <- dist(red[silIdx, ])

kValues <- c(10, 15, 20, 30, 50, 75, 100)
kSweep <- data.frame()
clusterAssignments <- list()
for (k in kValues) {
  tmp <- runCluster(obj, k = k, reduction = "irlba", method = "louvain", colname = "cluster_id")
  cl <- tmp@metadata[["cluster_id"]]
  names(cl) <- rownames(tmp@metadata)
  clusterAssignments[[as.character(k)]] <- cl
  sizes <- table(cl)
  sil <- silhouette(as.integer(factor(cl[silIdx])), silDist)
  kSweep <- rbind(kSweep, data.frame(k = k, n_clusters = length(unique(cl)),
      min_size = min(sizes), median_size = median(sizes), max_size = max(sizes),
      mean_silhouette = mean(sil[, "sil_width"])))
}

# Pick the k giving the best silhouette
minSizeThresh <- 0.005 * nrow(obj@metadata)
candidates <- kSweep %>% filter(n_clusters >= 4, n_clusters <= 25, min_size >= minSizeThresh)
if (nrow(candidates) == 0) candidates <- kSweep %>% filter(n_clusters >= 4, n_clusters <= 25)
if (nrow(candidates) == 0) candidates <- kSweep
chosenK <- candidates$k[which.max(candidates$mean_silhouette)]
obj@metadata[["cluster_id"]] <- clusterAssignments[[as.character(chosenK)]][rownames(obj@metadata)]

umapNb <- 30
umapDist <- 0.1
umapCoords <- runUmap(obj, neighbors = umapNb, dist = umapDist, method = "euclidean", reduction = "irlba")
rownames(umapCoords) <- rownames(obj@reductions[["irlba"]])
obj@reductions[["umap"]] <- umapCoords

finalDf <- as.data.frame(umapCoords)
colnames(finalDf) <- c("UMAP1", "UMAP2")
finalDf$cluster_id <- as.factor(obj@metadata[rownames(finalDf), "cluster_id"])
finalDf$sample <- obj@metadata[rownames(finalDf), "sample"]
finalDf$diagnosis <- diagnosisOf(finalDf$sample)

pdf(file.path(figDir, "umap.pdf"), width = 12, height = 9)
print(ggplot(finalDf, aes(UMAP1, UMAP2, color = cluster_id)) +
    geom_point(size = 0.4, alpha = 0.6) + theme_bw() +
    labs(title = sprintf("UMAP by cluster (k=%d)", chosenK)))
print(ggplot(finalDf, aes(UMAP1, UMAP2, color = sample)) +
    geom_point(size = 0.3, alpha = 0.5) + theme_bw() + labs(title = "UMAP by sample"))
print(ggplot(finalDf, aes(UMAP1, UMAP2, color = diagnosis)) +
    geom_point(size = 0.3, alpha = 0.5) + theme_bw() +
    scale_color_manual(values = c("Case (PASC)" = "#E05C5C", "Control" = "#4A90D9")) +
    labs(title = "UMAP by diagnosis"))
print(ggplot(finalDf, aes(UMAP1, UMAP2, color = cluster_id)) +
    geom_point(size = 0.2, alpha = 0.6) + facet_wrap(~sample) + theme_bw() +
    labs(title = "UMAP by cluster, faceted by sample"))
dev.off()


###############################
# Doublet detection & removal
###############################

dbobj <- makeDoubletObject(obj, simFraction = 0.25, threads = 32, genomeMatrices = c("CG.score"))
dbobj@reductions[["irlba"]] <- runIrlba(dbobj, genomeMatrices = c("CG.score"),dims = c(dimsCg), replaceNA = c(0))
dbModel <- buildDoubletModel(dbobj, method = "rf", reduction = "irlba")
dbobj <- predictDoubletScores(dbobj, model = dbModel$model, reduction = "irlba")
obj <- addDoubletScores(obj, dbobj)

meta <- as.data.frame(obj@metadata)
clusterSummary <- meta %>%
    group_by(cluster_id) %>%
    summarise(n = n(), mean_doublet_score = mean(doublet_score), pct_over_0.5 = mean(doublet_score > 0.5) * 100) %>%
    arrange(desc(mean_doublet_score))

pdf(file.path(figDir, "doublet_diagnostics.pdf"), width = 10, height = 6)
print(ggplot(meta, aes(x = doublet_score)) +
    geom_histogram(bins = 60, fill = "darkorange", color = "white") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
    theme_bw() + labs(title = "Doublet score distribution (real cells)"))
print(ggplot(meta, aes(x = factor(cluster_id), y = doublet_score)) +
    geom_boxplot(outlier.size = 0.3) + theme_bw() +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    labs(title = "Doublet score by cluster", x = "Cluster", y = "Doublet score"))
udf <- as.data.frame(obj@reductions[["umap"]])
colnames(udf) <- c("UMAP1", "UMAP2")
udf$doublet_score <- meta[rownames(udf), "doublet_score"]
print(ggplot(udf, aes(UMAP1, UMAP2, color = doublet_score)) +
    geom_point(size = 0.3, alpha = 0.6) + scale_color_viridis_c() + theme_bw() +
    labs(title = "UMAP colored by doublet score"))
dev.off()

# Flag individual cells over the RF boundary, plus any cluster where most cells
# score as doublets.
flaggedClusters <- clusterSummary$cluster_id[clusterSummary$pct_over_0.5 > 50]
isDoublet <- meta$doublet_score > 0.5 | meta$cluster_id %in% flaggedClusters
keepIds <- rownames(meta)[!isDoublet]

obj@metadata <- obj@metadata[keepIds, ]
obj@h5paths  <- obj@h5paths[obj@h5paths$barcode %in% keepIds, ]
for (m in names(obj@genomeMatrices)) {
  obj@genomeMatrices[[m]] <- obj@genomeMatrices[[m]][, intersect(colnames(obj@genomeMatrices[[m]]), keepIds)]
}
for (r in names(obj@reductions)) {
  obj@reductions[[r]] <- obj@reductions[[r]][intersect(rownames(obj@reductions[[r]]), keepIds), ]
}

# Re-cluster and re-embed the surviving cells.
obj <- runCluster(obj, k = chosenK, reduction = "irlba", method = "louvain", colname = "cluster_id")
umapCoords <- runUmap(obj, neighbors = umapNb, dist = umapDist, method = "euclidean", reduction = "irlba")
rownames(umapCoords) <- rownames(obj@reductions[["irlba"]])
obj@reductions[["umap"]] <- umapCoords

finalDf <- as.data.frame(umapCoords)
colnames(finalDf) <- c("UMAP1", "UMAP2")
finalDf$cluster_id <- as.factor(obj@metadata[rownames(finalDf), "cluster_id"])
finalDf$sample <- obj@metadata[rownames(finalDf), "sample"]
finalDf$diagnosis <- diagnosisOf(finalDf$sample)

pdf(file.path(figDir, "umap_postdoublet.pdf"), width = 12, height = 9)
print(ggplot(finalDf, aes(UMAP1, UMAP2, color = cluster_id)) +
    geom_point(size = 0.4, alpha = 0.6) + theme_bw() +
    labs(title = sprintf("Post-doublet-removal UMAP by cluster (k=%d)", chosenK)))
print(ggplot(finalDf, aes(UMAP1, UMAP2, color = sample)) +
    geom_point(size = 0.3, alpha = 0.5) + theme_bw() + labs(title = "Post-doublet-removal UMAP by sample"))
print(ggplot(finalDf, aes(UMAP1, UMAP2, color = diagnosis)) +
    geom_point(size = 0.3, alpha = 0.5) + theme_bw() +
    scale_color_manual(values = c("Case (PASC)" = "#E05C5C", "Control" = "#4A90D9")) +
    labs(title = "Post-doublet-removal UMAP by diagnosis"))
dev.off()


###############################
# Cell-type annotation (FlowSorted.Blood.EPIC reference)
# Honestly not sure if this is the best option
# But gotta play with what you've got
###############################

library(FlowSorted.Blood.EPIC)
library(minfi)
library(rtracklayer)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

epicRg <- libraryDataGet("FlowSorted.Blood.EPIC")
epicRg <- epicRg[, epicRg$CellType != "MIX"]
epicBeta <- getBeta(preprocessNoob(epicRg, verbose = FALSE))

cellTypesRef <- sort(unique(epicRg$CellType))
refByCt <- sapply(cellTypesRef, function(ct)
    rowMeans(epicBeta[, epicRg$CellType == ct, drop = FALSE], na.rm = TRUE))

epicAnno <- as.data.frame(getAnnotation(epicRg))
probesInBoth <- intersect(rownames(refByCt), rownames(epicAnno))
refByCt <- refByCt[probesInBoth,]
epicAnno <- epicAnno[probesInBoth,]

# Lift EPIC probes hg19 -> hg38 and bin them at 50kb to match the CG.score grid.
chain <- import.chain("/media/data/andy/GEO_Datasets/geo_data/predictions/data/hg19ToHg38.over.chain")
probeGr <- GRanges(seqnames = epicAnno$chr, ranges = IRanges(start = epicAnno$pos, width = 1), probe_id = rownames(epicAnno))
lifted <- liftOver(probeGr, chain)
probeGr38 <- unlist(lifted[lengths(lifted) == 1])

binSize <- 50000L
binStart <- (start(probeGr38) %/% binSize) * binSize
binName  <- paste(as.character(seqnames(probeGr38)), binStart, binStart + binSize, sep = "_")

cgmat <- obj@genomeMatrices[["CG.score"]]
validBins <- binName %in% rownames(cgmat)
refSub <- do.call(rbind, lapply(
    split(seq_along(binName[validBins]), binName[validBins]), function(idx)
        colMeans(refByCt[probeGr38$probe_id[validBins][idx], , drop = FALSE], na.rm = TRUE)))

annoBins <- intersect(rownames(refSub), rownames(cgmat))
refSub <- refSub[annoBins, ]
cgmatSub <- as.matrix(cgmat[annoBins,])

meta <- as.data.frame(obj@metadata)
clusters <- sort(unique(meta$cluster_id))
clusterProfiles <- matrix(NA, nrow = length(annoBins), ncol = length(clusters), dimnames = list(annoBins, clusters))
for (cl in clusters) {
  cellsInCl <- intersect(rownames(meta)[meta$cluster_id == cl], colnames(cgmatSub))
  subM <- cgmatSub[, cellsInCl, drop = FALSE]
  rowNz <- rowSums(subM != 0)
  clusterProfiles[, as.character(cl)] <- ifelse(rowNz > 0, rowSums(subM) / rowNz, NA)
}

# Gotta rescale CG.Score (-1, 1) to a beta-like 0..1 to correlate against EPIC betas
clusterBeta <- (clusterProfiles + 1) / 2
corMat <- matrix(NA, nrow = length(clusters), ncol = ncol(refSub), dimnames = list(clusters, colnames(refSub)))
for (cl in clusters) {
    clVals <- clusterBeta[, as.character(cl)]
    valid <- !is.na(clVals) & complete.cases(refSub)
    corMat[as.character(cl), ] <- cor(clVals[valid], refSub[valid, ], method = "pearson")
}

# Assign each cluster to its best-correlating reference type; Neu -> Gran.
remap <- c("Bcell" = "Bcell", "CD4T" = "CD4T", "CD8T" = "CD8T", "Mono" = "Mono", "Neu" = "Gran", "NK" = "NK")
bestTypeRaw <- colnames(refSub)[apply(corMat, 1, which.max)]
clusterAnnot <- data.frame(cluster_id = clusters,
  n_cells = as.integer(table(meta$cluster_id)[as.character(clusters)]),
  best_celltype_raw = bestTypeRaw,
  best_correlation = round(apply(corMat, 1, max), 3),
  cell_type = unname(remap[bestTypeRaw]),
  stringsAsFactors = FALSE)

meta$cell_type <- clusterAnnot$cell_type[match(meta$cluster_id, clusterAnnot$cluster_id)]
obj@metadata$cell_type <- meta[rownames(obj@metadata), "cell_type"]

pdf(file.path(figDir, "cluster_celltype_cor_heatmap.pdf"), width = 10, height = 8)
corDf <- as.data.frame(corMat)
corDf$cluster <- factor(rownames(corDf), levels = as.character(clusters))
corLong <- pivot_longer(corDf, -cluster, names_to = "cell_type", values_to = "correlation")
print(ggplot(corLong, aes(x = cell_type, y = cluster, fill = correlation)) +
    geom_tile() + geom_text(aes(label = round(correlation, 2)), size = 2.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_bw() + labs(title = "Cluster-celltype correlation (FlowSorted.Blood.EPIC)"))
dev.off()

umapDf <- as.data.frame(obj@reductions[["umap"]])
colnames(umapDf)[1:2] <- c("UMAP1", "UMAP2")
umapDf$cell_id <- rownames(umapDf)
umapDf$sample <- meta[umapDf$cell_id, "sample"]
umapDf$cell_type <- meta[umapDf$cell_id, "cell_type"]
umapDf$diagnosis <- diagnosisOf(umapDf$sample)

pdf(file.path(figDir, "umap_celltype.pdf"), width = 14, height = 10)
print(ggplot(umapDf, aes(UMAP1, UMAP2, color = cell_type)) +
    geom_point(size = 0.3, alpha = 0.6) + theme_bw() +
    labs(title = "UMAP by annotated cell type") +
    guides(colour = guide_legend(override.aes = list(size = 3))))
print(ggplot(umapDf, aes(UMAP1, UMAP2, color = cell_type)) +
    geom_point(size = 0.2, alpha = 0.5) + facet_wrap(~sample) + theme_bw() +
    labs(title = "UMAP per sample, colored by cell type") +
    guides(colour = guide_legend(override.aes = list(size = 3))))
set.seed(608)
umapShuf <- umapDf[sample(nrow(umapDf)), ]
print(ggplot(umapShuf, aes(UMAP1, UMAP2, color = diagnosis)) +
    geom_point(size = 0.3, alpha = 0.6) +
    scale_color_manual(values = c("Case (PASC)" = "#E05C5C", "Control" = "#4A90D9")) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_bw() + labs(title = "UMAP by diagnosis", color = "Diagnosis"))
print(ggplot(umapShuf, aes(UMAP1, UMAP2, color = diagnosis)) +
    geom_point(size = 0.2, alpha = 0.5) +
    scale_color_manual(values = c("Case (PASC)" = "#E05C5C", "Control" = "#4A90D9")) +
    facet_wrap(~cell_type) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_bw() + labs(title = "UMAP by diagnosis, per cell type", color = "Diagnosis"))
dev.off()

ctCounts <- as.data.frame(table(meta$cell_type))
colnames(ctCounts) <- c("cell_type", "n_cells")
ctCounts$pct <- round(100 * ctCounts$n_cells / sum(ctCounts$n_cells), 2)

ctSampleDf <- as.data.frame(table(meta$sample, meta$cell_type))
colnames(ctSampleDf) <- c("sample", "cell_type", "n_cells")
ctSamplePct <- ctSampleDf %>%
    group_by(sample) %>% mutate(pct = n_cells / sum(n_cells) * 100) %>% ungroup()

pdf(file.path(figDir, "celltype_composition.pdf"), width = 10, height = 6)
print(ggplot(ctCounts, aes(x = reorder(cell_type, -n_cells), y = pct, fill = cell_type)) +
    geom_bar(stat = "identity") + theme_bw() +
    labs(title = "Cell type composition (all samples)", x = "Cell type", y = "% cells") +
    geom_text(aes(label = paste0(n_cells, " (", pct, "%)")), vjust = -0.3, size = 3))
print(ggplot(ctSamplePct, aes(x = cell_type, y = sample, fill = pct)) +
    geom_tile() + geom_text(aes(label = round(pct, 1)), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue") + theme_bw() +
    labs(title = "Cell type % per sample (scDNAm)"))
dev.off()

# Cross-check proportions against CBC differentials.
cbc <- read_excel(file.path(resultsDir, "PASC_single-cell_metadata_for_andy.xlsx"))
cbc$sample <- paste0("S", seq_len(nrow(cbc)))
cbcLong <- cbc %>%
    transmute(sample,
        Gran = `Neutrophils (%)` + `Immature Granulocytes (%)` + `Eosinophils (%)` + `Basophils (%)`,
        Lymphocytes = `Lymphocytes (%)`,
        Mono = `Monocytes (%)`) %>%
    pivot_longer(-sample, names_to = "category", values_to = "cbc_pct")
scdnamMap <- ctSamplePct %>%
    mutate(category = case_when(
        cell_type %in% c("CD4T", "CD8T", "Bcell", "NK") ~ "Lymphocytes",
        cell_type == "Mono" ~ "Mono",
        cell_type == "Gran" ~ "Gran",
        TRUE ~ "Other")) %>%
      group_by(sample, category) %>% summarise(scdnam_pct = sum(pct), .groups = "drop")
crossCheck <- inner_join(cbcLong, scdnamMap, by = c("sample", "category"))

pdf(file.path(figDir, "cbc_crosscheck.pdf"), width = 8, height = 6)
print(ggplot(crossCheck, aes(x = cbc_pct, y = scdnam_pct, color = category, label = sample)) +
    geom_point(size = 2) + geom_text(vjust = -0.7, size = 3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_bw() + labs(title = "CBC differential % vs scDNAm cell-type %", x = "CBC %", y = "scDNAm %"))
dev.off()

# Cross-check against matched scRNA-seq proportions
rnaseqCsv <- "/media/data/andy/RNAseq/scRNAseq/Datasets/Albany/PASC/pasc_analysis/results/celltype_counts_per_sample.csv"
if (file.exists(rnaseqCsv)) {
    rna <- read.csv(rnaseqCsv) %>%
        mutate(category = case_when(
        cell_type %in% c("CD4+ T cells", "CD8+ T cells", "B cells", "NK cells") ~ "Lymphocytes",
        cell_type == "Monocytes" ~ "Mono",
        cell_type == "Neutrophils" ~ "Gran",
        TRUE ~ "Other")) %>%
    group_by(sample, category) %>% summarise(n = sum(n_cells), .groups = "drop") %>%
    group_by(sample) %>% mutate(rna_pct = 100 * n / sum(n)) %>% ungroup()
  write.csv(rna, file.path(resDir, "scrnaseq_proportions.csv"), row.names = FALSE)
}

write.csv(clusterAnnot, file.path(resDir, "cluster_annotations.csv"), row.names = FALSE)
write.csv(round(corMat, 4), file.path(resDir, "cluster_celltype_correlations.csv"))
write.csv(ctCounts, file.path(resDir, "celltype_overall_proportions.csv"), row.names = FALSE)
write.csv(as.data.frame(ctSamplePct), file.path(resDir, "celltype_per_sample_proportions.csv"), row.names = FALSE)
write.csv(crossCheck, file.path(resDir, "cbc_crosscheck.csv"), row.names = FALSE)
write.csv(as.data.frame(obj@metadata), file.path(resDir, "cell_metadata_annotated.csv"))


###############################
# Pseudobulk DMRs (case vs control, per cell type)
###############################

minCells <- 5
minSamples <- 3
pvalThresh <- 0.05
deltaThresh <- 0.10

meta <- as.data.frame(obj@metadata)
cellTypes <- sort(unique(meta$cell_type))
cgmat <- obj@genomeMatrices[["CG.score"]]

makePseudobulk <- function(cgmat, meta, cellTypeName) {
    cellsCt <- intersect(rownames(meta)[meta$cell_type == cellTypeName], colnames(cgmat))
    pbList <- list()
    for (samp in unique(meta[cellsCt, "sample"])) {
        sampCells <- cellsCt[meta[cellsCt, "sample"] == samp]
        if (length(sampCells) < minCells) next
        subM  <- cgmat[, sampCells, drop = FALSE]
        rowNz <- Matrix::rowSums(subM != 0)
        pbList[[samp]] <- ifelse(rowNz > 0, Matrix::rowSums(subM) / rowNz, NA)
    }
    if (length(pbList) == 0) return(NULL)
    pbMat <- do.call(cbind, pbList)
    rownames(pbMat) <- rownames(cgmat)
    pbMat
    }

allResults <- list()
dmrSummary <- data.frame()
for (ct in cellTypes) {
    pbMat <- makePseudobulk(cgmat, meta, ct)
    if (is.null(pbMat)) next

    ctrlSamps <- intersect(colnames(pbMat), controls)
    caseSamps <- intersect(colnames(pbMat), cases)
    if (length(ctrlSamps) < 2 || length(caseSamps) < 2) next

    nCtrlObs <- rowSums(!is.na(pbMat[, ctrlSamps, drop = FALSE]))
    nCaseObs <- rowSums(!is.na(pbMat[, caseSamps, drop = FALSE]))
    keepBins <- nCtrlObs >= minSamples & nCaseObs >= minSamples
    pbFilt <- pbMat[keepBins, c(ctrlSamps, caseSamps), drop = FALSE]
    group <- factor(ifelse(colnames(pbFilt) %in% cases, "case", "ctrl"), levels = c("ctrl", "case"))

    betaMat <- (pbFilt + 1) / 2
    eps <- 1e-4
    betaClamped <- pmin(pmax(betaMat, eps), 1 - eps)
    logitMat <- log(betaClamped / (1 - betaClamped))
    rowMu <- rowMeans(logitMat, na.rm = TRUE)
    for (j in seq_len(ncol(logitMat))) {
        naRows <- is.na(logitMat[, j])
        logitMat[naRows, j] <- rowMu[naRows]
    }

    fit <- eBayes(lmFit(logitMat, model.matrix(~group)))
    tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")

    meanCtrlBeta <- rowMeans(betaMat[, ctrlSamps, drop = FALSE], na.rm = TRUE)
    meanCaseBeta <- rowMeans(betaMat[, caseSamps, drop = FALSE], na.rm = TRUE)
    deltaBeta <- meanCaseBeta - meanCtrlBeta

    resultsDf <- data.frame(bin = rownames(tt),
        chr = sub("_[0-9]+_[0-9]+$", "", rownames(tt)),
        start = as.integer(sub(".*_([0-9]+)_[0-9]+$", "\\1", rownames(tt))),
        end   = as.integer(sub(".*_[0-9]+_([0-9]+)$", "\\1", rownames(tt))),
        mean_beta_ctrl = round(meanCtrlBeta[rownames(tt)], 4),
        mean_beta_case = round(meanCaseBeta[rownames(tt)], 4),
        delta_beta = round(deltaBeta[rownames(tt)], 4),
        logFC = round(tt$logFC, 4), AveExpr = round(tt$AveExpr, 4), t_stat = round(tt$t, 4),
        P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, cell_type = ct,
        n_ctrl_samples = nCtrlObs[rownames(tt)], n_case_samples = nCaseObs[rownames(tt)],
        stringsAsFactors = FALSE)
    resultsDf <- resultsDf[order(resultsDf$P.Value), ]
    write.csv(resultsDf, file.path(resDir, paste0("all_bins_tested_", ct, ".csv")), row.names = FALSE)

    dmrsNom    <- resultsDf[resultsDf$P.Value < pvalThresh, ]
    dmrsStrict <- dmrsNom[abs(dmrsNom$delta_beta) >= deltaThresh, ]
    if (nrow(dmrsNom) > 0)
        write.csv(dmrsNom, file.path(resDir, paste0("dmrs_", ct, "_nominal.csv")), row.names = FALSE)
    if (nrow(dmrsStrict) > 0)
        write.csv(dmrsStrict, file.path(resDir, paste0("dmrs_", ct, "_strict.csv")), row.names = FALSE)

    allResults[[ct]] <- resultsDf
    dmrSummary <- rbind(dmrSummary, data.frame(cell_type = ct, bins_tested = nrow(resultsDf),
      dmrs_nominal = nrow(dmrsNom), dmrs_strict = nrow(dmrsStrict),
      dmrs_bonferroni = sum(resultsDf$adj.P.Val < pvalThresh),
      min_pval = min(resultsDf$P.Value), min_adj_pval = min(resultsDf$adj.P.Val),
      stringsAsFactors = FALSE))

    vdf <- resultsDf
    vdf$sig <- ifelse(vdf$P.Value < pvalThresh & abs(vdf$delta_beta) >= deltaThresh,
        "DMR (p<0.05 & |delta_beta|>=0.10)",
        ifelse(vdf$P.Value < pvalThresh, "p<0.05", "NS"))
    sigCols <- c("DMR (p<0.05 & |delta_beta|>=0.10)" = "red", "p<0.05" = "orange", "NS" = "grey70")
    pdf(file.path(figDir, paste0("volcano_", ct, ".pdf")), width = 8, height = 6)
    print(ggplot(vdf, aes(x = delta_beta, y = -log10(P.Value), color = sig)) +
        geom_point(size = 0.5, alpha = 0.6) + scale_color_manual(values = sigCols) +
        geom_hline(yintercept = -log10(pvalThresh), linetype = "dashed", color = "grey40") +
        geom_vline(xintercept = c(-deltaThresh, deltaThresh), linetype = "dashed", color = "grey40") +
        theme_bw() + labs(title = paste("Volcano -", ct), x = "Delta beta (case - control)", y = "-log10(P)", color = ""))
  dev.off()

    vdf$avg_beta <- (vdf$mean_beta_ctrl + vdf$mean_beta_case) / 2
    pdf(file.path(figDir, paste0("MA_", ct, ".pdf")), width = 8, height = 6)
    print(ggplot(vdf, aes(x = avg_beta, y = delta_beta, color = sig)) +
        geom_point(size = 0.5, alpha = 0.6) + scale_color_manual(values = sigCols) +
        geom_hline(yintercept = 0, color = "black") +
        theme_bw() + labs(title = paste("MA plot -", ct), x = "Average beta", y = "Delta beta (case - control)", color = ""))
    dev.off()
}

if (length(allResults) > 0) {
    combined <- do.call(rbind, allResults)
    write.csv(combined, file.path(resDir, "all_bins_tested_all_celltypes.csv"), row.names = FALSE)
    write.csv(combined[combined$P.Value < pvalThresh, ],
    file.path(resDir, "all_dmrs_nominal_combined.csv"), row.names = FALSE)
}
write.csv(dmrSummary, file.path(resDir, "dmr_summary.csv"), row.names = FALSE)

if (nrow(dmrSummary) > 0) {
    pdf(file.path(figDir, "dmr_counts.pdf"), width = 8, height = 5)
    print(ggplot(dmrSummary, aes(x = cell_type, y = dmrs_nominal, fill = cell_type)) +
        geom_bar(stat = "identity") + geom_text(aes(label = dmrs_nominal), vjust = -0.3) +
          theme_bw() + labs(title = "Nominal DMR counts per cell type (p<0.05)", x = "Cell type", y = "# DMRs") + theme(legend.position = "none"))
  dev.off()

    pvalDf <- do.call(rbind, lapply(names(allResults), function(ct)
      data.frame(cell_type = ct, P.Value = allResults[[ct]]$P.Value)))
    pdf(file.path(figDir, "pvalue_histograms.pdf"), width = 12, height = 8)
    print(ggplot(pvalDf, aes(x = P.Value)) +
        geom_histogram(bins = 50, fill = "steelblue", color = "white") +
        facet_wrap(~cell_type, scales = "free_y") + theme_bw() +
        labs(title = "P-value distributions by cell type", x = "P-value", y = "Count"))
     dev.off()

    effDf <- do.call(rbind, lapply(names(allResults), function(ct)
        data.frame(cell_type = ct, delta_beta = allResults[[ct]]$delta_beta)))
    pdf(file.path(figDir, "effect_size_distribution.pdf"), width = 12, height = 8)
    print(ggplot(effDf, aes(x = delta_beta, fill = cell_type)) +
        geom_histogram(bins = 60, color = "white") +
        facet_wrap(~cell_type, scales = "free_y") + theme_bw() +
        labs(title = "Delta beta distributions by cell type", x = "Delta beta (case - control)", y = "Count") + theme(legend.position = "none"))
  dev.off()

  for (ct in cellTypes) {
    if (is.null(allResults[[ct]])) next
    dmrsCt <- allResults[[ct]][allResults[[ct]]$P.Value < pvalThresh & abs(allResults[[ct]]$delta_beta) >= deltaThresh, ]
    if (nrow(dmrsCt) < 1) next
    pbMat <- makePseudobulk(cgmat, meta, ct)
    if (is.null(pbMat)) next
    topBins <- intersect(head(dmrsCt$bin, 50), rownames(pbMat))
    if (length(topBins) < 1) next

    colOrder <- c(intersect(controls, colnames(pbMat)), intersect(cases, colnames(pbMat)))
    betaHm <- ((pbMat[topBins, , drop = FALSE] + 1) / 2)[, colOrder, drop = FALSE]
    annCol <- data.frame(Group = ifelse(colOrder %in% cases, "Case", "Control"), row.names = colOrder)
    pdf(file.path(figDir, paste0("heatmap_top_dmrs_", ct, ".pdf")),
        width = 8, height = max(4, min(15, 0.2 * length(topBins) + 3)))
    pheatmap::pheatmap(betaHm, annotation_col = annCol,
        cluster_rows = length(topBins) > 1, cluster_cols = FALSE,
        show_rownames = length(topBins) <= 30,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        main = paste("Top DMRs -", ct, "(p<0.05 & |delta_beta|>=0.10)"))
    dev.off()
  }
}

# save 'er on up
save(obj, file = file.path(ckptDir, "obj_final.rda"))

echo("See you, space cowboy . . .")
