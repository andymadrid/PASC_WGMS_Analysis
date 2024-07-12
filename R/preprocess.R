# Clean up and set working directory
rm(list=ls())
setwd("/media/data/WGBS/Datasets/Albany/Long_COVID/Differential_Analysis")

# load R packages used for analysis
suppressPackageStartupMessages({
library(DSS)
library(ggsci)
library(dmrseq)
library(ChIPseeker)
library(BiocParallel)
library(DMRichR)
library(clusterProfiler)
library(devtools)
library(wesanderson)
library(ggfortify)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)
library(fdrtool)
library(viridis)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
})
hgAnno <- getAnnot("hg38")

# load in phenotypic data
targets <- read.csv("samplesheet.csv", header = T)
cellCounts <- read.table("cellCounts.txt",header=T)
targets <- cbind(targets,cellCounts)
targets$Group <- rep(c("H","LC"),c(27,108))

# get coverage files
infile <- list.files(pattern="cov$")

# generate matrices
bs <- read.bismark(files = infile,rmZeroCov=TRUE,strandCollapse=T,verbose=T)
x <- colnames(bs)
x <- gsub(".CpG_report.merged_CpG_evidence.cov","",x)
colnames(bs) <- x

# remove sex chromosomes
bs <- chrSelectBSseq(bs,seqnames=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))

# remove samples from cohort 2
bs <- bs[,1:130]
targets <- targets[1:130,]

# filter CpGs based on read coverage
ix <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")<5) == 0)
bs.filtered <- bs[ix,]

# PCA of samples (top 5% most variable CpGs)
meth.mat <- bsseq::getMeth(bs.filtered,type="raw")
rVars <- meth.mat[order(rowVars(meth.mat),decreasing=T),]
rVars <- rVars[1:round(0.05*nrow(meth.mat)),]
library(ggfortify)
pc <- prcomp(t(rVars),center=T,scale=T)
pc2 <- as.data.frame(pc$x[,1:2])
x <- as.data.frame(pc$x)
var_explained <- pc$sdev^2/sum(pc$sdev^2)
#pdf("/media/data/pca.pdf")
#ggplot(data=x,aes(x=PC1,y=PC2,color=targets$Group)) + geom_point(size=-1) + theme_classic() + theme(legend.position="none") + theme(text=element_text(size=20)) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +scale_color_manual(values=c("firebrick3","darkblue")) + geom_text(label=colnames(bs.filtered))
#dev.off()

# add metadata to bsseq object
targets <- cbind(targets,pc2)
sampleCov <- getCoverage(bs.filtered,type="Cov")
covMean <- colMeans(sampleCov)
sampleMeth <- bsseq::getMeth(bs.filtered,type="raw")
methMean <- colMeans(sampleMeth)
targets$Cov <- covMean
targets$Meth <- methMean
pData(bs.filtered)$Condition <- targets$Group
pData(bs.filtered)$Age <- targets$Age
pData(bs.filtered)$Sex <- targets$Sex
pData(bs.filtered)$type1 <- targets$type1
pData(bs.filtered)$type2 <- targets$type2
pData(bs.filtered)$type3 <- targets$type3
pData(bs.filtered)$type4 <- targets$type4
pData(bs.filtered)$type5 <- targets$type5
pData(bs.filtered)$type6 <- targets$type6
pData(bs.filtered)$PC1 <- targets$PC1
pData(bs.filtered)$PC2 <- targets$PC2
pData(bs.filtered)$Cov <- targets$Cov
pData(bs.filtered)$Meth <- targets$Meth
pData(bs.filtered)$BMI <- targets$BMI
pData(bs.filtered)$Race <- targets$Ethnicity
