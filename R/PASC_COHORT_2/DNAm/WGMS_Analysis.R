# Whole Blood Whole-Genome Methylation Sequencing Samples
# PASC Cohort 2 versus Various Control Cohorts
# DMRs: qval < 0.05, >5% diff, 5+ CpGs
# Updated: June 23, 2026

# load packages
library(bsseq)
library(DSS)
library(dmrseq)
library(DMRichR)
library(fdrtool)
library(dplyr)
library(ggsci)
library(methylCC)
library(oliveR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(devtools)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
source_url("https://raw.githubusercontent.com/andymadrid/One-Off-Scripts/refs/heads/main/madridPal.R")
source_url("https://raw.githubusercontent.com/andymadrid/One-Off-Scripts/refs/heads/main/plotEmpiricalDistribution2.r")

# Main comparisons of interest:
# PASC 2 vs Albany Healthy controls
# PASC 2 vs Colorado Controls
# PASC 2 vs Israeli Controls

# load pregenerated datasets
load("combined.datasets.pasc.israeli.gse270454.pasc3.rda")
#   29011907 methylation loci
#   615 samples

# filter to standard, autosomes
bs <- chrSelectBSseq(bs,seqnames=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))

# cohort subsets
pasc1 <- which(pData(bs)$Cohort == "PASC 1")
pasc2 <- which(pData(bs)$Cohort == "PASC 2")
pascfu <- which(pData(bs)$Cohort == "PASC f/u")
pcFU <- which(pData(bs)$Cohort == "Acute f/u (no PASC)")
healthy <- which(pData(bs)$Cohort == "Healthy")
colorado.ctrls <- which(pData(bs)$Cohort=="WBC_Control")
israeli.ctrls <- which(pData(bs)$Cohort=="Control")
prepandemic <- grep("Pandemic", pData(bs)$Cohort)

# filter to only groups of interest here
bs <- bs[,c(pasc2, healthy, colorado.ctrls, israeli.ctrls, prepandemic)]

# filter out sex chromosomes and low coverage sites
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov") < 3) == 0)
bs <- bs[loci.idx,]
# 17,930,535 CpGs at 3x

# get new subgroups
pasc2 <- which(pData(bs)$Cohort == "PASC 2")
healthy <- which(pData(bs)$Cohort == "Healthy")
colorado.ctrls <- which(pData(bs)$Cohort=="WBC_Control")
israeli.ctrls <- which(pData(bs)$Cohort=="Control")
prepandemic <- grep("Pandemic", pData(bs)$Cohort)

# remove prepandemic from current analysis
bs.sub <- bs[,-c(prepandemic)]

# PCA
# subset samples
meth.mat <- bsseq::getMeth(bs.sub, type = "raw")
rVars <- meth.mat[order(rowVars(meth.mat), decreasing = T),]
rVars <- rVars[1:round(0.05 * nrow(meth.mat)),]
library(ggfortify)
pc <- prcomp(t(rVars), center = T, scale = T)
pc2 <- as.data.frame(pc$x[,1:10])
pData(bs.sub) <- cbind(pData(bs.sub), pc2)
x <- as.data.frame(pc$x)
var_explained <- pc$sdev^2/sum(pc$sdev^2)
colPal <- madridPal(4)
pData(bs.sub)$Cohort2 <- pData(bs.sub)$Cohort
pData(bs.sub)$Cohort2 <- gsub("Healthy", "A", pData(bs.sub)$Cohort2)
pData(bs.sub)$Cohort2 <- gsub("WBC_Control", "C", pData(bs.sub)$Cohort2)
pData(bs.sub)$Cohort2 <- gsub("Control", "I", pData(bs.sub)$Cohort2)
pData(bs.sub)$Cohort2 <- factor(pData(bs.sub)$Cohort2, levels = c("PASC 2", "A", "C", "I"))
pdf("/media/data/pca_cohort.pdf")
ggplot(data=x,aes(x=PC1,y=PC2,color=factor(pData(bs.sub)$Cohort2))) +
	geom_point(size=2) +
	theme_classic() +
	theme(text=element_text(size=20)) +
	theme(legend.position="right") +
	labs(color = "Cohort") +
	labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
	scale_color_manual(values=colPal)
dev.off()

# estimate cell counts
bs.hg19 <- liftBuild(bs, current = "hg38", new = "hg19")

set.seed(714)
e <- estimatecc(bs.hg19)
cc <- cell_counts(e)
pData(bs) <- cbind(pData(bs), cc)

# Differential Methylation of PASC2 vs All Controls
message("Working on PASC2 vs All Controls")
group <- pData(bs)$Cohort
group[!grepl("PASC", group)] <- "Control"
pData(bs)$Cohort2 <- group
cores <- 20
regions <- dmrseq(bs, testCovariate = "Cohort2", maxPerms = 20, minNumRegion = 5, cutoff = 0.05, BPPARAM= BiocParallel::MulticoreParam(workers = cores), adjustCovariate = c("Age", "Sex", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "PC1", "PC2"))
save(regions, bs, file="dmrseqResults.PASC2vAllControls.rdata")

# Differential Methylation of PASC2 vs Colorado Ctrls
message("Working on PASC2 vs Colorado Controls")
bs.sub <- bs[,c(pasc2, colorado.ctrls)]
pData(bs.sub)$Cohort <- gsub("WBC_Control", "Control", pData(bs.sub)$Cohort)
cores <- 20
regions <- dmrseq(bs.sub, testCovariate = "Cohort", maxPerms = 20, minNumRegion = 5, cutoff = 0.05, BPPARAM= BiocParallel::MulticoreParam(workers = cores), adjustCovariate = c("Age", "Sex", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "PC1", "PC2"))
save(regions, bs.sub, file="dmrseqResults.PASC2vColoradoCtrls.rdata")

# Differential Methylation of PASC2 vs Israeli Ctrls
message("Working on PASC2 vs Israeli Controls")
bs.sub <- bs[,c(pasc2, israeli.ctrls)]
cores <- 20
regions <- dmrseq(bs.sub, testCovariate = "Cohort", maxPerms = 20, minNumRegion = 5, cutoff = 0.05, BPPARAM= BiocParallel::MulticoreParam(workers = cores), adjustCovariate = c("Age", "Sex", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "PC1", "PC2"))
save(regions, bs.sub, file="dmrseqResults.PASC2vIsraeliCtrls.rdata")

# Differential Methylation of PASC2 vs Healthy
message("Working on PASC2 vs Albany Controls")
bs.sub <- bs[,c(pasc2, healthy)]
cores <- 20
regions <- dmrseq(bs.sub, testCovariate = "Cohort", maxPerms = 20, minNumRegion = 5, cutoff = 0.05, BPPARAM= BiocParallel::MulticoreParam(workers = cores), adjustCovariate = c("Age", "Sex", "Bcell", "CD4T", "CD8T", "Gran", "Mono", "PC1", "PC2"))
save(regions, bs.sub, file="dmrseqResults.PASC2vH.rdata")

###########################################
# Find DMRs that are conserved across Ctrl cohorts
###########################################

library(GenomicRanges)
setwd("/media/data/andy/WGBS/Datasets/Albany/Long_COVID/Differential_Analysis/Tertiary_Abalysis")
load("pasc1.h.dmrs.rda")
regions.pacs1.v.h <- regions.pacs1.v.h[which(regions.pacs1.v.h$qval<0.05)]
load("justDMRs.pasc2.v.albany.rda")
regions.al <- regions
load("justDMRs.pasc2.v.is.rda")
regions.is <- regions
load("justDMRs.pasc2.v.co.rda")
regions.co <- regions

regions <- regions.al
o <- as.data.frame(findOverlaps(regions, regions.pacs1.v.h))
o1 <- mcols(regions[o[,1]])$beta
o2 <- mcols(regions.pacs1.v.h[o[,2]])$beta
up <- which(o1 > 0 & o2 > 0)
down <- which(o1 < 0 & o2 < 0)
x <- c(up, down)
regions.pacs1.v.h.sub <- regions.pacs1.v.h[o[x,2]]

regions <- regions.co
o <- as.data.frame(findOverlaps(regions, regions.pacs1.v.h.sub))
o1 <- mcols(regions[o[,1]])$beta
o2 <- mcols(regions.pacs1.v.h.sub[o[,2]])$beta
up <- which(o1 > 0 & o2 > 0)
down <- which(o1 < 0 & o2 < 0)
x <- c(up, down)
regions.pacs1.v.h.sub <- regions.pacs1.v.h.sub[o[x,2]]

regions <- regions.is
o <- as.data.frame(findOverlaps(regions, regions.pacs1.v.h.sub))
o1 <- mcols(regions[o[,1]])$beta
o2 <- mcols(regions.pacs1.v.h.sub[o[,2]])$beta
up <- which(o1 > 0 & o2 > 0)
down <- which(o1 < 0 & o2 < 0)
x <- c(up, down)
regions.pacs1.v.h.sub <- regions.pacs1.v.h.sub[o[x,2]]

regions <- reduce(regions.pacs1.v.h.sub)
save(regions, file = "pasc1.dmrs.found.in.ctrl.cohorts.rda")

############################################
# Figure-Making Time
############################################

### Manhattan plots
makeManhattan <- function(regions, cols, pdf) {
	
	# load packages needed
	library(dplyr)
	library(tidyverse)
	library(ggsci)

	regions <- as.data.frame(regions)
	
	# Order chromosomes accordingly
	regions$seqnames <- factor(regions$seqnames, levels=c(paste0("chr",1:22)))
	regions <- regions[order(regions$seqnames),]
	regions$chr <- regions$seqnames

	# Calculate size of chromosomes
	toPlot2 <- regions %>% 
		group_by(chr) %>% 
		summarise(chr_len=max(end)) %>%

	# Calculate cumulative bases per chromosome
		mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
		dplyr::select(-chr_len) %>%

	# Add data back to dataset
	left_join(regions, ., by=c("chr"="chr")) %>%
	arrange(chr, end) %>%
	mutate( startC=end+tot)

	# Prepare x-axis
	newAxis <- toPlot2 %>% group_by(chr) %>% summarize(center=( max(startC) + min(startC) ) / 2 )

	# Change lFDR direction depending on if hyper or hypomethylated (slow moving honey)
   	toPlot2 <- toPlot2 %>% 
		dplyr::mutate(
        new.lfdr = 
	case_when(beta > 0 ~ -log10(qval),
		beta < 0 ~ log10(qval)))

	# Prepare color palette for chromosomes
	colPal <- cols
        colPal <- rep(cols, times = 11)

	# East side plot it out…north side plot it out
	pdf(pdf, width = 12, height = 5)
	man <- ggplot(toPlot2, aes(x=startC, y=new.lfdr)) +
		geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
		scale_color_manual(values = colPal) +
		scale_x_continuous( label = newAxis$chr, breaks= newAxis$center ) +
		scale_y_continuous(expand = c(0, 0) ) + 
		xlab("Chromosome") + 
		ylab("-log10(qval)") + 
		ylim(-2,2) +
		theme_bw() +
		theme( 
     			legend.position="none",
     			panel.border = element_blank(),
      			panel.grid.major.x = element_blank(),
      			panel.grid.minor.x = element_blank(),
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
		geom_hline(yintercept = -log10(0.05), color = "red") +
		geom_hline(yintercept = log10(0.05), color = "red")
	print(man)
	dev.off()
}
setwd("/media/data/WGBS/Datasets/Albany/Long_COVID/Differential_Analysis/Tertiary_Abalysis")
load("dmrseqResults.PASC2vColoradoCtrls.rdata")
makeManhattan(regions, c("cadetblue4", "cadetblue2"), "/media/data/manhattan_colorado.pdf")
load("dmrseqResults.PASC2vIsraeliCtrls.rdata")
makeManhattan(regions, c("darkgoldenrod4", "darkgoldenrod2"), "/media/data/manhattan_israeli.pdf")
load("dmrseqResults.PASC2vH.rdata")
makeManhattan(regions, c("darkorange4", "darkorange2"), "/media/data/manhattan_albany.pdf")

############################################ 
### Volcano plots of PASC 2 vs Various Control Cohorts
############################################ 
library(GenomicRanges)
load("pasc1.h.dmrs.rda")
regions.pacs1.v.h <- regions.pacs1.v.h[which(regions.pacs1.v.h$qval<0.05)]

makeVolcano <- function(regions, pdf) {
    regions.sig <- regions[which(regions$qval<0.05)]
    o <- as.data.frame(findOverlaps(regions.sig, regions.pacs1.v.h))
    o1 <- mcols(regions.sig[o[,1]])$beta
    o2 <- mcols(regions.pacs1.v.h[o[,2]])$beta
    down <- which(o1 < 0 & o2 < 0)
    up <- which(o1 > 0 & o2 > 0)
    mcols(regions)$col <- "lightgrey"
    mcols(regions)$size <- 0.5
    mcols(regions.sig)$col <- "lightgrey"
    mcols(regions.sig)$size <- 0.5
    if (length(up) > 0) {
        mcols(regions.sig[o[up,1]])$col <- "firebrick2"
        mcols(regions.sig[o[up,1]])$size <- 2
     }
    if (length(down) > 0) {
        mcols(regions.sig[o[down,1]])$col <- "blue"
        mcols(regions.sig[o[down,1]])$size <- 2
    }
    regions2 <- c(regions, regions.sig)
    regions2 <- regions2[order(mcols(regions2)$col, decreasing=T)]
    pdf(pdf)
         plot(regions2$beta, -log10(regions2$qval),
             pch = 21, bg = mcols(regions2)$col,
             xlab = "Effect Size (Beta)",
             ylab = "-log10(qval)", cex = mcols(regions2)$size,
             main = "PASC 2 vs Albany Controls")
         abline(v = 0, lty = 2, lwd = 3)
         abline(h = -log10(0.05), lty = 2, lwd = 3)
     dev.off()
}
load("dmrseqResults.PASC2vColoradoCtrls.rdata")
makeVolcano(regions, "/media/data/volcanoCO.pdf")
load("dmrseqResults.PASC2vIsraeliCtrls.rdata")
makeVolcano(regions, "/media/data/volcanoIS.pdf")
load("dmrseqResults.PASC2vH.rdata")
makeVolcano(regions, "/media/data/volcanoAlbany.pdf")

############################################ 
# Pie charts of genomic distribution
############################################ 

library(ggplot2)
library(devtools)
library(scales)
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
source_url("https://raw.githubusercontent.com/andymadrid/One-Off-Scripts/refs/heads/main/madridPal.R")
pal <- madridPal(7)
pal <- pal[-6]

makePie <- function(regions, pdf) {
    regions.sig <- regions[which(regions$qval < 0.05)]
    peaks <- annotatePeak(regions.sig, tssRegion=c(-2000, 0), level = "gene", annoDb="org.Hs.eg.db", TxDb=txdb)
    peaks <- as.data.frame(peaks)
    toUse <- peaks
    totalNum <- nrow(toUse)
    intronPer <- 100*(length(grep("Intron", peaks$annotation))/totalNum)
    exonPer <- 100*(length(grep("Exon", peaks$annotation))/totalNum)
    promoterPer <- 100*(length(grep("Promoter", peaks$annotation))/totalNum)
    utr5Per <- 100*(length(grep("5' UTR", peaks$annotation))/totalNum)
    utr3Per <- 100*(length(grep("3' UTR", peaks$annotation))/totalNum)
    downstreamPer <- 100*(length(grep("Downstream", peaks$annotation))/totalNum)
    interPer <- 100*(length(grep("Intergenic", peaks$annotation))/totalNum)
    if(intronPer > 1) { intronPer <- round(intronPer)}
    if(exonPer > 1) { exonPer <- round(exonPer)}
    if(promoterPer > 1) { promoterPer <- round(promoterPer)}
    if(utr5Per > 1) { utr5Per <- round(utr5Per)}
    if(utr3Per > 1) { utr3Per <- round(utr3Per)}
    if(downstreamPer > 1) { downstreamPer <- round(downstreamPer)}
    if(interPer > 1) { interPer <- round(interPer)}
    if(intronPer < 1 & intronPer > 0) { intronPer <- 1}
    if(exonPer < 1 & exonPer > 0) { exonPer <- 1}
    if(promoterPer < 1 & promoterPer > 0) { promoterPer <- 1}
    if(utr5Per < 1 & utr5Per > 0) { utr5Per <- 1}
    if(utr3Per < 1 & utr3Per > 0) { utr3Per <- 1}
    if(downstreamPer < 1 & downstreamPer > 0) { downstreamPer <- 1}
    if(interPer < 1 & interPer > 0) { interPer <- 1}
    x <- c("Promoter" = promoterPer, "5' UTR" = utr5Per, "Exon" = exonPer, "Intron" = intronPer, "3' UTR" = utr3Per, "Intergenic" = interPer)
    diff <- sum(x)-100
    x[6] <- x[6]-diff
    data <- data.frame(group = names(x), value = x)
    data$show <- data$value
    data$group <- factor(data$group, levels = c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Intergenic"))
    pdf(pdf)
    p <- ggplot(data, aes(x="", y=value, fill = group)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) + 
        theme_void() + 
        geom_text(aes(label = show), position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values = pal)
    print(p)
    dev.off()
}
load("dmrseqResults.PASC2vColoradoCtrls.rdata")
makePie(regions, "/media/data/pieCO.pdf")
load("dmrseqResults.PASC2vIsraeliCtrls.rdata")
makePie(regions, "/media/data/pieIS.pdf")
load("dmrseqResults.PASC2vH.rdata")
makePie(regions, "/media/data/pieAlbany.pdf")

############################################ 
### Annotate and write out results
############################################ 
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

load("pasc1.h.dmrs.rda")
regions.pacs1.v.h <- regions.pacs1.v.h[which(regions.pacs1.v.h$qval<0.05)]

load("dmrseqResults.PASC2vColoradoCtrls.rdata")
regions <- regions[which(regions$qval < 0.05)]
peaks <- annotatePeak(regions, tssRegion=c(-2000, 0), level = "gene", annoDb="org.Hs.eg.db", TxDb=txdb)
peaks <- as.data.frame(peaks)
peaks <- peaks[,-c(12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 25)]
peaks[grep("Distal", peaks$annotation), "SYMBOL"] <- NA
o <- as.data.frame(findOverlaps(regions, regions.pacs1.v.h))
o1 <- mcols(regions[o[,1]])$beta
o2 <- mcols(regions.pacs1.v.h[o[,2]])$beta
up <- which(o1 > 0 & o2 > 0)
down <- which(o1 < 0 & o2 < 0)
mcols(regions)$overlap <- "no"
if (length(up > 0)) {
    peaks[o[up,1], "overlap"] <- "YES"
}
if (length(down > 0)) {
    peaks[o[down,1], "overlap"] <- "YES"
}
write.table(peaks, file = "/media/data/dmrs.colorado.txt", quote = F, row.names = F, sep ='\t')

load("dmrseqResults.PASC2vIsraeliCtrls.rdata")
regions <- regions[which(regions$qval < 0.05)]
peaks <- annotatePeak(regions, tssRegion=c(-2000, 0), level = "gene", annoDb="org.Hs.eg.db", TxDb=txdb)
peaks <- as.data.frame(peaks)
peaks <- peaks[,-c(12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 25)]
peaks[grep("Distal", peaks$annotation), "SYMBOL"] <- NA
o <- as.data.frame(findOverlaps(regions, regions.pacs1.v.h))
o1 <- mcols(regions[o[,1]])$beta
o2 <- mcols(regions.pacs1.v.h[o[,2]])$beta
up <- which(o1 > 0 & o2 > 0)
down <- which(o1 < 0 & o2 < 0)
mcols(regions)$overlap <- "no"
if (length(up > 0)) {
    peaks[o[up,1], "overlap"] <- "YES"
}
if (length(down > 0)) {
    peaks[o[down,1], "overlap"] <- "YES"
}
write.table(peaks, file = "/media/data/dmrs.israeli.txt", quote = F, row.names = F, sep ='\t')

load("dmrseqResults.PASC2vH.rdata")
regions <- regions[which(regions$qval < 0.05)]
peaks <- annotatePeak(regions, tssRegion=c(-2000, 0), level = "gene", annoDb="org.Hs.eg.db", TxDb=txdb)
peaks <- as.data.frame(peaks)
peaks <- peaks[,-c(12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 25)]
peaks[grep("Distal", peaks$annotation), "SYMBOL"] <- NA
o <- as.data.frame(findOverlaps(regions, regions.pacs1.v.h))
o1 <- mcols(regions[o[,1]])$beta
o2 <- mcols(regions.pacs1.v.h[o[,2]])$beta
up <- which(o1 > 0 & o2 > 0)
down <- which(o1 < 0 & o2 < 0)
mcols(regions)$overlap <- "no"
if (length(up > 0)) {
    peaks[o[up,1], "overlap"] <- "YES"
}
if (length(down > 0)) {
    peaks[o[down,1], "overlap"] <- "YES"
}
write.table(peaks, file = "/media/data/dmrs.albany.txt", quote = F, row.names = F, sep ='\t')

load("dmrseqResults.PASC2vPrePandemicHealthy.rdata")
regions <- regions[which(regions$qval < 0.05)]
peaks <- annotatePeak(regions, tssRegion=c(-2000, 0), level = "gene", annoDb="org.Hs.eg.db", TxDb=txdb)
peaks <- as.data.frame(peaks)
peaks <- peaks[,-c(12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 25)]
peaks[grep("Distal", peaks$annotation), "SYMBOL"] <- NA
o <- as.data.frame(findOverlaps(regions, regions.pacs1.v.h))
o1 <- mcols(regions[o[,1]])$beta
o2 <- mcols(regions.pacs1.v.h[o[,2]])$beta
up <- which(o1 > 0 & o2 > 0)
down <- which(o1 < 0 & o2 < 0)
mcols(regions)$overlap <- "no"
if (length(up > 0)) {
    peaks[o[up,1], "overlap"] <- "YES"
}
if (length(down > 0)) {
    peaks[o[down,1], "overlap"] <- "YES"
}
write.table(peaks, file = "/media/data/dmrs.prepan.txt", quote = F, row.names = F, sep ='\t')

echo("See you, space cowboy . . .")
