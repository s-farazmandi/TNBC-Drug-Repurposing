# Loading Necessary Libraries & Setting Working Directory
pkgs<- c("stringr", "biomaRt", "DESeq2", "BiocParallel", "clusterProfiler",
         "DOSE", "ReactomePA", "circlize", "RColorBrewer", "corrplot",
         "org.Hs.eg.db", "data.table", "tidyverse", "Hmisc", "AnnotationDbi", "ComplexHeatmap")
lapply(pkgs, require, character.only = TRUE)
setwd("C://Users/Shimbill/Desktop/work/R")



# Generating Count Matrix
library(data.table)
generate_count_matrix <- function(path, pattern) {
  files = list.files(path, pattern, full.names = TRUE, recursive=TRUE, include.dirs=TRUE)
  count_matrix = as.data.frame(do.call(cbind, lapply(files, function(x) fread(x, stringsAsFactors = FALSE))))
  count_matrix <- count_matrix[-c(1:4),]
  rownames(count_matrix) = count_matrix[,1]
  count_matrix = as.data.frame(count_matrix[, seq(4, ncol(count_matrix), 9)])
  return(count_matrix)
}


## raw counts
counts <- generate_count_matrix("tnbc", "\\.rna_seq.augmented_star_gene_counts.tsv$")



# Fetching Metadata
library(tidyverse)
file_names <- list.files("tnbc", "\\.rna_seq.augmented_star_gene_counts.tsv$", full.names = FALSE, recursive = TRUE, include.dirs = FALSE)
file_names <- sub(".*/", "", file_names)
samplesheet <- read.table("tnbc/meta.tsv", header=T, sep="\t")
samplesheet <- samplesheet[match(file_names, samplesheet$File.Name),]
colnames(counts) <- samplesheet$Sample.ID
meta <- subset(samplesheet, select=c(Sample.ID, Sample.Type,Subtype , Age, Race))
rownames(meta) <- NULL
meta <- column_to_rownames(meta, var="Sample.ID")

meta$Sample.Type <- gsub("Primary Tumor", "Tumor", meta$Sample.Type,)
meta$Sample.Type <- gsub("Solid Tissue Normal", "Normal", meta$Sample.Type)
meta$Sample.Type <- as.factor(meta$Sample.Type)


# Check if order of columns in count data frame is similar to that of sample annotation data frame
# This should return TRUE
table(colnames(counts) == rownames(meta))

# Check the tumor and normal frequency
table(meta[,"Sample.Type"])

# Check how many genes are not expressed in any of the samples?
rawCntDat<- counts
isZero<- rowSums(rawCntDat) == 0
zeroCnt<- sum(isZero, na.rm = TRUE)
print(zeroCnt)

# What percentage of the genes are zero expressed?
zeroPerc<- zeroCnt/nrow(rawCntDat)*100
print(zeroPerc)

# Filtering the genes that are not expressed in any sample
filterInd<- isZero
counts <- rawCntDat[!filterInd,]

# Data normalization
library(DESeq2)
vstcnt<- vst(as.matrix(counts))

# Plot VSTnormalization_sd_vs_mean for vstcnt and counts
pdf("VSTnormalization_sd_vs_mean.pdf")
par(mfrow=c(2,1), col="#66CCFF")
plot(rowMeans(counts), apply(counts, 1, sd), main="Not Normalized", pch = 16, ylab="SD", xlab="Mean")
plot(rowMeans(vstcnt), apply(vstcnt, 1, sd), ylim=c(0, 50), main="VST Normalized", pch = 16, ylab="SD", xlab="Mean")
dev.off()

#PCA analysis
# Run PCA on transpose of the vstcnt
genVstcntPca<- prcomp(t(vstcnt), scale=TRUE)

#Measure the variance and plot the variance of the data that
# is explained by each PC dimension
pcVar<- genVstcntPca$sdev^2
pdf("PC_Var.pdf")
barplot(pcVar, main="Variance Explained by each PC", col=rainbow(10), names=colnames(genVstcntPca$x), ylab="Var")
dev.off()

# Measure the variance percentage and
# plot the percentage of the variance within the data that
# is explained by each PC dimension
varPerc<- 100*pcVar/sum(pcVar)
pdf("PC_Var_Perc.pdf")
barplot(varPerc, main="Variance % Explained by each PC", col=rainbow(10), names=colnames(genVstcntPca$x), ylab="Var %")
dev.off()

#Now visualize to see how well the "Sample.Type" sample annotations in the "meta" (accessible
#using meta$Sample.Type) separate the expression data. Does it separate the data well?
pdf("pairs_sample_type.pdf")
pairs(genVstcntPca$x[,1:6], main = "Sample Type", col=c("black", "red")[as.numeric(as.factor(meta$Sample.Type))],  pch = 20, labels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"))
par(xpd = TRUE)
legend(0.05, 1, title="Sample Type", cex=0.4, fill = c("red", "black"), legend = (as.vector(unique(meta$Sample.Type))))
dev.off()
#only PC1 & PC2
pdf("pairs_sample_type2.pdf")
pairs(genVstcntPca$x[,1:2], main = "Sample Type", col=c("black", "red")[as.numeric(as.factor(meta$Sample.Type))],  pch = 20, labels = c("PC1", "PC2"))
par(xpd = TRUE)
legend(0.05, 1, title="Sample Type", cex=0.4, fill = c("red", "black"), legend = (as.vector(unique(meta$Sample.Type))))
dev.off()

#Visualize to see how well the "Subtype" in the "meta" (accessible using meta$Subtype)
#separate the expression data.
pdf("pairs_subtype.pdf")
pairs(genVstcntPca$x[,1:6], main = "Subtype", col=c("red", "black", "lightblue", "gold", "green", "purple")[as.numeric(as.factor(meta$Subtype))],  pch = 20, labels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"))
par(xpd = TRUE)
legend(0.05, 1, title="Subtype", cex=0.19, fill = c("red", "purple", "green", "lightblue", "black", "gold"), legend = (as.vector(unique(meta$Subtype))))
dev.off()
#only PC1 & PC2
pdf("pairs_subtype2.pdf")
pairs(genVstcntPca$x[,1:2], main = "Subtype", col=c("red", "black", "lightblue", "gold", "green", "purple")[as.numeric(as.factor(meta$Subtype))],  pch = 20, labels = c("PC1", "PC2"))
par(xpd = TRUE)
legend(0.05, 1, title="Subtype", cex=0.4, fill = c("red", "purple", "green", "lightblue", "black", "gold"), legend = (as.vector(unique(meta$Subtype))))
dev.off()

#Visualize to see how well the "Race" in the "meta" (accessible using meta$Race)
#separate the expression data.
pdf("pairs_race.pdf")
pairs(genVstcntPca$x[,1:6], main = "Race", col=c("red", "black", "lightblue")[as.numeric(as.factor(meta$Race))],  pch = 20, labels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"))
par(xpd = TRUE)
legend(0.05, 1, title="Race", cex=0.31, fill = c("lightblue", "black", "red"), legend = (as.vector(unique(meta$Race))))
dev.off()


#Visualize to see how well the "Age" in the "meta" (accessible using meta$Age)
#separate the expression data.
library(Hmisc)
age_cat <- cut2(meta$Age, c(45, 60, 75))
summary(age_cat)
meta <- cbind(meta, age_cat)

pdf("pairs_age.pdf")
pairs(genVstcntPca$x[,1:6], main = "Age", col=c("red", "black", "lightblue", "green")[as.numeric(as.factor(meta$age_cat))],  pch = 20, labels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6"))
par(xpd = TRUE)
legend(0.05, 1, title="Age", cex=0.26, fill = c("red", "black", "lightblue", "green"), legend = (as.vector(unique(meta$age_cat))))
dev.off()

# Differential gene expression analysis for all 116 tumors samples vs 113 normal ones using DESeq2
dds<- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData =
                                       meta[,c("Race", "age_cat", "Sample.Type")], design = ~Race+age_cat+Sample.Type)
##filtering low counts
##113 is the minimum number of samples (n LAR = 20)
dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds) >= 10) >= 20
dds <- dds[keep,]

##DeSeq
library(BiocParallel)
dds$Sample.Type <- relevel(dds$Sample.Type, ref = "Normal")
dds<- DESeq2::DESeq(dds, BPPARAM = SnowParam(workers=2))

##plotPCA
vsdata <- vst(dds)

###Sample.Type PCA plot
pdf("plotPCA_Sample_Type.pdf")
plotPCA(vsdata, intgroup = "Sample.Type")
dev.off()

###Race PCA plot
library(ggplot2)
pdf("plotPCA_Race.pdf")
pcaData <- plotPCA(vsdata, intgroup=c("Sample.Type", "Race"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Race, shape=Sample.Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

###Age PCA plot
pdf("plotPCA_age_cat.pdf")
pcaData <- plotPCA(vsdata, intgroup=c("Sample.Type", "age_cat"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=age_cat, shape=Sample.Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

##results
ddsDiff<- DESeq2::results(dds, name="Sample.Type_Tumor_vs_Normal", pAdjustMethod =
                            "BH", parallel = TRUE, BPPARAM = SnowParam(workers=2))

## Adjust and correct Log Fold change relative to mean expression (ASHR)
ddsAshr <- DESeq2::lfcShrink(dds, coef="Sample.Type_Tumor_vs_Normal", type="ashr",
                             parallel = TRUE, BPPARAM = SnowParam(workers=2))

##Plot LogFC_vs_mean_expr_Unadjusted
pdf("LogFC_vs_mean_expr_Unadjusted.pdf")
plotMA(main="LogFC vs mean Unadjusted", object = ddsDiff, ylim = c(-10, 10))
dev.off()

##Plot LogFC_vs_mean_expr_ASHR_adjusted
pdf("LogFC_vs_mean_expr_ASHR_adjusted.pdf")
plotMA(main="LogFC vs mean ASHR adjusted", object = ddsAshr, ylim = c(-10, 10))
dev.off()

##significantly differentially expressed genes
sigs <- na.omit(ddsAshr)
sigs <- sigs[which(sigs$padj < 0.05 & abs(sigs$log2FoldChange) > 1),]

##annotating the genes
library(biomaRt)
library(stringr)
gene_ids <- str_replace(rownames(sigs),
                        pattern = ".[0-9]+$",
                        replacement = "")
ensembl<- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
annoGen<- biomaRt::getBM( attributes=c("ensembl_gene_id",
                                       "external_gene_name", "gene_biotype",
                                       "chromosome_name", "start_position", "end_position"),
                          filters = "ensembl_gene_id",
                          values = unique(as.character(gene_ids)),
                          mart=ensembl, uniqueRows = TRUE)

##deleting rows with missing ensemble_ids
missingrows <- as.vector(gene_ids[ !gene_ids %in% annoGen$ensembl_gene_id ])
rownames(sigs) <- gene_ids
sigs <- sigs[!(row.names(sigs) %in% missingrows),]

##binding genes anotation data & significantly differentially expressed genes
outDf<- cbind(annoGen, as.data.frame(sigs))
write.csv(outDf, file = "sigs.csv", row.names = FALSE)



# Differential gene expression analysis (subtypes vs normal) DESeq2
dds2<- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData =
                                       meta[,c("Race", "age_cat", "Subtype")], design = ~Race+age_cat+Subtype)
##filtering low counts
##20 is the minimum number of samples (n LAR = 20)
dds2 <- estimateSizeFactors(dds2)
keep <- rowSums(counts(dds2) >= 10) >= 20
dds2 <- dds2[keep,]

##DeSeq
dds2$Subtype <- relevel(dds2$Subtype, ref = "N")
library(BiocParallel)
dds2<- DESeq2::DESeq(dds2, BPPARAM = SnowParam(workers=2))

##plotPCA
vsdata2 <- vst(dds2)

###Subtype PCA plot
library(ggplot2)
pdf("plotPCA_Subtype.pdf")
pcaData <- plotPCA(vsdata2, intgroup="Subtype", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Subtype, shape=Subtype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()

###Age & Subtype PCA plot
library(ggplot2)
pdf("plotPCA_Age&Subtype.pdf")
pcaData <- plotPCA(vsdata2, intgroup=c("Subtype", "age_cat"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=age_cat, shape=Subtype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()


##results
###checking resultsnames
resultsNames(dds2)

###BL1 Subtype
BL1_ddsDiff<- DESeq2::results(dds2, name="Subtype_BL1_vs_N", pAdjustMethod =
                            "BH", parallel = TRUE, BPPARAM = SnowParam(workers=2))

#### Adjust and correct Log Fold change relative to mean expression (ASHR)
BL1_ddsAshr <- DESeq2::lfcShrink(dds2, coef="Subtype_BL1_vs_N", type="ashr",
                             parallel = TRUE, BPPARAM = SnowParam(workers=2))

####Plot LogFC_vs_mean_expr_Unadjusted
pdf("BL1_LogFC_vs_mean_expr_Unadjusted.pdf")
plotMA(main="LogFC vs mean Unadjusted", object = BL1_ddsDiff, ylim = c(-10, 10))
dev.off()

####Plot LogFC_vs_mean_expr_ASHR_adjusted
pdf("BL1_LogFC_vs_mean_expr_ASHR_adjusted.pdf")
plotMA(main="LogFC vs mean ASHR adjusted", object = BL1_ddsAshr, ylim = c(-10, 10))
dev.off()

####significantly differentially expressed genes
BL1_sigs <- na.omit(BL1_ddsAshr)
BL1_sigs <- BL1_sigs[which(BL1_sigs$padj < 0.05 & abs(BL1_sigs$log2FoldChange) > 1),]

####annotating the genes
library(biomaRt)
library(stringr)
gene_ids <- str_replace(rownames(BL1_sigs),
                        pattern = ".[0-9]+$",
                        replacement = "")
ensembl<- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
annoGen<- biomaRt::getBM( attributes=c("ensembl_gene_id",
                                       "external_gene_name", "gene_biotype",
                                       "chromosome_name", "start_position", "end_position"),
                          filters = "ensembl_gene_id",
                          values = unique(as.character(gene_ids)),
                          mart=ensembl, uniqueRows = TRUE)

####deleting rows with missing ensemble_ids
missingrows <- as.vector(gene_ids[ !gene_ids %in% annoGen$ensembl_gene_id ])
rownames(BL1_sigs) <- gene_ids
BL1_sigs <- BL1_sigs[!(row.names(BL1_sigs) %in% missingrows),]

####binding genes anotation data & significantly differentially expressed genes
outDf<- cbind(annoGen, as.data.frame(BL1_sigs))
write.csv(outDf, file = "BL1_sigs.csv", row.names = FALSE)


###BL2 Subtype
BL2_ddsDiff<- DESeq2::results(dds2, name="Subtype_BL2_vs_N", pAdjustMethod =
                                "BH", parallel = TRUE, BPPARAM = SnowParam(workers=2))

#### Adjust and correct Log Fold change relative to mean expression (ASHR)
BL2_ddsAshr <- DESeq2::lfcShrink(dds2, coef="Subtype_BL2_vs_N", type="ashr",
                                 parallel = TRUE, BPPARAM = SnowParam(workers=2))

####Plot LogFC_vs_mean_expr_Unadjusted
pdf("BL2_LogFC_vs_mean_expr_Unadjusted.pdf")
plotMA(main="LogFC vs mean Unadjusted", object = BL2_ddsDiff, ylim = c(-10, 10))
dev.off()

####Plot LogFC_vs_mean_expr_ASHR_adjusted
pdf("BL2_LogFC_vs_mean_expr_ASHR_adjusted.pdf")
plotMA(main="LogFC vs mean ASHR adjusted", object = BL2_ddsAshr, ylim = c(-10, 10))
dev.off()

####significantly differentially expressed genes
BL2_sigs <- na.omit(BL2_ddsAshr)
BL2_sigs <- BL2_sigs[which(BL2_sigs$padj < 0.05 & abs(BL2_sigs$log2FoldChange) > 1),]

####annotating the genes
library(biomaRt)
library(stringr)
gene_ids <- str_replace(rownames(BL2_sigs),
                        pattern = ".[0-9]+$",
                        replacement = "")
ensembl<- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
annoGen<- biomaRt::getBM( attributes=c("ensembl_gene_id",
                                       "external_gene_name", "gene_biotype",
                                       "chromosome_name", "start_position", "end_position"),
                          filters = "ensembl_gene_id",
                          values = unique(as.character(gene_ids)),
                          mart=ensembl, uniqueRows = TRUE)

####deleting rows with missing ensemble_ids
missingrows <- as.vector(gene_ids[ !gene_ids %in% annoGen$ensembl_gene_id ])
rownames(BL2_sigs) <- gene_ids
BL2_sigs <- BL2_sigs[!(row.names(BL2_sigs) %in% missingrows),]

####binding genes anotation data & significantly differentially expressed genes
outDf<- cbind(annoGen, as.data.frame(BL2_sigs))
write.csv(outDf, file = "BL2_sigs.csv", row.names = FALSE)


###LAR Subtype
LAR_ddsDiff<- DESeq2::results(dds2, name="Subtype_LAR_vs_N", pAdjustMethod =
                                "BH", parallel = TRUE, BPPARAM = SnowParam(workers=2))

#### Adjust and correct Log Fold change relative to mean expression (ASHR)
LAR_ddsAshr <- DESeq2::lfcShrink(dds2, coef="Subtype_LAR_vs_N", type="ashr",
                                 parallel = TRUE, BPPARAM = SnowParam(workers=2))

####Plot LogFC_vs_mean_expr_Unadjusted
pdf("LAR_LogFC_vs_mean_expr_Unadjusted.pdf")
plotMA(main="LogFC vs mean Unadjusted", object = LAR_ddsDiff, ylim = c(-10, 10))
dev.off()

####Plot LogFC_vs_mean_expr_ASHR_adjusted
pdf("LAR_LogFC_vs_mean_expr_ASHR_adjusted.pdf")
plotMA(main="LogFC vs mean ASHR adjusted", object = LAR_ddsAshr, ylim = c(-10, 10))
dev.off()

####significantly differentially expressed genes
LAR_sigs <- na.omit(LAR_ddsAshr)
LAR_sigs <- LAR_sigs[which(LAR_sigs$padj < 0.05 & abs(LAR_sigs$log2FoldChange) > 1),]

####annotating the genes
library(biomaRt)
library(stringr)
gene_ids <- str_replace(rownames(LAR_sigs),
                        pattern = ".[0-9]+$",
                        replacement = "")
ensembl<- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
annoGen<- biomaRt::getBM( attributes=c("ensembl_gene_id",
                                       "external_gene_name", "gene_biotype",
                                       "chromosome_name", "start_position", "end_position"),
                          filters = "ensembl_gene_id",
                          values = unique(as.character(gene_ids)),
                          mart=ensembl, uniqueRows = TRUE)

####deleting rows with missing ensemble_ids
missingrows <- as.vector(gene_ids[ !gene_ids %in% annoGen$ensembl_gene_id ])
rownames(LAR_sigs) <- gene_ids
LAR_sigs <- LAR_sigs[!(row.names(LAR_sigs) %in% missingrows),]

####binding genes anotation data & significantly differentially expressed genes
outDf<- cbind(annoGen, as.data.frame(LAR_sigs))
write.csv(outDf, file = "LAR_sigs.csv", row.names = FALSE)


###M Subtype
M_ddsDiff<- DESeq2::results(dds2, name="Subtype_M_vs_N", pAdjustMethod =
                                "BH", parallel = TRUE, BPPARAM = SnowParam(workers=2))

#### Adjust and correct Log Fold change relative to mean expression (ASHR)
M_ddsAshr <- DESeq2::lfcShrink(dds2, coef="Subtype_M_vs_N", type="ashr",
                                 parallel = TRUE, BPPARAM = SnowParam(workers=2))

####Plot LogFC_vs_mean_expr_Unadjusted
pdf("M_LogFC_vs_mean_expr_Unadjusted.pdf")
plotMA(main="LogFC vs mean Unadjusted", object = M_ddsDiff, ylim = c(-10, 10))
dev.off()

####Plot LogFC_vs_mean_expr_ASHR_adjusted
pdf("M_LogFC_vs_mean_expr_ASHR_adjusted.pdf")
plotMA(main="LogFC vs mean ASHR adjusted", object = M_ddsAshr, ylim = c(-10, 10))
dev.off()

####significantly differentially expressed genes
M_sigs <- na.omit(M_ddsAshr)
M_sigs <- M_sigs[which(M_sigs$padj < 0.05 & abs(M_sigs$log2FoldChange) > 1),]

####annotating the genes
library(biomaRt)
library(stringr)
gene_ids <- str_replace(rownames(M_sigs),
                        pattern = ".[0-9]+$",
                        replacement = "")
ensembl<- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
annoGen<- biomaRt::getBM( attributes=c("ensembl_gene_id",
                                       "external_gene_name", "gene_biotype",
                                       "chromosome_name", "start_position", "end_position"),
                          filters = "ensembl_gene_id",
                          values = unique(as.character(gene_ids)),
                          mart=ensembl, uniqueRows = TRUE)

####deleting rows with missing ensemble_ids
missingrows <- as.vector(gene_ids[ !gene_ids %in% annoGen$ensembl_gene_id ])
rownames(M_sigs) <- gene_ids
M_sigs <- M_sigs[!(row.names(M_sigs) %in% missingrows),]

####binding genes anotation data & significantly differentially expressed genes
outDf<- cbind(annoGen, as.data.frame(M_sigs))
write.csv(outDf, file = "M_sigs.csv", row.names = FALSE)



# Gene Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DOSE)
library(ReactomePA)

## TNBC
sigs_up <- rownames(sigs[sigs$log2FoldChange > 1,])
sigs_down <- rownames(sigs[sigs$log2FoldChange < 1,])


### Gene Ontology

#### Biological Processes
tnbc_bp_up <- enrichGO(gene = sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_bp_up.pdf", height = 12)
plot(barplot(tnbc_bp_up, showCategory = 20))
dev.off()

tnbc_bp_down <- enrichGO(gene = sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_bp_down.pdf", height = 12)
plot(barplot(tnbc_bp_down, showCategory = 20))
dev.off()

#### Molecular Functions
tnbc_mf_up <- enrichGO(gene = sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                       ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_mf_up.pdf", height = 12)
plot(barplot(tnbc_mf_up, showCategory = 20))
dev.off()

tnbc_mf_down <- enrichGO(gene = sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                         ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_mf_down.pdf", height = 12)
plot(barplot(tnbc_mf_down, showCategory = 20))
dev.off()

#### Cellular Components

tnbc_cc_up <- enrichGO(gene = sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                       ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_cc_up.pdf", height = 12)
plot(barplot(tnbc_cc_up, showCategory = 20))
dev.off()

tnbc_cc_down <- enrichGO(gene = sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                         ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_cc_down.pdf", height = 12)
plot(barplot(tnbc_cc_down, showCategory = 20))
dev.off()

### Pathway Enrichment

tnbc_entrez_up <- mapIds(x = org.Hs.eg.db,
                         keys = sigs_up,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

tnbc_entrez_down <- mapIds(x = org.Hs.eg.db,
                         keys = sigs_down,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

#### KEGG

tnbc_kegg_up <- enrichKEGG(gene = tnbc_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_kegg_up.pdf", height = 12)
plot(barplot(tnbc_kegg_up, showCategory = 20))
dev.off()

tnbc_kegg_down <- enrichKEGG(gene = tnbc_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_kegg_down.pdf", height = 12)
plot(barplot(tnbc_kegg_down, showCategory = 20))
dev.off()

#### Reactome

library(ReactomePA)
tnbc_reactome_up <- enrichPathway(gene = tnbc_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_reactome_up.pdf", height = 12)
plot(barplot(tnbc_reactome_up, showCategory = 20))
dev.off()

tnbc_reactome_down <- enrichPathway(gene = tnbc_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_reactome_down.pdf", height = 12)
plot(barplot(tnbc_reactome_down, showCategory = 20))
dev.off()

#### Disease Ontology

library(DOSE)
tnbc_do_up <- enrichDO(gene = tnbc_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_do_up.pdf", height = 12)
plot(barplot(tnbc_do_up, showCategory = 20))
dev.off()


## BL1
bl1_sigs_up <- rownames(BL1_sigs[BL1_sigs$log2FoldChange > 1,])
bl1_sigs_down <- rownames(BL1_sigs[BL1_sigs$log2FoldChange < 1,])


### Gene Ontology

#### Biological Processes
bl1_bp_up <- enrichGO(gene = bl1_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_bp_up.pdf", height = 12)
plot(barplot(bl1_bp_up, showCategory = 20))
dev.off()

bl1_bp_down <- enrichGO(gene = bl1_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                         ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_bp_down.pdf", height = 12)
plot(barplot(bl1_bp_down, showCategory = 20))
dev.off()

#### Molecular Functions
bl1_mf_up <- enrichGO(gene = bl1_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                       ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_mf_up.pdf", height = 12)
plot(barplot(bl1_mf_up, showCategory = 20))
dev.off()

bl1_mf_down <- enrichGO(gene = bl1_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                         ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_mf_down.pdf", height = 12)
plot(barplot(bl1_mf_down, showCategory = 20))
dev.off()

#### Cellular Components

bl1_cc_up <- enrichGO(gene = bl1_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                       ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_cc_up.pdf", height = 12)
plot(barplot(bl1_cc_up, showCategory = 20))
dev.off()

bl1_cc_down <- enrichGO(gene = bl1_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                         ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_cc_down.pdf", height = 12)
plot(barplot(bl1_cc_down, showCategory = 20))
dev.off()

### Pathway Enrichment

bl1_entrez_up <- mapIds(x = org.Hs.eg.db,
                         keys = bl1_sigs_up,
                         column = "ENTREZID",
                         keytype = "ENSEMBL",
                         multiVals = "first")

bl1_entrez_down <- mapIds(x = org.Hs.eg.db,
                           keys = bl1_sigs_down,
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")

#### KEGG

bl1_kegg_up <- enrichKEGG(gene = bl1_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_kegg_up.pdf", height = 12)
plot(barplot(bl1_kegg_up, showCategory = 20))
dev.off()

bl1_kegg_down <- enrichKEGG(gene = bl1_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_kegg_down.pdf", height = 12)
plot(barplot(bl1_kegg_down, showCategory = 20))
dev.off()

#### Reactome

bl1_reactome_up <- enrichPathway(gene = bl1_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_reactome_up.pdf", height = 12)
plot(barplot(bl1_reactome_up, showCategory = 20))
dev.off()

bl1_reactome_down <- enrichPathway(gene = bl1_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_reactome_down.pdf", height = 13)
plot(barplot(bl1_reactome_down, showCategory = 20))
dev.off()

#### Disease Ontology

bl1_do_up <- enrichDO(gene = bl1_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_do_up.pdf", height = 12)
plot(barplot(bl1_do_up, showCategory = 20))
dev.off()


## BL2
bl2_sigs_up <- rownames(BL2_sigs[BL2_sigs$log2FoldChange > 1,])
bl2_sigs_down <- rownames(BL2_sigs[BL2_sigs$log2FoldChange < 1,])


### Gene Ontology

#### Biological Processes
bl2_bp_up <- enrichGO(gene = bl2_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_bp_up.pdf", height = 12)
plot(barplot(bl2_bp_up, showCategory = 20))
dev.off()

bl2_bp_down <- enrichGO(gene = bl2_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_bp_down.pdf", height = 12)
plot(barplot(bl2_bp_down, showCategory = 20))
dev.off()

#### Molecular Functions
bl2_mf_up <- enrichGO(gene = bl2_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_mf_up.pdf", height = 12)
plot(barplot(bl2_mf_up, showCategory = 20))
dev.off()

bl2_mf_down <- enrichGO(gene = bl2_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_mf_down.pdf", height = 12)
plot(barplot(bl2_mf_down, showCategory = 20))
dev.off()

#### Cellular Components

bl2_cc_up <- enrichGO(gene = bl2_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_cc_up.pdf", height = 12)
plot(barplot(bl2_cc_up, showCategory = 20))
dev.off()

bl2_cc_down <- enrichGO(gene = bl2_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_cc_down.pdf", height = 12)
plot(barplot(bl2_cc_down, showCategory = 20))
dev.off()

### Pathway Enrichment

bl2_entrez_up <- mapIds(x = org.Hs.eg.db,
                        keys = bl2_sigs_up,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")

bl2_entrez_down <- mapIds(x = org.Hs.eg.db,
                          keys = bl2_sigs_down,
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")

#### KEGG

bl2_kegg_up <- enrichKEGG(gene = bl2_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_kegg_up.pdf", height = 12)
plot(barplot(bl2_kegg_up, showCategory = 20))
dev.off()

bl2_kegg_down <- enrichKEGG(gene = bl2_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_kegg_down.pdf", height = 12)
plot(barplot(bl2_kegg_down, showCategory = 20))
dev.off()

#### Reactome

bl2_reactome_up <- enrichPathway(gene = bl2_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_reactome_up.pdf", height = 12)
plot(barplot(bl2_reactome_up, showCategory = 20))
dev.off()

bl2_reactome_down <- enrichPathway(gene = bl2_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_reactome_down.pdf", height = 12)
plot(barplot(bl2_reactome_down, showCategory = 20))
dev.off()

#### Disease Ontology

bl2_do_up <- enrichDO(gene = bl2_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_do_up.pdf", height = 12)
plot(barplot(bl2_do_up, showCategory = 20))
dev.off()


## M
m_sigs_up <- rownames(M_sigs[M_sigs$log2FoldChange > 1,])
m_sigs_down <- rownames(M_sigs[M_sigs$log2FoldChange < 1,])


### Gene Ontology

#### Biological Processes
m_bp_up <- enrichGO(gene = m_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_bp_up.pdf", height = 12)
plot(barplot(m_bp_up, showCategory = 20))
dev.off()

m_bp_down <- enrichGO(gene = m_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_bp_down.pdf", height = 12)
plot(barplot(m_bp_down, showCategory = 20))
dev.off()

#### Molecular Functions
m_mf_up <- enrichGO(gene = m_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_mf_up.pdf", height = 12)
plot(barplot(m_mf_up, showCategory = 20))
dev.off()

m_mf_down <- enrichGO(gene = m_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_mf_down.pdf", height = 12)
plot(barplot(m_mf_down, showCategory = 20))
dev.off()

#### Cellular Components

m_cc_up <- enrichGO(gene = m_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_cc_up.pdf", height = 12)
plot(barplot(m_cc_up, showCategory = 20))
dev.off()

m_cc_down <- enrichGO(gene = m_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_cc_down.pdf", height = 12)
plot(barplot(m_cc_down, showCategory = 20))
dev.off()

### Pathway Enrichment

m_entrez_up <- mapIds(x = org.Hs.eg.db,
                        keys = m_sigs_up,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")

m_entrez_down <- mapIds(x = org.Hs.eg.db,
                          keys = m_sigs_down,
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")

#### KEGG

m_kegg_up <- enrichKEGG(gene = m_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_kegg_up.pdf", height = 12)
plot(barplot(m_kegg_up, showCategory = 20))
dev.off()

m_kegg_down <- enrichKEGG(gene = m_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_kegg_down.pdf", height = 12)
plot(barplot(m_kegg_down, showCategory = 20))
dev.off()

#### Reactome

m_reactome_up <- enrichPathway(gene = m_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_reactome_up.pdf", height = 12)
plot(barplot(m_reactome_up, showCategory = 20))
dev.off()

m_reactome_down <- enrichPathway(gene = m_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_reactome_down.pdf", height = 13)
plot(barplot(m_reactome_down, showCategory = 20))
dev.off()

#### Disease Ontology

m_do_up <- enrichDO(gene = m_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_do_up.pdf", height = 12)
plot(barplot(m_do_up, showCategory = 20))
dev.off()


## LAR
lar_sigs_up <- rownames(LAR_sigs[LAR_sigs$log2FoldChange > 1,])
lar_sigs_down <- rownames(LAR_sigs[LAR_sigs$log2FoldChange < 1,])


### Gene Ontology

#### Biological Processes
lar_bp_up <- enrichGO(gene = lar_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_bp_up.pdf", height = 12)
plot(barplot(lar_bp_up, showCategory = 20))
dev.off()

lar_bp_down <- enrichGO(gene = lar_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_bp_down.pdf", height = 12)
plot(barplot(lar_bp_down, showCategory = 20))
dev.off()

#### Molecular Functions
lar_mf_up <- enrichGO(gene = lar_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_mf_up.pdf", height = 12)
plot(barplot(lar_mf_up, showCategory = 20))
dev.off()

lar_mf_down <- enrichGO(gene = lar_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_mf_down.pdf", height = 12)
plot(barplot(lar_mf_down, showCategory = 20))
dev.off()

#### Cellular Components

lar_cc_up <- enrichGO(gene = lar_sigs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_cc_up.pdf", height = 12)
plot(barplot(lar_cc_up, showCategory = 20))
dev.off()

lar_cc_down <- enrichGO(gene = lar_sigs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_cc_down.pdf", height = 12)
plot(barplot(lar_cc_down, showCategory = 20))
dev.off()

### Pathway Enrichment

lar_entrez_up <- mapIds(x = org.Hs.eg.db,
                        keys = lar_sigs_up,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")

lar_entrez_down <- mapIds(x = org.Hs.eg.db,
                          keys = lar_sigs_down,
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")

#### KEGG

lar_kegg_up <- enrichKEGG(gene = lar_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_kegg_up.pdf", height = 12)
plot(barplot(lar_kegg_up, showCategory = 20))
dev.off()

lar_kegg_down <- enrichKEGG(gene = lar_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_kegg_down.pdf", height = 12)
plot(barplot(lar_kegg_down, showCategory = 20))
dev.off()

#### Reactome

lar_reactome_up <- enrichPathway(gene = lar_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_reactome_up.pdf", height = 12)
plot(barplot(lar_reactome_up, showCategory = 20))
dev.off()

lar_reactome_down <- enrichPathway(gene = lar_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_reactome_down.pdf", height = 12)
plot(barplot(lar_reactome_down, showCategory = 20))
dev.off()

#### Disease Ontology

lar_do_up <- enrichDO(gene = lar_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_do_up.pdf", height = 12)
plot(barplot(lar_do_up, showCategory = 20))
dev.off()


## TNBC Hub Genes
tnbc_hubs_up <- read.csv("tnbc_up.csv")
tnbc_hubs_up <- tnbc_hubs_up$id

tnbc_hubs_down <- read.csv("tnbc_down.csv")
tnbc_hubs_down <- tnbc_hubs_down$id


### Gene Ontology

#### Biological Processes
tnbc_hub_bp_up <- enrichGO(gene = tnbc_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_hub_bp_up.pdf", height = 8)
plot(dotplot(tnbc_hub_bp_up, showCategory = 10))
dev.off()

tnbc_hub_bp_down <- enrichGO(gene = tnbc_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_hub_bp_down.pdf", height = 8)
plot(dotplot(tnbc_hub_bp_down, showCategory = 10))
dev.off()

#### Molecular Functions
tnbc_hub_mf_up <- enrichGO(gene = tnbc_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_hub_mf_up.pdf", height = 8)
plot(dotplot(tnbc_hub_mf_up, showCategory = 10))
dev.off()

tnbc_hub_mf_down <- enrichGO(gene = tnbc_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_hub_mf_down.pdf", height = 12)
plot(dotplot(tnbc_hub_mf_down, showCategory = 10))
dev.off()

#### Cellular Components

tnbc_hub_cc_up <- enrichGO(gene = tnbc_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                      ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_hub_cc_up.pdf", height = 8)
plot(dotplot(tnbc_hub_cc_up, showCategory = 10))
dev.off()

tnbc_hub_cc_down <- enrichGO(gene = tnbc_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                        ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("tnbc_hub_cc_down.pdf", height = 8)
plot(dotplot(tnbc_hub_cc_down, showCategory = 10))
dev.off()

### Pathway Enrichment

tnbc_hub_entrez_up <- mapIds(x = org.Hs.eg.db,
                        keys = tnbc_hubs_up,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")

tnbc_hub_entrez_down <- mapIds(x = org.Hs.eg.db,
                          keys = tnbc_hubs_down,
                          column = "ENTREZID",
                          keytype = "ENSEMBL",
                          multiVals = "first")

#### KEGG

tnbc_hub_kegg_up <- enrichKEGG(gene = tnbc_hub_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_hub_kegg_up.pdf", height = 8)
plot(dotplot(tnbc_hub_kegg_up, showCategory = 10))
dev.off()

tnbc_hub_kegg_down <- enrichKEGG(gene = tnbc_hub_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_hub_kegg_down.pdf", height = 8)
plot(dotplot(tnbc_hub_kegg_down, showCategory = 10))
dev.off()

#### Reactome

tnbc_hub_reactome_up <- enrichPathway(gene = tnbc_hub_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_hub_reactome_up.pdf", height = 8)
plot(dotplot(tnbc_hub_reactome_up, showCategory = 10))
dev.off()

tnbc_hub_reactome_down <- enrichPathway(gene = tnbc_hub_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_hub_reactome_down.pdf", height = 8)
plot(dotplot(tnbc_hub_reactome_down, showCategory = 10))
dev.off()

#### Disease Ontology

tnbc_hub_do_up <- enrichDO(gene = tnbc_hub_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("tnbc_hub_do_up.pdf", height = 8)
plot(dotplot(tnbc_hub_do_up, showCategory = 10))
dev.off()


## BL1 Hub Genes
bl1_hubs_up <- read.csv("bl1_up.csv")
bl1_hubs_up <- bl1_hubs_up$id

bl1_hubs_down <- read.csv("bl1_down.csv")
bl1_hubs_down <- bl1_hubs_down$id


### Gene Ontology

#### Biological Processes
bl1_hub_bp_up <- enrichGO(gene = bl1_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_hub_bp_up.pdf", height = 8)
plot(dotplot(bl1_hub_bp_up, showCategory = 10))
dev.off()

bl1_hub_bp_down <- enrichGO(gene = bl1_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_hub_bp_down.pdf", height = 8)
plot(dotplot(bl1_hub_bp_down, showCategory = 10))
dev.off()

#### Molecular Functions
bl1_hub_mf_up <- enrichGO(gene = bl1_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_hub_mf_up.pdf", height = 8)
plot(dotplot(bl1_hub_mf_up, showCategory = 10))
dev.off()

bl1_hub_mf_down <- enrichGO(gene = bl1_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_hub_mf_down.pdf", height = 12)
plot(dotplot(bl1_hub_mf_down, showCategory = 10))
dev.off()

#### Cellular Components

bl1_hub_cc_up <- enrichGO(gene = bl1_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_hub_cc_up.pdf", height = 8)
plot(dotplot(bl1_hub_cc_up, showCategory = 10))
dev.off()

bl1_hub_cc_down <- enrichGO(gene = bl1_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl1_hub_cc_down.pdf", height = 8)
plot(dotplot(bl1_hub_cc_down, showCategory = 10))
dev.off()

### Pathway Enrichment

bl1_hub_entrez_up <- mapIds(x = org.Hs.eg.db,
                             keys = bl1_hubs_up,
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals = "first")

bl1_hub_entrez_down <- mapIds(x = org.Hs.eg.db,
                               keys = bl1_hubs_down,
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")

#### KEGG

bl1_hub_kegg_up <- enrichKEGG(gene = bl1_hub_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_hub_kegg_up.pdf", height = 8)
plot(dotplot(bl1_hub_kegg_up, showCategory = 10))
dev.off()

bl1_hub_kegg_down <- enrichKEGG(gene = bl1_hub_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_hub_kegg_down.pdf", height = 8)
plot(dotplot(bl1_hub_kegg_down, showCategory = 10))
dev.off()

#### Reactome

bl1_hub_reactome_up <- enrichPathway(gene = bl1_hub_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_hub_reactome_up.pdf", height = 8)
plot(dotplot(bl1_hub_reactome_up, showCategory = 10))
dev.off()

bl1_hub_reactome_down <- enrichPathway(gene = bl1_hub_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_hub_reactome_down.pdf", height = 8)
plot(dotplot(bl1_hub_reactome_down, showCategory = 10))
dev.off()

#### Disease Ontology

bl1_hub_do_up <- enrichDO(gene = bl1_hub_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl1_hub_do_up.pdf", height = 8)
plot(dotplot(bl1_hub_do_up, showCategory = 10))
dev.off()


## BL2 Hub Genes
bl2_hubs_up <- read.csv("bl2_up.csv")
bl2_hubs_up <- bl2_hubs_up$id

bl2_hubs_down <- read.csv("bl2_down.csv")
bl2_hubs_down <- bl2_hubs_down$id


### Gene Ontology

#### Biological Processes
bl2_hub_bp_up <- enrichGO(gene = bl2_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_hub_bp_up.pdf", height = 8)
plot(dotplot(bl2_hub_bp_up, showCategory = 10))
dev.off()

bl2_hub_bp_down <- enrichGO(gene = bl2_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_hub_bp_down.pdf", height = 8)
plot(dotplot(bl2_hub_bp_down, showCategory = 10))
dev.off()

#### Molecular Functions
bl2_hub_mf_up <- enrichGO(gene = bl2_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_hub_mf_up.pdf", height = 8)
plot(dotplot(bl2_hub_mf_up, showCategory = 10))
dev.off()

bl2_hub_mf_down <- enrichGO(gene = bl2_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_hub_mf_down.pdf", height = 11)
plot(dotplot(bl2_hub_mf_down, showCategory = 10))
dev.off()

#### Cellular Components

bl2_hub_cc_up <- enrichGO(gene = bl2_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_hub_cc_up.pdf", height = 8)
plot(dotplot(bl2_hub_cc_up, showCategory = 10))
dev.off()

bl2_hub_cc_down <- enrichGO(gene = bl2_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("bl2_hub_cc_down.pdf", height = 8)
plot(dotplot(bl2_hub_cc_down, showCategory = 10))
dev.off()

### Pathway Enrichment

bl2_hub_entrez_up <- mapIds(x = org.Hs.eg.db,
                             keys = bl2_hubs_up,
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals = "first")

bl2_hub_entrez_down <- mapIds(x = org.Hs.eg.db,
                               keys = bl2_hubs_down,
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")

#### KEGG

bl2_hub_kegg_up <- enrichKEGG(gene = bl2_hub_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_hub_kegg_up.pdf", height = 8)
plot(dotplot(bl2_hub_kegg_up, showCategory = 10))
dev.off()

bl2_hub_kegg_down <- enrichKEGG(gene = bl2_hub_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_hub_kegg_down.pdf", height = 8)
plot(dotplot(bl2_hub_kegg_down, showCategory = 10))
dev.off()

#### Reactome

bl2_hub_reactome_up <- enrichPathway(gene = bl2_hub_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_hub_reactome_up.pdf", height = 8)
plot(dotplot(bl2_hub_reactome_up, showCategory = 10))
dev.off()

bl2_hub_reactome_down <- enrichPathway(gene = bl2_hub_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_hub_reactome_down.pdf", height = 8)
plot(dotplot(bl2_hub_reactome_down, showCategory = 10))
dev.off()

#### Disease Ontology

bl2_hub_do_up <- enrichDO(gene = bl2_hub_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("bl2_hub_do_up.pdf", height = 8)
plot(dotplot(bl2_hub_do_up, showCategory = 10))
dev.off()


## M Hub Genes
m_hubs_up <- read.csv("m_up.csv")
m_hubs_up <- m_hubs_up$id

m_hubs_down <- read.csv("m_down.csv")
m_hubs_down <- m_hubs_down$id


### Gene Ontology

#### Biological Processes
m_hub_bp_up <- enrichGO(gene = m_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_hub_bp_up.pdf", height = 8)
plot(dotplot(m_hub_bp_up, showCategory = 10))
dev.off()

m_hub_bp_down <- enrichGO(gene = m_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_hub_bp_down.pdf", height = 8)
plot(dotplot(m_hub_bp_down, showCategory = 10))
dev.off()

#### Molecular Functions
m_hub_mf_up <- enrichGO(gene = m_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_hub_mf_up.pdf", height = 8)
plot(dotplot(m_hub_mf_up, showCategory = 10))
dev.off()

m_hub_mf_down <- enrichGO(gene = m_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_hub_mf_down.pdf", height = 8)
plot(dotplot(m_hub_mf_down, showCategory = 10))
dev.off()

#### Cellular Components

m_hub_cc_up <- enrichGO(gene = m_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_hub_cc_up.pdf", height = 8)
plot(dotplot(m_hub_cc_up, showCategory = 10))
dev.off()

m_hub_cc_down <- enrichGO(gene = m_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("m_hub_cc_down.pdf", height = 8)
plot(dotplot(m_hub_cc_down, showCategory = 10))
dev.off()

### Pathway Enrichment

m_hub_entrez_up <- mapIds(x = org.Hs.eg.db,
                             keys = m_hubs_up,
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals = "first")

m_hub_entrez_down <- mapIds(x = org.Hs.eg.db,
                               keys = m_hubs_down,
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")

#### KEGG

m_hub_kegg_up <- enrichKEGG(gene = m_hub_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_hub_kegg_up.pdf", height = 8)
plot(dotplot(m_hub_kegg_up, showCategory = 10))
dev.off()

m_hub_kegg_down <- enrichKEGG(gene = m_hub_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_hub_kegg_down.pdf", height = 8)
plot(dotplot(m_hub_kegg_down, showCategory = 10))
dev.off()

#### Reactome

m_hub_reactome_up <- enrichPathway(gene = m_hub_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_hub_reactome_up.pdf", height = 8)
plot(dotplot(m_hub_reactome_up, showCategory = 10))
dev.off()

m_hub_reactome_down <- enrichPathway(gene = m_hub_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_hub_reactome_down.pdf", height = 8)
plot(dotplot(m_hub_reactome_down, showCategory = 10))
dev.off()

#### Disease Ontology

m_hub_do_up <- enrichDO(gene = m_hub_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("m_hub_do_up.pdf", height = 8)
plot(dotplot(m_hub_do_up, showCategory = 10))
dev.off()


## LAR Hub Genes
lar_hubs_up <- read.csv("lar_up.csv")
lar_hubs_up <- lar_hubs_up$id

lar_hubs_down <- read.csv("lar_down.csv")
lar_hubs_down <- lar_hubs_down$id


### Gene Ontology

#### Biological Processes
lar_hub_bp_up <- enrichGO(gene = lar_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_hub_bp_up.pdf", height = 8)
plot(dotplot(lar_hub_bp_up, showCategory = 10))
dev.off()

lar_hub_bp_down <- enrichGO(gene = lar_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_hub_bp_down.pdf", height = 8)
plot(dotplot(lar_hub_bp_down, showCategory = 10))
dev.off()

#### Molecular Functions
lar_hub_mf_up <- enrichGO(gene = lar_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_hub_mf_up.pdf", height = 8)
plot(dotplot(lar_hub_mf_up, showCategory = 10))
dev.off()

lar_hub_mf_down <- enrichGO(gene = lar_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_hub_mf_down.pdf", height = 8)
plot(dotplot(lar_hub_mf_down, showCategory = 10))
dev.off()

#### Cellular Components

lar_hub_cc_up <- enrichGO(gene = lar_hubs_up, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                           ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_hub_cc_up.pdf", height = 8)
plot(dotplot(lar_hub_cc_up, showCategory = 10))
dev.off()

lar_hub_cc_down <- enrichGO(gene = lar_hubs_down, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL",
                             ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("lar_hub_cc_down.pdf", height = 8)
plot(dotplot(lar_hub_cc_down, showCategory = 10))
dev.off()

### Pathway Enrichment

lar_hub_entrez_up <- mapIds(x = org.Hs.eg.db,
                             keys = lar_hubs_up,
                             column = "ENTREZID",
                             keytype = "ENSEMBL",
                             multiVals = "first")

lar_hub_entrez_down <- mapIds(x = org.Hs.eg.db,
                               keys = lar_hubs_down,
                               column = "ENTREZID",
                               keytype = "ENSEMBL",
                               multiVals = "first")

#### KEGG

lar_hub_kegg_up <- enrichKEGG(gene = lar_hub_entrez_up, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_hub_kegg_up.pdf", height = 8)
plot(dotplot(lar_hub_kegg_up, showCategory = 10))
dev.off()

lar_hub_kegg_down <- enrichKEGG(gene = lar_hub_entrez_down, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_hub_kegg_down.pdf", height = 8)
plot(dotplot(lar_hub_kegg_down, showCategory = 10))
dev.off()

#### Reactome

lar_hub_reactome_up <- enrichPathway(gene = lar_hub_entrez_up, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_hub_reactome_up.pdf", height = 8)
plot(dotplot(lar_hub_reactome_up, showCategory = 10))
dev.off()

lar_hub_reactome_down <- enrichPathway(gene = lar_hub_entrez_down, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_hub_reactome_down.pdf", height = 8)
plot(dotplot(lar_hub_reactome_down, showCategory = 10))
dev.off()

#### Disease Ontology

lar_hub_do_up <- enrichDO(gene = lar_hub_entrez_up, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("lar_hub_do_up.pdf", height = 8)
plot(dotplot(lar_hub_do_up, showCategory = 10))
dev.off()



# Hub Genes Heatmap
rownames(vsdata2) <- str_replace(rownames(vsdata2),
                                 pattern = ".[0-9]+$",
                                 replacement = "")


## TNBC
tnbc_hub_sigs_up <- sigs[as.vector(tnbc_hubs_up),]
tnbc_hubs_up_symbols <- read.csv("tnbc_up.csv") 
tnbc_hubs_up_symbols <- tnbc_hubs_up_symbols$symbol  
tnbc_hub_sigs_up$symbol <- tnbc_hubs_up_symbols
###top 100 upregulated hub genes ###changed to 20 after re-run
tnbc_hub_sigs_heatmap <- tnbc_hub_sigs_up[1:20,]

tnbc_hub_sigs_down <- sigs[as.vector(tnbc_hubs_down),]
tnbc_hubs_down_symbols <- read.csv("tnbc_down.csv") 
tnbc_hubs_down_symbols <- tnbc_hubs_down_symbols$symbol  
tnbc_hub_sigs_down$symbol <- tnbc_hubs_down_symbols
###top 10 downregulated hub genes ###changed to 6 after re_run
tnbc_hub_sigs_heatmap <- rbind(tnbc_hub_sigs_heatmap,tnbc_hub_sigs_down[1:6,])

tnbc_hub_sigs_heatmap <- tnbc_hub_sigs_heatmap[order(tnbc_hub_sigs_heatmap$log2FoldChange, decreasing = TRUE),]

###using vst normalized count data from before
##### alternatively can use rlog with this command: rlog_out <- rlog(dds2, blind = FALSE) 
mat_tnbc <- assay(vsdata2)[rownames(tnbc_hub_sigs_heatmap), rownames(meta)] ###hub genes x samples
colnames(mat_tnbc) <- rownames(meta)
base_mean_tnbc <- rowMeans(mat_tnbc)
###center and scale each column (Z-score) then transpose
mat_tnbc_scaled <- t(apply(mat_tnbc, 1, scale))
colnames(mat_tnbc_scaled) <- colnames(mat_tnbc)

###getting log2foldchange for each gene
l2fc_val_tnbc <- as.matrix(tnbc_hub_sigs_heatmap$log2FoldChange)
colnames(l2fc_val_tnbc) <- "logFC"

###getting mean value for each gene
mean_tnbc <- as.matrix(tnbc_hub_sigs_heatmap$baseMean)
colnames(mean_tnbc) <- "AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2fc values
col_logFC_tnbc <- colorRamp2(c(min(l2fc_val_tnbc), 0, max(l2fc_val_tnbc)), c("blue", "white", "red"))

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr_tnbc <- colorRamp2(c(quantile(mean_tnbc)[1], quantile(mean_tnbc)[4]), c("yellow", "orange"))




ha_main_tnbc <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2) ,
                                                   labels = c("Normal", "TNBC"),
                                                   labels_gp = gpar(col = "black", fontsize = 20)))

ha_l2fc_tnbc <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill =2),
                                                         height = unit(1, "cm")))

h1_tnbc <- Heatmap(mat_tnbc_scaled, cluster_rows = F, top_annotation = ha_main_tnbc, show_column_names = F,
                   column_km = 2, name = "Z-score", cluster_columns = F)

h2_tnbc <- Heatmap(l2fc_val_tnbc, row_labels = tnbc_hub_sigs_heatmap$symbol,
                   cluster_rows = F, name = "logFC", top_annotation = ha_l2fc_tnbc, col = col_logFC_tnbc,
                   cell_fun = function(j, i, x, y, w, h, col) {  #add text to each grid
                     grid.text(round(l2fc_val_tnbc[i, j],2), x, y)
                   })

h3_tnbc <- Heatmap(mean_tnbc, row_labels = tnbc_hub_sigs_heatmap$symbol,
                   cluster_rows = F, name = "AveExpr", col = col_AveExpr_tnbc,
                   cell_fun = function(j, i, x, y, w, h, col){   #add text to each grid
                     grid.text(round(mean_tnbc[i, j],2), x, y)
                   })

pdf("tnbc_heatmap_100x10.pdf", width = 125, height = 55)
h_tnbc <- h1_tnbc + h2_tnbc + h3_tnbc
h_tnbc
dev.off()

## re-run but with 20x6 genes
pdf("tnbc_heatmap_20x6.pdf", width = 125, height = 15)
h_tnbc <- h1_tnbc + h2_tnbc + h3_tnbc
h_tnbc
dev.off()



## BL1
bl1_hub_sigs_up <- BL1_sigs[as.vector(bl1_hubs_up),]
bl1_hubs_up_symbols <- read.csv("bl1_up.csv") 
bl1_hubs_up_symbols <- bl1_hubs_up_symbols$symbol  
bl1_hub_sigs_up$symbol <- bl1_hubs_up_symbols
###top 100 upregulated hub genes ###changed to 20 after re-run
bl1_hub_sigs_heatmap <- bl1_hub_sigs_up[1:20,]

bl1_hub_sigs_down <- BL1_sigs[as.vector(bl1_hubs_down),]
bl1_hubs_down_symbols <- read.csv("bl1_down.csv") 
bl1_hubs_down_symbols <- bl1_hubs_down_symbols$symbol  
bl1_hub_sigs_down$symbol <- bl1_hubs_down_symbols
###top 10 downregulated hub genes ###changed to 6 after re_run
bl1_hub_sigs_heatmap <- rbind(bl1_hub_sigs_heatmap,bl1_hub_sigs_down[1:6,])

bl1_hub_sigs_heatmap <- bl1_hub_sigs_heatmap[order(bl1_hub_sigs_heatmap$log2FoldChange, decreasing = TRUE),]

###using vst normalized count data from before
#### alternatively can use rlog with this command: rlog_out <- rlog(dds2, blind = FALSE) 

###hub genes x samples
mat_bl1 <- assay(vsdata2)[rownames(bl1_hub_sigs_heatmap), rownames(meta[which(meta$Subtype == c("BL1", "N")),])]
colnames(mat_bl1) <- rownames(meta[which(meta$Subtype == c("BL1", "N")),])
base_mean_bl1 <- rowMeans(mat_bl1)
###center and scale each column (Z-score) then transpose
mat_bl1_scaled <- t(apply(mat_bl1, 1, scale))
colnames(mat_bl1_scaled) <- colnames(mat_bl1)

###getting log2foldchange for each gene
l2fc_val_bl1 <- as.matrix(bl1_hub_sigs_heatmap$log2FoldChange)
colnames(l2fc_val_bl1) <- "logFC"

###getting mean value for each gene
mean_bl1 <- as.matrix(bl1_hub_sigs_heatmap$baseMean)
colnames(mean_bl1) <- "AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2fc values
col_logFC_bl1 <- colorRamp2(c(min(l2fc_val_bl1), 0, max(l2fc_val_bl1)), c("blue", "white", "red"))

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr_bl1 <- colorRamp2(c(quantile(mean_bl1)[1], quantile(mean_bl1)[4]), c("yellow", "orange"))




ha_main_bl1 <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2) ,
                                 labels = c("Normal", "BL1"),
                                 labels_gp = gpar(col = "black", fontsize = 20)))

ha_l2fc_bl1 <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill =2),
                                                        height = unit(1, "cm")))

h1_bl1 <- Heatmap(mat_bl1_scaled, cluster_rows = F, top_annotation = ha_main_bl1, show_column_names = F,
                  column_km = 2, name = "Z-score", cluster_columns = F)

h2_bl1 <- Heatmap(l2fc_val_bl1, row_labels = bl1_hub_sigs_heatmap$symbol,
                  cluster_rows = F, name = "logFC", top_annotation = ha_l2fc_bl1, col = col_logFC_bl1,
                  cell_fun = function(j, i, x, y, w, h, col) {  #add text to each grid
                    grid.text(round(l2fc_val_bl1[i, j],2), x, y)
                  })

h3_bl1 <- Heatmap(mean_bl1, row_labels = bl1_hub_sigs_heatmap$symbol,
                  cluster_rows = F, name = "AveExpr", col = col_AveExpr_bl1,
                  cell_fun = function(j, i, x, y, w, h, col){   #add text to each grid
                    grid.text(round(mean_bl1[i, j],2), x, y)
                  })

pdf("bl1_heatmap_100x10.pdf", width = 40, height = 55)
h_bl1 <- h1_bl1 + h2_bl1 + h3_bl1
h_bl1
dev.off()

## re-run but with 20x6 genes
pdf("bl1_heatmap_20x6.pdf", width = 40, height = 15)
h_bl1 <- h1_bl1 + h2_bl1 + h3_bl1
h_bl1
dev.off()



## BL2
bl2_hub_sigs_up <- BL2_sigs[as.vector(bl2_hubs_up),]
bl2_hubs_up_symbols <- read.csv("bl2_up.csv") 
bl2_hubs_up_symbols <- bl2_hubs_up_symbols$symbol  
bl2_hub_sigs_up$symbol <- bl2_hubs_up_symbols
###top 100 upregulated hub genes ###changed to 20 after re-run
bl2_hub_sigs_heatmap <- bl2_hub_sigs_up[1:20,]

bl2_hub_sigs_down <- BL2_sigs[as.vector(bl2_hubs_down),]
bl2_hubs_down_symbols <- read.csv("bl2_down.csv") 
bl2_hubs_down_symbols <- bl2_hubs_down_symbols$symbol  
bl2_hub_sigs_down$symbol <- bl2_hubs_down_symbols
###top 10 downregulated hub genes ###changed to 6 after re_run
bl2_hub_sigs_heatmap <- rbind(bl2_hub_sigs_heatmap,bl2_hub_sigs_down[1:6,])

bl2_hub_sigs_heatmap <- bl2_hub_sigs_heatmap[order(bl2_hub_sigs_heatmap$log2FoldChange, decreasing = TRUE),]

###using vst normalized count data from before
#### alternatively can use rlog with this command: rlog_out <- rlog(dds2, blind = FALSE) 

###hub genes x samples
mat_bl2 <- assay(vsdata2)[rownames(bl2_hub_sigs_heatmap), rownames(meta[which(meta$Subtype == c("BL2", "N")),])]
colnames(mat_bl2) <- rownames(meta[which(meta$Subtype == c("BL2", "N")),])
base_mean_bl2 <- rowMeans(mat_bl2)
###center and scale each column (Z-score) then transpose
mat_bl2_scaled <- t(apply(mat_bl2, 1, scale))
colnames(mat_bl2_scaled) <- colnames(mat_bl2)

###getting log2foldchange for each gene
l2fc_val_bl2 <- as.matrix(bl2_hub_sigs_heatmap$log2FoldChange)
colnames(l2fc_val_bl2) <- "logFC"

###getting mean value for each gene
mean_bl2 <- as.matrix(bl2_hub_sigs_heatmap$baseMean)
colnames(mean_bl2) <- "AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2fc values
col_logFC_bl2 <- colorRamp2(c(min(l2fc_val_bl2), 0, max(l2fc_val_bl2)), c("blue", "white", "red"))

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr_bl2 <- colorRamp2(c(quantile(mean_bl2)[1], quantile(mean_bl2)[4]), c("yellow", "orange"))




ha_main_bl2 <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2) ,
                                                  labels = c("Normal", "BL2"),
                                                  labels_gp = gpar(col = "black", fontsize = 20)))

ha_l2fc_bl2 <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill =2),
                                                        height = unit(1, "cm")))

h1_bl2 <- Heatmap(mat_bl2_scaled, cluster_rows = F, top_annotation = ha_main_bl2, show_column_names = F,
                  column_km = 2, name = "Z-score", cluster_columns = F)

h2_bl2 <- Heatmap(l2fc_val_bl2, row_labels = bl2_hub_sigs_heatmap$symbol,
                  cluster_rows = F, name = "logFC", top_annotation = ha_l2fc_bl2, col = col_logFC_bl2,
                  cell_fun = function(j, i, x, y, w, h, col) {  #add text to each grid
                    grid.text(round(l2fc_val_bl2[i, j],2), x, y)
                  })

h3_bl2 <- Heatmap(mean_bl2, row_labels = bl2_hub_sigs_heatmap$symbol,
                  cluster_rows = F, name = "AveExpr", col = col_AveExpr_bl2,
                  cell_fun = function(j, i, x, y, w, h, col){   #add text to each grid
                    grid.text(round(mean_bl2[i, j],2), x, y)
                  })

pdf("bl2_heatmap_100x10.pdf", width = 40, height = 55)
h_bl2 <- h1_bl2 + h2_bl2 + h3_bl2
h_bl2
dev.off()

## re-run but with 20x6 genes
pdf("bl2_heatmap_20x6.pdf", width = 40, height = 15)
h_bl2 <- h1_bl2 + h2_bl2 + h3_bl2
h_bl2
dev.off()



## M
m_hub_sigs_up <- M_sigs[as.vector(m_hubs_up),]
m_hubs_up_symbols <- read.csv("m_up.csv") 
m_hubs_up_symbols <- m_hubs_up_symbols$symbol  
m_hub_sigs_up$symbol <- m_hubs_up_symbols
###top 100 upregulated hub genes ###changed to 20 after re-run
m_hub_sigs_heatmap <- m_hub_sigs_up[1:20,]

m_hub_sigs_down <- M_sigs[as.vector(m_hubs_down),]
m_hubs_down_symbols <- read.csv("m_down.csv") 
m_hubs_down_symbols <- m_hubs_down_symbols$symbol  
m_hub_sigs_down$symbol <- m_hubs_down_symbols
###top 10 downregulated hub genes ###changed to 6 after re_run
m_hub_sigs_heatmap <- rbind(m_hub_sigs_heatmap,m_hub_sigs_down[1:6,])

m_hub_sigs_heatmap <- m_hub_sigs_heatmap[order(m_hub_sigs_heatmap$log2FoldChange, decreasing = TRUE),]

###using vst normalized count data from before
#### alternatively can use rlog with this command: rlog_out <- rlog(dds2, blind = FALSE) 

###hub genes x samples
mat_m <- assay(vsdata2)[rownames(m_hub_sigs_heatmap), rownames(meta[which(meta$Subtype == c("M", "N")),])]
colnames(mat_m) <- rownames(meta[which(meta$Subtype == c("M", "N")),])
base_mean_m <- rowMeans(mat_m)
###center and scale each column (Z-score) then transpose
mat_m_scaled <- t(apply(mat_m, 1, scale))
colnames(mat_m_scaled) <- colnames(mat_m)

###getting log2foldchange for each gene
l2fc_val_m <- as.matrix(m_hub_sigs_heatmap$log2FoldChange)
colnames(l2fc_val_m) <- "logFC"

###getting mean value for each gene
mean_m <- as.matrix(m_hub_sigs_heatmap$baseMean)
colnames(mean_m) <- "AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2fc values
col_logFC_m <- colorRamp2(c(min(l2fc_val_m), 0, max(l2fc_val_m)), c("blue", "white", "red"))

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr_m <- colorRamp2(c(quantile(mean_m)[1], quantile(mean_m)[4]), c("yellow", "orange"))




ha_main_m <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2) ,
                                                  labels = c("Normal", "M"),
                                                  labels_gp = gpar(col = "black", fontsize = 20)))

ha_l2fc_m <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill =2),
                                                        height = unit(1, "cm")))

h1_m <- Heatmap(mat_m_scaled, cluster_rows = F, top_annotation = ha_main_m, show_column_names = F,
                  column_km = 2, name = "Z-score", cluster_columns = F)

h2_m <- Heatmap(l2fc_val_m, row_labels = m_hub_sigs_heatmap$symbol,
                  cluster_rows = F, name = "logFC", top_annotation = ha_l2fc_m, col = col_logFC_m,
                  cell_fun = function(j, i, x, y, w, h, col) {  #add text to each grid
                    grid.text(round(l2fc_val_m[i, j],2), x, y)
                  })

h3_m <- Heatmap(mean_m, row_labels = m_hub_sigs_heatmap$symbol,
                  cluster_rows = F, name = "AveExpr", col = col_AveExpr_m,
                  cell_fun = function(j, i, x, y, w, h, col){   #add text to each grid
                    grid.text(round(mean_m[i, j],2), x, y)
                  })

pdf("m_heatmap_100x10.pdf", width = 40, height = 55)
h_m <- h1_m + h2_m + h3_m
h_m
dev.off()

## re-run but with 20x6 genes
pdf("m_heatmap_20x6.pdf", width = 40, height = 15)
h_m <- h1_m + h2_m + h3_m
h_m
dev.off()



## LAR
lar_hub_sigs_up <- LAR_sigs[as.vector(lar_hubs_up),]
lar_hubs_up_symbols <- read.csv("lar_up.csv") 
lar_hubs_up_symbols <- lar_hubs_up_symbols$symbol  
lar_hub_sigs_up$symbol <- lar_hubs_up_symbols
###top 100 upregulated hub genes ###changed to 20 after re-run
lar_hub_sigs_heatmap <- lar_hub_sigs_up[1:20,]

lar_hub_sigs_down <- LAR_sigs[as.vector(lar_hubs_down),]
lar_hubs_down_symbols <- read.csv("lar_down.csv") 
lar_hubs_down_symbols <- lar_hubs_down_symbols$symbol  
lar_hub_sigs_down$symbol <- lar_hubs_down_symbols
###top 10 downregulated hub genes ###changed to 6 after re_run
lar_hub_sigs_heatmap <- rbind(lar_hub_sigs_heatmap,lar_hub_sigs_down[1:6,])

lar_hub_sigs_heatmap <- lar_hub_sigs_heatmap[order(lar_hub_sigs_heatmap$log2FoldChange, decreasing = TRUE),]

###using vst normalized count data from before
#### alternatively can use rlog with this command: rlog_out <- rlog(dds2, blind = FALSE) 

###hub genes x samples
mat_lar <- assay(vsdata2)[rownames(lar_hub_sigs_heatmap), rownames(meta[which(meta$Subtype == c("LAR", "N")),])]
colnames(mat_lar) <- rownames(meta[which(meta$Subtype == c("LAR", "N")),])
base_mean_lar <- rowMeans(mat_lar)
###center and scale each column (Z-score) then transpose
mat_lar_scaled <- t(apply(mat_lar, 1, scale))
colnames(mat_lar_scaled) <- colnames(mat_lar)

###getting log2foldchange for each gene
l2fc_val_lar <- as.matrix(lar_hub_sigs_heatmap$log2FoldChange)
colnames(l2fc_val_lar) <- "logFC"

###getting mean value for each gene
mean_lar <- as.matrix(lar_hub_sigs_heatmap$baseMean)
colnames(mean_lar) <- "AveExpr"

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#maps values between b/w/r for min and max l2fc values
col_logFC_lar <- colorRamp2(c(min(l2fc_val_lar), 0, max(l2fc_val_lar)), c("blue", "white", "red"))

#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr_lar <- colorRamp2(c(quantile(mean_lar)[1], quantile(mean_lar)[4]), c("yellow", "orange"))




ha_main_lar <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:2) ,
                                                  labels = c("Normal", "LAR"),
                                                  labels_gp = gpar(col = "black", fontsize = 20)))

ha_l2fc_lar <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill =2),
                                                        height = unit(1, "cm")))

h1_lar <- Heatmap(mat_lar_scaled, cluster_rows = F, top_annotation = ha_main_lar, show_column_names = F,
                  column_km = 2, name = "Z-score", cluster_columns = F)

h2_lar <- Heatmap(l2fc_val_lar, row_labels = lar_hub_sigs_heatmap$symbol,
                  cluster_rows = F, name = "logFC", top_annotation = ha_l2fc_lar, col = col_logFC_lar,
                  cell_fun = function(j, i, x, y, w, h, col) {  #add text to each grid
                    grid.text(round(l2fc_val_lar[i, j],2), x, y)
                  })

h3_lar <- Heatmap(mean_lar, row_labels = lar_hub_sigs_heatmap$symbol,
                  cluster_rows = F, name = "AveExpr", col = col_AveExpr_lar,
                  cell_fun = function(j, i, x, y, w, h, col){   #add text to each grid
                    grid.text(round(mean_lar[i, j],2), x, y)
                  })

pdf("lar_heatmap_100x10.pdf", width = 40, height = 55)
h_lar <- h1_lar + h2_lar + h3_lar
h_lar
dev.off()

## re-run but with 20x6 genes
pdf("lar_heatmap_20x6.pdf", width = 40, height = 15)
h_lar <- h1_lar + h2_lar + h3_lar
h_lar
dev.off()




# Structural analysis of the drugs
library(RColorBrewer)

## correlation matrix
library(corrplot)
dragon <- read.csv("dragon.csv")
dragon <- dragon[,-1]
colnames(dragon) <- dragon[2,]
dragon <- dragon[-c(1,2),]
rownames(dragon) <- dragon[,1]
dragon <- dragon[,-1]
### converting data type from chr to numeric
dragon <- apply(dragon, c(1,2) ,as.numeric)

cor_mat <- round(cor(t(dragon)),2)

color_data <- read.csv("colors.csv")
rownames(color_data) <- color_data[,1]

pdf("correlation_matrix.pdf", width = 45, height = 45)
corrplot(cor_mat, method = "color", type = "lower", order = "AOE", cl.ratio = 0.07,
         tl.col = color_data$color, tl.cex = 2, outline = T, addgrid.col = "darkgray", tl.srt = 45,
         addCoef.col = "white", col = colorRampPalette(c("darkred","grey","midnightblue"))(100))
dev.off()


## distance dendrogram

### Compute Euclidean distance between samples
dist_mat <- dist(dragon , diag=TRUE)
### hierarchical clustering
hc <- hclust(dist_mat)
dhc <- as.dendrogram(hc)


### setting leaf attributes
i=0
colLab<<-function(n){
  if(is.leaf(n)){
    
    ####taking the current attributes
    a=attributes(n)
    
    ####setting color.
    ligne=match(attributes(n)$label, rownames(color_data))
    color = color_data[ligne,3];
   
    ####Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5, lab.cex=1, pch=20,
                                        col=color, lab.col=color, lab.font=1, lab.cex=1))
  }
  return(n)
}

### applying leaf attributes to dendrogram
dL <- dendrapply(dhc, colLab)

### plot
pdf("euclidean_distance_dendrogram.pdf", width = 20, height = 31)
plot(dL)
dev.off()


# Drugs target enrichment analysis

drugs_target <- read.csv("drug_targets.csv")
drugs_target_entrez <- mapIds(x = org.Hs.eg.db,
                           keys = drugs_target[,1],
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")

## KEGG

drugs_target_kegg <- enrichKEGG(gene = drugs_target_entrez, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("drugs_target_kegg.pdf", height = 12)
plot(barplot(drugs_target_kegg, showCategory = 20))
dev.off()

write.csv(drugs_target_kegg, "drugs_target_kegg.csv")


## Reactome

drugs_target_reactome <- enrichPathway(gene = drugs_target_entrez, organism = "human", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("drugs_target_reactome.pdf", height = 12)
plot(barplot(drugs_target_reactome, showCategory = 20))
dev.off()
write.csv(drugs_target_reactome, "drugs_target_reactome.csv")


## Disease Ontology

drugs_target_do <- enrichDO(gene = drugs_target_entrez, ont = "DO", pAdjustMethod = "BH", pvalueCutoff = 0.05)

pdf("drugs_target_do.pdf", height = 12)
plot(barplot(drugs_target_do, showCategory = 20))
dev.off()
write.csv(drugs_target_do, "drugs_target_do.csv")


## Gene Ontology

### Biological Processes
drugs_target_bp <- enrichGO(gene = drugs_target_entrez, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID",
                    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("drugs_target_bp.pdf", height = 12)
plot(barplot(drugs_target_bp, showCategory = 20))
dev.off()
write.csv(drugs_target_bp, "drugs_target_bp.csv")

### Molecular Functions
drugs_target_mf <- enrichGO(gene = drugs_target_entrez, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID",
                    ont = "MF" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("drugs_target_mf.pdf", height = 12)
plot(barplot(drugs_target_mf, showCategory = 20))
dev.off()
write.csv(drugs_target_mf, "drugs_target_mf.csv")

### Cellular Components

drugs_target_cc <- enrichGO(gene = drugs_target_entrez, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID",
                    ont = "CC" , pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

pdf("drugs_target_cc.pdf", height = 12)
plot(barplot(drugs_target_cc, showCategory = 20))
dev.off()
write.csv(drugs_target_cc, "drugs_target_cc.csv")