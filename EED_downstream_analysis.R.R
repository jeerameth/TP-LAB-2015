## EED_downstream_analysis.R
## Jeerameth Klomsing
## Created for EED Project (Atsadang Boonmee, Tanapat Lab 2015)
## Feb 2022
## 1) Differential expression analysis from .count file
## 2) Plot volcano plots according to DE analysis results
## 3) Plot PCA from .count files
## SCRIPTS VERSION
## 15/02/22: Version 1.0
## 28/08/22: Cleaned 
###########################################################
# Set working directory
setwd('/home/atsadang/EED/jkwashere/results/aligned/counts')
directory <- '/home/atsadang/EED/jkwashere/results/aligned/counts'
# Retrieve files with .count 
sampleFiles <- grep(".count",list.files(),value=TRUE)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles)
# Extract condition & genotype from sample name
metaTable <- sampleTable
metaTable$condition <- metaTable$sampleName
metaTable$condition <- gsub("[1-6]-", "", metaTable$condition)
metaTable$condition <- gsub("-.*", "", metaTable$condition)
metaTable$genotype <- metaTable$sampleName
metaTable$genotype <- gsub(".*-", "", metaTable$genotype)
metaTable$genotype <- gsub("[1-4].*", "", metaTable$genotype)
# Analyze with DESeq2
# Take into account genotype, condition, and their interaction
library("DESeq2")
dds <- DESeqDataSetFromHTSeqCount(sampleTable = metaTable,
                                       directory = directory,
                                       design= ~condition + genotype + condition:genotype)     
# Set reference (baseline) condition
# Set to LPS100 and WT
dds$condition = relevel( dds$condition, "LPS100")
dds$genotype = relevel( dds$genotype, "WT")
# Analyze
dds <- DESeq(dds)
resultsNames(dds)
##> resultsNames(dds)
##[1] "Intercept"                  "condition_Tol_vs_LPS100"
##[3] "condition_Unstim_vs_LPS100" "genotype_KO_vs_WT"
##[5] "conditionTol.genotypeKO"         "conditionUnstim.genotypeKO"
# Set new working directory for saving files
setwd('/home/atsadang/EED/jkwashere/results/aligned')
# Get result of each pair of DE analysis and write into file
# LPS100 KO vs WT
res <- results(dds, contrast=c("genotype","KO","WT"))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/LPS100_KO_vs_WT.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
# Unstim KO vs WT
res <- results(dds, list( c("genotype_KO_vs_WT","conditionUnstim.genotypeKO")))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/Unstim_KO_vs_WT.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
# Tol KO vs WT
res <- results(dds, list( c("genotype_KO_vs_WT","conditionTol.genotypeKO")))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/Tol_KO_vs_WT.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
# Unstim vs LPS100 WT
res <- results(dds, contrast=c("condition","Unstim","LPS100"))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/Unstim_vs_LPS100_WT.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
# Unstim vs LPS100 KO
res <- results(dds, list( c("condition_Unstim_vs_LPS100","conditionUnstim.genotypeKO")))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/Unstim_vs_LPS100_KO.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
# TOL vs LPS100 WT
res <- results(dds, contrast=c("condition","Tol","LPS100"))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/Tol_vs_LPS100_WT.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
# TOL vs LPS100 KO
res <- results(dds, list( c("condition_Tol_vs_LPS100","conditionTol.genotypeKO")))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/Tol_vs_LPS100_KO.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
# Set new reference condition (for another set of analysis) and reanalyze
# Set to Unstim and WT
dds$condition = relevel( dds$condition, "Unstim")
dds$genotype = relevel( dds$genotype, "WT")
dds <- DESeq(dds)
resultsNames(dds)
## > resultsNames(dds)
## [1] "Intercept"                  "condition_LPS100_vs_Unstim"
## [3] "condition_Tol_vs_Unstim"    "genotype_KO_vs_WT"
## [5] "conditionLPS100.genotypeKO" "conditionTol.genotypeKO"
# Set new working directory for saving files
setwd('/home/atsadang/EED/jkwashere/results/aligned')
# Get result of each pair of DE analysis and write into file
# Tol vs Unstim WT
res <- results(dds, contrast=c("condition","Tol","Unstim"))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/Tol_vs_Unstim_WT.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
# Unstim vs LPS100 KO
res <- results(dds, list( c("condition_Tol_vs_Unstim","conditionTol.genotypeKO")))
res <- res[order(res$padj),]
res <- na.omit(res)
resdf <- as.data.frame(res)
write.table(resdf , file = 'DEGs/Tol_vs_Unstim_KO.txt', quote = FALSE, sep = '\t',row.names=TRUE) 
#############################
# (Optional) Plot volcano plot
# Setting parameters
FCcut = 2 # cutoff for log2 fold change
pCut = 0.01 # cutoff for p-value
# Set new working directory for retrieving file
setwd('/home/atsadang/EED/jkwashere/results/aligned/DEGs')
# Get list of DE analysis result file to loop plotting
sampleFiles <- grep(".txt",list.files(),value=TRUE)
# Import libraries
library(EnhancedVolcano)
library(grid)
# Declare plotting function
fun <- function(file) {
    df <- read.table(file,sep="\t",header=TRUE)
    df <- as.data.frame(df)
    library(data.table)
    setDT(df, keep.rownames = TRUE)[]
    df <- as.data.frame(df)
    colnames(df) <- c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
    EnhancedVolcano(df,
        lab = df$gene,  # label with gene name
        x = 'log2FoldChange',
        y = 'pvalue',
        title = gsub(".txt","",file),
        FCcutoff = FCcut,
        pCutoff = pCut,
    )
}
# Initiate PDF
pdf('volcanoplots.pdf',width=8, height=8)
# Start loop plotting through list of files
lapply(sampleFiles,fun)
# Add a page for parameter logging
grid.newpage()
grid.text(paste('FCcutoff=',FCcut,' , pCutoff=',pCut,sep=''),gp=gpar(fontsize=30))
dev.off()
#############################
# (Optional) PCA Plot
# Set new working directory for retrieving file
setwd('/home/atsadang/EED/jkwashere/results/aligned/counts')
# Get list of .count file to loop plotting
sampleFiles <- grep(".count",list.files(),value=TRUE)
# Merging count file into single table and annotate with file name
fun <- function(file) {
    dfa <- read.table(file,sep="\t",header=FALSE)
    dfa <- as.data.frame(dfa)
    colnames(dfa) <- c("gene",gsub(".count","",file))
    if(file == sampleFiles[1]){
        df <<- dfa
    } else {
        df <<- merge(df,dfa,by.x='gene',by.y='gene')
    }   
}
lapply(sampleFiles,fun)
# Write merged count table into file
write.table(df , file = 'mergedcount.txt', quote = FALSE, sep = '\t',row.names=FALSE) 
# Prepare merged table for PC analysis(transverse,filter,extract and annotate group)
df <- na.omit(df)
n <- df$gene
dft <- as.data.frame(t(df[,-1]))
colnames(dft) <- n
dft <- dft[,colSums(dft)>0]
dft$condition <- rownames(dft)
dft$condition <- gsub("[1-6]-", "", dft$condition)
dft$condition <- gsub("-.*", "", dft$condition)
dft$genotype <- rownames(dft)
dft$genotype <- gsub(".*-", "", dft$genotype)
dft$genotype <- gsub("[1-4].*", "", dft$genotype)
dft$group <- paste(dft$condition,'-',dft$genotype,sep='')
# PC analysis
library(ggfortify)
dfpca <- dft[,1:(ncol(dft)-3)]
pca_res <- prcomp(dfpca, scale. = TRUE)
# Initiate PDF
pdf('PCA.pdf',width=5, height=3)
# Plot PCA with different label targets
print(autoplot(pca_res, data = dft, colour = 'condition'))
print(autoplot(pca_res, data = dft, colour = 'genotype'))
print(autoplot(pca_res, data = dft, colour = 'group'))
print(autoplot(pca_res, data = dft, colour = "genotype", shape="condition")+
  scale_colour_manual(values=c("red","blue"))+
  scale_fill_manual(values=c("red","blue"))+
  scale_shape_manual(values=c(25,22,23)))
dev.off()