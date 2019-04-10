
#deseq2-analysis-template.R

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")
BiocManager::install("XML", version = "3.8")
library(XML)
library(DESeq2)
install.packages("RColorBrewer")
library(RColorBrewer)
biocLite()
source("https://bioconductor.org/biocLite.R")
biocLite("data.table")


#deseq2-analysis-template.R



load(file='~/koulu/LST project/BRCA_tumor-normal.RData')
counts_tumor <- as.data.frame(rnaseq.brca.normal$primary.tumor)
counts_normal <- as.data.frame(rnaseq.brca.normal$normal)
colnames(counts_tumor) <- paste("tumor", colnames(counts_tumor), sep = "_")
colnames(counts_normal) <- paste("normal", colnames(counts_normal), sep = "_")
colnames(counts_normal) <- rep(x = 'normal', length(colnames(counts_normal)))
colnames(counts_tumor) <- rep(x = 'tumor', length(colnames(counts_normal)))
countdata <- cbind(counts_tumor, counts_normal)

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition (first 97 are tumor, second 97 contain the normals)
(condition <- factor(c(rep("tumor", 97), rep("normal", 97))))

# Analysis with DESeq2 ----------------------------------------------------

# Create a coldata frame and instantiate the DESeqDataSet. 
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds)
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
# rld <- rlogTransformation(dds)
rld <- vst(dds)
head(assay(rld))
hist(assay(rld))
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
plotPCA(rld, intgroup="condition")


png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()


# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

## Order by adjusted p-value
res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")
dev.off()
## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", 
     xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
plotMA(dds, ylim=c(-1,1))
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()

Results <- read.csv("diffexpr-results.csv")
install.packages("rlist")
library(rlist)

