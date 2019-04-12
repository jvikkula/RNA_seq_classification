load(file='~/LST/BRCA_tumor-normal.RData')
counts_tumor <- as.data.frame(rnaseq.brca.normal$primary.tumor)
counts_normal <- as.data.frame(rnaseq.brca.normal$normal)
colnames(counts_normal) <- rep(x = 'normal', length(colnames(counts_normal)))
colnames(counts_tumor) <- rep(x = 'tumor', length(colnames(counts_normal)))
count.table <- cbind(counts_tumor, counts_normal)

##########################################################
##########################################################
##########################################################
# Remove non-expressing genes 
trh <- count.table > 0
keep <- rowSums(trh) > 0
y <- count.table[keep,]

##########################################################
##########################################################
##########################################################
# More filtering based on cpm
y <- DGEList(count.table)
CPM <- cpm(y, keep.lib.sizes=FALSE)
trh <- CPM > 1
keep <- rowSums(trh) >= 1
y <- DGEList(y[keep, ])

##########################################################
##########################################################
##########################################################
# Transform to pseudocounts to ease out the computations
pseudoCounts <- log2(y$counts+1)
colnames(pseudoCounts) <- colnames(count.table)

##########################################################
##########################################################
##########################################################
# Distributions of the samples..
# boxplot(pseudoCounts,
#         xlab="Log2 counts",
#         ylab="",
#         las=2,
#         horizontal = T)
# 
# abline(v=median(pseudoCounts),col="blue")

##########################################################
##########################################################
##########################################################
# MDS plot
col.samples <- c("dark green","grey")[factor(colnames(pseudoCounts))]
plotMDS.DGEList(pseudoCounts, col = col.samples)


##########################################################
##########################################################
##########################################################
# PCA
pca <- prcomp(pseudoCounts, scale. = TRUE)
library(ggplot2)
groups <- rownames(pca$rotation)
ggplot(as.data.frame(pca$rotation), aes(x=PC1,y=PC2,col=groups)) +
  geom_text(label=groups) +
  scale_color_manual(values=c("dark green","grey"))

# Summary of first ten principal components
summary(pca)$importance[,1:10]

# Screeplot of the components
screeplot(pca, type='l')
abline(h = 1, col="red", lty=5)

# Top 500 gene loadings of PC1 and PC2
loadings <- as.data.frame(pca$x)
loadingsPC1 <- loadings[with(loadings, order(-PC1)),]
loadingsPC1 <- loadingsPC1[1:500, 1:2]

loadingsPC2 <- loadings[with(loadings, order(-PC2)),]
loadingsPC2 <- loadingsPC2[1:500, 1:2]
##########################################################
##########################################################
##########################################################
# DE analysis

# Create grouping for the samples
groups <- colnames(count.table)
# Normalize
y <- calcNormFactors(y)

#Construct a design matrix and define column names
design <- model.matrix(~0+groups, data=y$samples)
colnames(design) <- levels(factor(groups))

#Construct a contrast for comparing groups "normal" and "tumor"
TumorVsNormal <- makeContrasts(tumor-normal, levels=design)

# Estimate the abundance-dispersions
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design, method="auto")
y <- estimateGLMTagwiseDisp(y, design)

#Fit a negative binomial generalized log-linear model to the read counts for each gene.
fit <- glmFit(y, design, prior.count = 1)

#Run the likelihood ratio test using the fit object and the contrast object
lrt <- glmLRT(fit, contrast = TumorVsNormal)

# Find significant genes
toptags <- topTags(lrt, p.value = 0.05, n=Inf)
pp <- toptags$table

# Filter the original count table to contain only DE genes
ductal_brca <- read.table(file = '~/LST/data_brca_tsv/subtypes_brca_counts/ductal_brca_raw.csv')
lobular_brca <- read.table(file = '~/LST/data_brca_tsv/subtypes_brca_counts/lobular_brca_raw.csv')
ductal_brca_DE <- ductal_brca[rownames(ductal_brca) %in% rownames(pp), ]
lobular_brca_DE <- lobular_brca[rownames(lobular_brca) %in% rownames(pp), ]
write.table(ductal_brca_DE, file="~/LST/ductal_brca_DE_raw.csv", sep = "\t")
write.table(lobular_brca_DE, file="~/LST/lobular_brca_DE_raw.csv", sep = "\t")

# lol <- read.table(file = '~/LST/data_brca_tsv/subtypes_brca_DE_counts/ductal_brca_DE_raw.csv')

