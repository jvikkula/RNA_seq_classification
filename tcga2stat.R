library(TCGA2STAT)

###############################################################################
# Get RSEM normalized data from TCGA
rnaseq.brca <- getTCGA(disease="BRCA", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.brca$clinical)

rnaseq.lgg <- getTCGA(disease="LGG", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.lgg$clinical)

rnaseq.coadread <- getTCGA(disease="COADREAD", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.coadread$clinical)

rnaseq.kipan <- getTCGA(disease="KIPAN", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.kipan$clinical)

rnaseq.meso <- getTCGA(disease="MESO", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.meso$clinical)

rnaseq.sarc <- getTCGA(disease="SARC", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.sarc$clinical)

rnaseq.thca <- getTCGA(disease="THCA", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.thca$clinical)

rnaseq.thym <- getTCGA(disease="THYM", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.thym$clinical)

rnaseq.ucec <- getTCGA(disease="UCEC", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.ucec$clinical)

rnaseq.ucs <- getTCGA(disease="UCS", data.type="RNASeq2", type="count", clinical=TRUE)
summary(rnaseq.ucs$clinical)

###############################################################################
# Get RSEM normalized tumor-normal sample pairs data from TCGA

















