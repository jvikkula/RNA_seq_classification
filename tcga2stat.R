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
brca.tum.norm <- TumorNormalMatch(rnaseq.brca$dat)
lgg.tum.norm <- TumorNormalMatch(rnaseq.lgg$dat)
coadread.tum.norm <- TumorNormalMatch(rnaseq.coadread$dat)
kipan.tum.norm <- TumorNormalMatch(rnaseq.kipan$dat)
meso.tum.norm <- TumorNormalMatch(rnaseq.meso$dat)
sarc.tum.norm <- TumorNormalMatch(rnaseq.sarc$dat)
thca.tum.norm <- TumorNormalMatch(rnaseq.thca$dat)
thym.tum.norm <- TumorNormalMatch(rnaseq.thym$dat)
ucec.tum.norm <- TumorNormalMatch(rnaseq.ucec$dat)
ucs.tum.norm <- TumorNormalMatch(rnaseq.ucs$dat)
save(brca.tum.norm, file = "brca.tum.norm.RData")
save(coadread.tum.norm, file = "coadread.tum.norm.RData")
save(kipan.tum.norm, file = "kipan.tum.norm.RData")
save(sarc.tum.norm, file = "sarc.tum.norm.RData")
save(thca.tum.norm, file = "thca.tum.norm.RData")
save(thym.tum.norm, file = "thym.tum.norm.RData")
save(ucec.tum.norm, file = "ucec.tum.norm.RData")













