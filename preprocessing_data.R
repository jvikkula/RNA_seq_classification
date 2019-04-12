library(dplyr, quietly = T)
library(tibble)
library(edgeR, quietly=T)


process <- function(df_counts, df_clinic){
  counts <- df_counts
  clinic <- df_clinic
  names_clin <- rownames(clinic)
  colnames(counts) <- substring(colnames(counts),1,12)
  matches <- match(names_clin, colnames(counts))
  matches <- matches[!is.na(matches)]
  counts <- counts[, matches]
  clinic <- clinic[which(names_clin %in% colnames(counts)), ]
  return(list(counts, clinic))
}

remove_rowzeros <- function(counts){
  trh <- counts > 0
  keep <- rowSums(trh) > 0
  counts <- counts[keep, ]
  counts
}

rowCPM_above_one <- function(counts){
  y <- DGEList(counts)
  CPM <- cpm(y, keep.lib.sizes=FALSE)
  trh <- CPM > 1
  keep <- rowSums(trh) >= 20
  y <- DGEList(y[keep, ])
  as.data.frame(y$counts)
}

log2_pseudo_counts <- function(counts){
  counts <- log2(counts+1)
  counts
}

filter_classes <- function(counts, clinical, clinical_filter){
  clin_count_bind <- rbind(counts, t(clinical))
  data <- as.data.frame(t(clin_count_bind))
  
  filtered_counts <- data %>% 
    rownames_to_column('Sample') %>% 
    filter(histologicaltype==clinical_filter) %>% #switch the filter value based on what youre filtering
    select(-seq(ncol(data)-ncol(clinical)+2, ncol(data)+1)) %>%
    column_to_rownames('Sample')
  filtered_counts <- as.data.frame(t(filtered_counts))
  filtered_counts
}

driver <- function(tcgadata, preprocess=T, log2 = F){
  full_counts <- as.data.frame(tcgadata$dat)
  clinical_full <- tcgadata$clinical
  processed <- process(full_counts, clinical_full)
  counts <- processed[[1]]
  
  if(preprocess == T){
    counts <- remove_rowzeros(counts)
    counts <- rowCPM_above_one(counts)
  }
  
  if(log2 == T){
    counts <- log2_pseudo_counts(counts)
  }
  
  clinical <- processed[[2]]
  return(list(counts, clinical))
}


#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# BRCA DATASETS
# WITH PREPROCESSING
brca_data <- driver(rnaseq.brca)
brca_counts <- brca_data[[1]]
brca_clinic <- brca_data[[2]]
head(brca_counts[, 1:20])
summary(brca_clinic)

brca1 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='infiltrating ductal carcinoma')
brca2 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='infiltrating lobular carcinoma')
brca3 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='mucinous carcinoma')

write.table(brca1, file="~/LST/brca_ductal_rsem_processed.csv", sep = "\t")
write.table(brca2, file="~/LST/brca_lobular_rsem_processed.csv", sep = "\t")
write.table(brca3, file="~/LST/brca_mucinous_rsem_processed.csv", sep = "\t")

#################################################################
# LOG2
brca_data <- driver(rnaseq.brca, log2 = T)
brca_counts <- brca_data[[1]]
brca_clinic <- brca_data[[2]]
head(brca_counts[, 1:20])
summary(brca_clinic)

brca1 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='infiltrating ductal carcinoma')
brca2 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='infiltrating lobular carcinoma')
brca3 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='mucinous carcinoma')

write.table(brca1, file="~/LST/brca_ductal_rsem_processed_log2.csv", sep = "\t")
write.table(brca2, file="~/LST/brca_lobular_rsem_processed_log2.csv", sep = "\t")
write.table(brca3, file="~/LST/brca_mucinous_rsem_processed_log2.csv", sep = "\t")
#################################################################
# NO PREPROCESSING
brca_data <- driver(rnaseq.brca, preprocess = F)
brca_counts <- brca_data[[1]]
brca_clinic <- brca_data[[2]]
head(brca_counts[, 1:20])
summary(brca_clinic)

brca1 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='infiltrating ductal carcinoma')
brca2 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='infiltrating lobular carcinoma')
brca3 <- filter_classes(counts=brca_counts, clinical=brca_clinic, clinical_filter='mucinous carcinoma')

write.table(brca1, file="~/LST/brca_ductal_rsem_unprocessed.csv", sep = "\t")
write.table(brca2, file="~/LST/brca_lobular_rsem_unprocessed.csv", sep = "\t")
write.table(brca3, file="~/LST/brca_mucinous_rsem_unprocessed.csv", sep = "\t")
#################################################################
#################################################################
#################################################################
#################################################################
#COADREAD DATASETS
# WITH PREPROCESSING
coadread_data <- driver(rnaseq.coadread)
coadread_counts <- coadread_data[[1]]
coadread_clinic <- coadread_data[[2]]
head(coadread_counts[, 1:20])
summary(coadread_clinic)

coadread1 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='colon adenocarcinoma')
coadread2 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='colon mucinous adenocarcinoma')
coadread3 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='rectal adenocarcinoma')
coadread4 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='rectal mucinous adenocarcinoma')

write.table(coadread1, file="~/LST/colon_adenocarcinoma_rsem_processed.csv", sep = "\t")
write.table(coadread2, file="~/LST/colon_mucinous_adenocarcinoma_rsem_processed.csv", sep = "\t")
write.table(coadread3, file="~/LST/rectal_adenocarcinoma_rsem_processed.csv", sep = "\t")
write.table(coadread4, file="~/LST/rectal_mucinous_adenocarcinoma_rsem_processed.csv", sep = "\t")


#################################################################
# LOG2
coadread_data <- driver(rnaseq.coadread, log2 = T)
coadread_counts <- coadread_data[[1]]
coadread_clinic <- coadread_data[[2]]
head(coadread_counts[, 1:20])
summary(coadread_clinic)

coadread1 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='colon adenocarcinoma')
coadread2 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='colon mucinous adenocarcinoma')
coadread3 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='rectal adenocarcinoma')
coadread4 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='rectal mucinous adenocarcinoma')

write.table(coadread1, file="~/LST/colon_adenocarcinoma_rsem_processed_log2.csv", sep = "\t")
write.table(coadread2, file="~/LST/colon_mucinous_adenocarcinoma_rsem_processed_log2.csv", sep = "\t")
write.table(coadread3, file="~/LST/rectal_adenocarcinoma_rsem_processed_log2.csv", sep = "\t")
write.table(coadread4, file="~/LST/rectal_mucinous_adenocarcinoma_rsem_processed_log2.csv", sep = "\t")
#################################################################
# NOPREPROCESSING
coadread_data <- driver(rnaseq.coadread, preprocess = F)
coadread_counts <- coadread_data[[1]]
coadread_clinic <- coadread_data[[2]]
head(coadread_counts[, 1:20])
summary(coadread_clinic)

coadread1 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='colon adenocarcinoma')
coadread2 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='colon mucinous adenocarcinoma')
coadread3 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='rectal adenocarcinoma')
coadread4 <- filter_classes(counts=coadread_counts, clinical=coadread_clinic, clinical_filter='rectal mucinous adenocarcinoma')

write.table(coadread1, file="~/LST/colon_adenocarcinoma_rsem_unprocessed.csv", sep = "\t")
write.table(coadread2, file="~/LST/colon_mucinous_adenocarcinoma_rsem_unprocessed.csv", sep = "\t")
write.table(coadread3, file="~/LST/rectal_adenocarcinoma_rsem_unprocessed.csv", sep = "\t")
write.table(coadread4, file="~/LST/rectal_mucinous_adenocarcinoma_rsem_unprocessed.csv", sep = "\t")
#################################################################
#################################################################
#################################################################
# KIPAN DATASETS
# WITH PREPROCESSING
kipan_data <- driver(rnaseq.kipan, preprocess = F)
kipan_counts <- kipan_data[[1]]
kipan_clinic <- kipan_data[[2]]
head(kipan_counts[, 1:20])
summary(kipan_clinic)

kipan1 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney chromophobe')
kipan2 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney clear cell renal carcinoma')
kipan3 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney papillary renal cell carcinoma')

write.table(kipan1, file="~/LST/kidney_chromophobe_rsem_processed.csv", sep = "\t")
write.table(kipan2, file="~/LST/kidney_clear_cell_renal_carcinoma_rsem_processed.csv", sep = "\t")
write.table(kipan3, file="~/LST/kidney_papillary_renal_cell_carcinoma_rsem_processed.csv", sep = "\t")
#################################################################
# LOG2
kipan_data <- driver(rnaseq.kipan, log2 = T)
kipan_counts <- kipan_data[[1]]
kipan_clinic <- kipan_data[[2]]
head(kipan_counts[, 1:20])
summary(kipan_clinic)

kipan1 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney chromophobe')
kipan2 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney clear cell renal carcinoma')
kipan3 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney papillary renal cell carcinoma')

write.table(kipan1, file="~/LST/kidney_chromophobe_rsem_processed_log2.csv", sep = "\t")
write.table(kipan2, file="~/LST/kidney_clear_cell_renal_carcinoma_rsem_processed_log2.csv", sep = "\t")
write.table(kipan3, file="~/LST/kidney_papillary_renal_cell_carcinoma_rsem_processed_log2.csv", sep = "\t")

#################################################################
# NO PREPROCESSING
kipan_data <- driver(rnaseq.kipan)
kipan_counts <- kipan_data[[1]]
kipan_clinic <- kipan_data[[2]]
head(kipan_counts[, 1:20])
summary(kipan_clinic)

kipan1 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney chromophobe')
kipan2 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney clear cell renal carcinoma')
kipan3 <- filter_classes(counts=kipan_counts, clinical=kipan_clinic, clinical_filter='kidney papillary renal cell carcinoma')

write.table(kipan1, file="~/LST/kidney_chromophobe_rsem_unprocessed.csv", sep = "\t")
write.table(kipan2, file="~/LST/kidney_clear_cell_renal_carcinoma_rsem_unprocessed.csv", sep = "\t")
write.table(kipan3, file="~/LST/kidney_papillary_renal_cell_carcinoma_rsem_unprocessed.csv", sep = "\t")
#################################################################
#################################################################
#################################################################
# LGG DATASETS
# WITH PREPROCESSING
lgg_data <- driver(rnaseq.lgg)
lgg_counts <- lgg_data[[1]]
lgg_clinic <- lgg_data[[2]]
head(lgg_counts[, 1:20])
summary(lgg_clinic)

lgg1 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='astrocytoma')
lgg2 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='oligoastrocytoma')
lgg3 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='oligodendroglioma')

write.table(lgg1, file="~/LST/astrocytoma_rsem_processed.csv", sep = "\t")
write.table(lgg2, file="~/LST/oligoastrocytoma_rsem_processed.csv", sep = "\t")
write.table(lgg3, file="~/LST/oligodendroglioma_rsem_processed.csv", sep = "\t")

#################################################################
# LOG2
lgg_data <- driver(rnaseq.lgg, log2 = T)
lgg_counts <- lgg_data[[1]]
lgg_clinic <- lgg_data[[2]]
head(lgg_counts[, 1:20])
summary(lgg_clinic)

lgg1 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='astrocytoma')
lgg2 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='oligoastrocytoma')
lgg3 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='oligodendroglioma')

write.table(lgg1, file="~/LST/astrocytoma_rsem_processed_log2.csv", sep = "\t")
write.table(lgg2, file="~/LST/oligoastrocytoma_rsem_processed_log2.csv", sep = "\t")
write.table(lgg3, file="~/LST/oligodendroglioma_rsem_processed_log2.csv", sep = "\t")

#################################################################
# NO PREPROCESSING
lgg_data <- driver(rnaseq.lgg, preprocess = F)
lgg_counts <- lgg_data[[1]]
lgg_clinic <- lgg_data[[2]]
head(lgg_counts[, 1:20])
summary(lgg_clinic)

lgg1 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='astrocytoma')
lgg2 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='oligoastrocytoma')
lgg3 <- filter_classes(counts=lgg_counts, clinical=lgg_clinic, clinical_filter='oligodendroglioma')

write.table(lgg1, file="~/LST/astrocytoma_rsem_unprocessed.csv", sep = "\t")
write.table(lgg2, file="~/LST/oligoastrocytoma_rsem_unprocessed.csv", sep = "\t")
write.table(lgg3, file="~/LST/oligodendroglioma_rsem_unprocessed.csv", sep = "\t")
#################################################################
#################################################################
#################################################################
# MESO DATASETS
# WITH PREPROCESSING
meso_data <- driver(rnaseq.meso)
meso_counts <- meso_data[[1]]
meso_clinic <- meso_data[[2]]
head(meso_counts[, 1:20])
summary(meso_clinic)

meso1 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='biphasic mesothelioma')
meso2 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='diffuse malignant mesothelioma - nos')
meso3 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='epithelioid mesothelioma')
meso4 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='sarcomatoid mesothelioma')

write.table(meso1, file="~/LST/biphasic_mesothelioma_rsem_processed.csv", sep = "\t")
write.table(meso2, file="~/LST/diffuse_malignant_mesothelioma_rsem_processed.csv", sep = "\t")
write.table(meso3, file="~/LST/epithelioid_mesothelioma_rsem_processed.csv", sep = "\t")
write.table(meso4, file="~/LST/sarcomatoid_mesothelioma_rsem_processed.csv", sep = "\t")
#################################################################
# LOG2
meso_data <- driver(rnaseq.meso, log2 = T)
meso_counts <- meso_data[[1]]
meso_clinic <- meso_data[[2]]
head(meso_counts[, 1:20])
summary(meso_clinic)

meso1 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='biphasic mesothelioma')
meso2 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='diffuse malignant mesothelioma - nos')
meso3 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='epithelioid mesothelioma')
meso4 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='sarcomatoid mesothelioma')

write.table(meso1, file="~/LST/biphasic_mesothelioma_rsem_processed_log2.csv", sep = "\t")
write.table(meso2, file="~/LST/diffuse_malignant_mesothelioma_rsem_processed_log2.csv", sep = "\t")
write.table(meso3, file="~/LST/epithelioid_mesothelioma_rsem_processed_log2.csv", sep = "\t")
write.table(meso4, file="~/LST/sarcomatoid_mesothelioma_rsem_processed_log2.csv", sep = "\t")

#################################################################
# NO PREPROCESSING
meso_data <- driver(rnaseq.meso, preprocess = F)
meso_counts <- meso_data[[1]]
meso_clinic <- meso_data[[2]]
head(meso_counts[, 1:20])
summary(meso_clinic)

meso1 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='biphasic mesothelioma')
meso2 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='diffuse malignant mesothelioma - nos')
meso3 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='epithelioid mesothelioma')
meso4 <- filter_classes(counts=meso_counts, clinical=meso_clinic, clinical_filter='sarcomatoid mesothelioma')

write.table(meso1, file="~/LST/biphasic_mesothelioma_rsem_unprocessed.csv", sep = "\t")
write.table(meso2, file="~/LST/diffuse_malignant_mesothelioma_rsem_unprocessed.csv", sep = "\t")
write.table(meso3, file="~/LST/epithelioid_mesothelioma_rsem_unprocessed.csv", sep = "\t")
write.table(meso4, file="~/LST/sarcomatoid_mesothelioma_rsem_unprocessed.csv", sep = "\t")
#################################################################
#################################################################
#################################################################
# SARKOMA DATASETS
# WITH PREPROCESSING
sarc_data <- driver(rnaseq.sarc)
sarc_counts <- sarc_data[[1]]
sarc_clinic <- sarc_data[[2]]
head(sarc_counts[, 1:20])
summary(sarc_clinic)

sarc1 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='leiomyosarcoma (lms)')
sarc2 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='dedifferentiated liposarcoma')
sarc3 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter="pleomorphic 'mfh' / undifferentiated pleomorphic sarcoma")
sarc4 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='myxofibrosarcoma')
sarc5 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='undifferentiated pleomorphic sarcoma (ups)')
sarc6 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='malignant peripheral nerve sheath tumors (mpnst)')

write.table(sarc1, file="~/LST/leiomyosarcoma_rsem_processed.csv", sep = "\t")
write.table(sarc2, file="~/LST/dedifferentiated_liposarcoma_rsem_processed.csv", sep = "\t")
write.table(sarc3, file="~/LST/pleomorphic_mfh_rsem_processed.csv", sep = "\t")
write.table(sarc4, file="~/LST/myxofibrosarcoma_rsem_processed.csv", sep = "\t")
write.table(sarc5, file="~/LST/undifferentiated_pleomorphic_sarcoma_rsem_processed.csv", sep = "\t")
write.table(sarc6, file="~/LST/malignant_peripheral_nerve_sheath_tumors_rsem_processed.csv", sep = "\t")

#################################################################
# LOG2
sarc_data <- driver(rnaseq.sarc, log2 = T)
sarc_counts <- sarc_data[[1]]
sarc_clinic <- sarc_data[[2]]
head(sarc_counts[, 1:20])
summary(sarc_clinic)

sarc1 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='leiomyosarcoma (lms)')
sarc2 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='dedifferentiated liposarcoma')
sarc3 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter="pleomorphic 'mfh' / undifferentiated pleomorphic sarcoma")
sarc4 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='myxofibrosarcoma')
sarc5 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='undifferentiated pleomorphic sarcoma (ups)')
sarc6 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='malignant peripheral nerve sheath tumors (mpnst)')

write.table(sarc1, file="~/LST/leiomyosarcoma_rsem_processed_log2.csv", sep = "\t")
write.table(sarc2, file="~/LST/dedifferentiated_liposarcoma_rsem_processed_log2.csv", sep = "\t")
write.table(sarc3, file="~/LST/pleomorphic_mfh_rsem_processed_log2.csv", sep = "\t")
write.table(sarc4, file="~/LST/myxofibrosarcoma_rsem_processed_log2.csv", sep = "\t")
write.table(sarc5, file="~/LST/undifferentiated_pleomorphic_sarcoma_rsem_processed_log2.csv", sep = "\t")
write.table(sarc6, file="~/LST/malignant_peripheral_nerve_sheath_tumors_rsem_processed_log2.csv", sep = "\t")

#################################################################
# NO PREPROCESSING
sarc_data <- driver(rnaseq.sarc, preprocess = F)
sarc_counts <- sarc_data[[1]]
sarc_clinic <- sarc_data[[2]]
head(sarc_counts[, 1:20])
summary(sarc_clinic)

sarc1 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='leiomyosarcoma (lms)')
sarc2 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='dedifferentiated liposarcoma')
sarc3 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter="pleomorphic 'mfh' / undifferentiated pleomorphic sarcoma")
sarc4 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='myxofibrosarcoma')
sarc5 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='undifferentiated pleomorphic sarcoma (ups)')
sarc6 <- filter_classes(counts=sarc_counts, clinical=sarc_clinic, clinical_filter='malignant peripheral nerve sheath tumors (mpnst)')

write.table(sarc1, file="~/LST/leiomyosarcoma_rsem_unprocessed.csv", sep = "\t")
write.table(sarc2, file="~/LST/dedifferentiated_liposarcoma_rsem_unprocessed.csv", sep = "\t")
write.table(sarc3, file="~/LST/pleomorphic_mfh_rsem_unprocessed.csv", sep = "\t")
write.table(sarc4, file="~/LST/myxofibrosarcoma_rsem_unprocessed.csv", sep = "\t")
write.table(sarc5, file="~/LST/undifferentiated_pleomorphic_sarcoma_rsem_unprocessed.csv", sep = "\t")
write.table(sarc6, file="~/LST/malignant_peripheral_nerve_sheath_tumors_rsem_unprocessed.csv", sep = "\t")
#################################################################
#################################################################
#################################################################
# THCA DATASETS
# WITH PREPROCESSING
thca_data <- driver(rnaseq.thca)
thca_counts <- thca_data[[1]]
thca_clinic <- thca_data[[2]]
head(thca_counts[, 1:20])
summary(thca_clinic)

thca1 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter='thyroid papillary carcinoma - classical/usual')
thca2 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter='thyroid papillary carcinoma - follicular (>= 99% follicular patterned)')
thca3 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter="thyroid papillary carcinoma - tall cell (>= 50% tall cell features)")

write.table(thca1, file="~/LST/thyroid_papillary_carcinoma_classical_rsem_processed.csv", sep = "\t")
write.table(thca2, file="~/LST/thyroid_papillary_carcinoma_follicular_rsem_processed.csv", sep = "\t")
write.table(thca3, file="~/LST/thyroid_papillary_carcinoma_tall_cell_rsem_processed.csv", sep = "\t")

#################################################################
# LOG2
thca_data <- driver(rnaseq.thca, log2 = T)
thca_counts <- thca_data[[1]]
thca_clinic <- thca_data[[2]]
head(thca_counts[, 1:20])
summary(thca_clinic)

thca1 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter='thyroid papillary carcinoma - classical/usual')
thca2 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter='thyroid papillary carcinoma - follicular (>= 99% follicular patterned)')
thca3 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter="thyroid papillary carcinoma - tall cell (>= 50% tall cell features)")

write.table(thca1, file="~/LST/thyroid_papillary_carcinoma_classical_rsem_processed_log2.csv", sep = "\t")
write.table(thca2, file="~/LST/thyroid_papillary_carcinoma_follicular_rsem_processed_log2.csv", sep = "\t")
write.table(thca3, file="~/LST/thyroid_papillary_carcinoma_tall_cell_rsem_processed_log2.csv", sep = "\t")

#################################################################
# NO PREPROCESSING
thca_data <- driver(rnaseq.thca, preprocess = F)
thca_counts <- thca_data[[1]]
thca_clinic <- thca_data[[2]]
head(thca_counts[, 1:20])
summary(thca_clinic)

thca1 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter='thyroid papillary carcinoma - classical/usual')
thca2 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter='thyroid papillary carcinoma - follicular (>= 99% follicular patterned)')
thca3 <- filter_classes(counts=thca_counts, clinical=thca_clinic, clinical_filter="thyroid papillary carcinoma - tall cell (>= 50% tall cell features)")

write.table(thca1, file="~/LST/thyroid_papillary_carcinoma_classical_rsem_unprocessed.csv", sep = "\t")
write.table(thca2, file="~/LST/thyroid_papillary_carcinoma_follicular_rsem_unprocessed.csv", sep = "\t")
write.table(thca3, file="~/LST/thyroid_papillary_carcinoma_tall_cell_rsem_unprocessed.csv", sep = "\t")
#################################################################
#################################################################
#################################################################
# THYM DATASETS
# WITH PREPROCESSING
thym_data <- driver(rnaseq.thym)
thym_counts <- thym_data[[1]]
thym_clinic <- thym_data[[2]]
head(thym_counts[, 1:20])
summary(thym_clinic)

thym1 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter='thymoma; type a')
thym2 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter='thymoma; type ab')
thym3 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b1")
thym4 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b2")
thym5 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b3")
thym6 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type c")

write.table(thym1, file="~/LST/thymoma_type_a_rsem_processed.csv", sep = "\t")
write.table(thym2, file="~/LST/thymoma_type_ab_rsem_processed.csv", sep = "\t")
write.table(thym3, file="~/LST/thymoma_type_b1_rsem_processed.csv", sep = "\t")
write.table(thym4, file="~/LST/thymoma_type_b2_rsem_processed.csv", sep = "\t")
write.table(thym5, file="~/LST/thymoma_type_b3_rsem_processed.csv", sep = "\t")
write.table(thym6, file="~/LST/thymoma_type_c_rsem_processed.csv", sep = "\t")

#################################################################
# LOG2
thym_data <- driver(rnaseq.thym, log2 = T)
thym_counts <- thym_data[[1]]
thym_clinic <- thym_data[[2]]
head(thym_counts[, 1:20])
summary(thym_clinic)

thym1 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter='thymoma; type a')
thym2 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter='thymoma; type ab')
thym3 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b1")
thym4 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b2")
thym5 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b3")
thym6 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type c")

write.table(thym1, file="~/LST/thymoma_type_a_rsem_processed_log2.csv", sep = "\t")
write.table(thym2, file="~/LST/thymoma_type_ab_rsem_processed_log2.csv", sep = "\t")
write.table(thym3, file="~/LST/thymoma_type_b1_rsem_processed_log2.csv", sep = "\t")
write.table(thym4, file="~/LST/thymoma_type_b2_rsem_processed_log2.csv", sep = "\t")
write.table(thym5, file="~/LST/thymoma_type_b3_rsem_processed_log2.csv", sep = "\t")
write.table(thym6, file="~/LST/thymoma_type_c_rsem_processed_log2.csv", sep = "\t")

#################################################################
# NO PREPROCESSING
thym_data <- driver(rnaseq.thym, preprocess = F)
thym_counts <- thym_data[[1]]
thym_clinic <- thym_data[[2]]
head(thym_counts[, 1:20])
summary(thym_clinic)

thym1 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter='thymoma; type a')
thym2 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter='thymoma; type ab')
thym3 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b1")
thym4 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b2")
thym5 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type b3")
thym6 <- filter_classes(counts=thym_counts, clinical=thym_clinic, clinical_filter="thymoma; type c")

write.table(thym1, file="~/LST/thymoma_type_a_rsem_unprocessed.csv", sep = "\t")
write.table(thym2, file="~/LST/thymoma_type_ab_rsem_unprocessed.csv", sep = "\t")
write.table(thym3, file="~/LST/thymoma_type_b1_rsem_unprocessed.csv", sep = "\t")
write.table(thym4, file="~/LST/thymoma_type_b2_rsem_unprocessed.csv", sep = "\t")
write.table(thym5, file="~/LST/thymoma_type_b3_rsem_unprocessed.csv", sep = "\t")
write.table(thym6, file="~/LST/thymoma_type_c_rsem_unprocessed.csv", sep = "\t")
#################################################################
#################################################################
#################################################################
# UCEC DATASETS
# WITH PREPROCESSING
ucec_data <- driver(rnaseq.ucec)
ucec_counts <- ucec_data[[1]]
ucec_clinic <- ucec_data[[2]]
head(ucec_counts[, 1:20])
summary(ucec_clinic)

ucec1 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter='endometrioid endometrial adenocarcinoma')
ucec2 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter='mixed serous and endometrioid')
ucec3 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter="serous endometrial adenocarcinoma")

write.table(ucec1, file="~/LST/endometrioid_endometrial_adenocarcinoma_rsem_processed.csv", sep = "\t")
write.table(ucec2, file="~/LST/mixed_serous_and_endometrioid_rsem_processed.csv", sep = "\t")
write.table(ucec3, file="~/LST/serous_endometrial_adenocarcinoma_rsem_processed.csv", sep = "\t")

#################################################################
# LOG2
ucec_data <- driver(rnaseq.ucec, log2 = T)
ucec_counts <- ucec_data[[1]]
ucec_clinic <- ucec_data[[2]]
head(ucec_counts[, 1:20])
summary(ucec_clinic)

ucec1 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter='endometrioid endometrial adenocarcinoma')
ucec2 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter='mixed serous and endometrioid')
ucec3 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter="serous endometrial adenocarcinoma")

write.table(ucec1, file="~/LST/endometrioid_endometrial_adenocarcinoma_rsem_processed_log2.csv", sep = "\t")
write.table(ucec2, file="~/LST/mixed_serous_and_endometrioid_rsem_processed_log2.csv", sep = "\t")
write.table(ucec3, file="~/LST/serous_endometrial_adenocarcinoma_rsem_processed_log2.csv", sep = "\t")

#################################################################
# NO PREPROCESSING
ucec_data <- driver(rnaseq.ucec, preprocess = F)
ucec_counts <- ucec_data[[1]]
ucec_clinic <- ucec_data[[2]]
head(ucec_counts[, 1:20])
summary(ucec_clinic)

ucec1 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter='endometrioid endometrial adenocarcinoma')
ucec2 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter='mixed serous and endometrioid')
ucec3 <- filter_classes(counts=ucec_counts, clinical=ucec_clinic, clinical_filter="serous endometrial adenocarcinoma")

write.table(ucec1, file="~/LST/endometrioid_endometrial_adenocarcinoma_rsem_unprocessed.csv", sep = "\t")
write.table(ucec2, file="~/LST/mixed_serous_and_endometrioid_rsem_unprocessed.csv", sep = "\t")
write.table(ucec3, file="~/LST/serous_endometrial_adenocarcinoma_rsem_unprocessed.csv", sep = "\t")

#################################################################
#################################################################
#################################################################
# UCS DATASETS
# WITH PREPROCESSING
ucs_data <- driver(rnaseq.ucs)
ucs_counts <- ucs_data[[1]]
ucs_clinic <- ucs_data[[2]]
head(ucs_counts[, 1:20])
summary(ucs_clinic)

ucs1 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter='uterine carcinosarcoma/ malignant mixed mullerian tumor (mmmt)')
ucs2 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter='uterine carcinosarcoma/ mmmt: heterologous type')
ucs3 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter="uterine carcinosarcoma/mmmt: homologous type")

write.table(ucs1, file="~/LST/malignant_mixed_mullerian_tumor_rsem_processed.csv", sep = "\t")
write.table(ucs2, file="~/LST/mmmt_heterologous_type_rsem_processed.csv", sep = "\t")
write.table(ucs3, file="~/LST/mmmt_homologous_type_rsem_processed.csv", sep = "\t")

#################################################################
# LOG2
ucs_data <- driver(rnaseq.ucs, log2 = T)
ucs_counts <- ucs_data[[1]]
ucs_clinic <- ucs_data[[2]]
head(ucs_counts[, 1:20])
summary(ucs_clinic)

ucs1 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter='uterine carcinosarcoma/ malignant mixed mullerian tumor (mmmt)')
ucs2 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter='uterine carcinosarcoma/ mmmt: heterologous type')
ucs3 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter="uterine carcinosarcoma/mmmt: homologous type")

write.table(ucs1, file="~/LST/malignant_mixed_mullerian_tumor_rsem_processed_log2.csv", sep = "\t")
write.table(ucs2, file="~/LST/mmmt_heterologous_type_rsem_processed_log2.csv", sep = "\t")
write.table(ucs3, file="~/LST/mmmt_homologous_type_rsem_processed_log2.csv", sep = "\t")

#################################################################
# NO PREPROCESSING
ucs_data <- driver(rnaseq.ucs, preprocess = F)
ucs_counts <- ucs_data[[1]]
ucs_clinic <- ucs_data[[2]]
head(ucs_counts[, 1:20])
summary(ucs_clinic)

ucs1 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter='uterine carcinosarcoma/ malignant mixed mullerian tumor (mmmt): nos')
ucs2 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter='uterine carcinosarcoma/ mmmt: heterologous type')
ucs3 <- filter_classes(counts=ucs_counts, clinical=ucs_clinic, clinical_filter="uterine carcinosarcoma/mmmt: homologous type")

write.table(ucs1, file="~/LST/malignant_mixed_mullerian_tumor_rsem_unprocessed.csv", sep = "\t")
write.table(ucs2, file="~/LST/mmmt_heterologous_type_rsem_unprocessed.csv", sep = "\t")
write.table(ucs3, file="~/LST/mmmt_homologous_type_rsem_unprocessed.csv", sep = "\t")


