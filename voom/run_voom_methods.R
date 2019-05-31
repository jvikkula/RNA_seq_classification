library(useful)
library(ggplot2)
library(dplyr)
library(limma)
library(caret)
library(pamr)

###############################################################
# BRCA
###############################################################

ductal <- read.csv("rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/brca_ductal_rsem_processed.csv",
                   sep = "\t")
lobular <- read.csv("rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/brca_lobular_rsem_processed.csv",
                    sep = "\t")

brca_sets <- list()
brca_sets[[1]] <- ductal
brca_sets[[2]] <- lobular

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  brca <- multiclass_preprocess(brca_sets, c("Ductal", "Lobular"))
  if (k == 1){
    folds <- do_folds(brca$sets, brca$data)
  }
  brca_params <- multiclass_train(brca$sets, brca$data, folds)
  priors <- c((782/985), (205/985))
  brca_preds <- multiclass_predict(brca_params$testdata, brca$sets, brca_params$zs,
                                   brca_params$ss, priors, c("Ductal", "Lobular"))
  conf <- confusionMatrix(factor(brca_preds$preds[1,]), factor(brca_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, conf$byClass[[1]])
  spef <- c(spef, conf$byClass[[2]])
  
}
mean(accs) 
mean(sens)
mean(spef)

###############################################################
# KIPAN
###############################################################

chrom <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_chromophobe_rsem_processed.csv",
                  sep = "\t")
clear <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_clear_cell_renal_carcinoma_rsem_processed.csv",
                  sep = "\t")


kipan_sets <- list()
kipan_sets[[1]] <- chrom
kipan_sets[[2]] <- clear

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  kipan <- multiclass_preprocess(kipan_sets, c("Chrom", "Clear", "Papil"))
  if (k == 1){
    folds <- do_folds(kipan$sets, kipan$data)
  }
  kipan_params <- multiclass_train(kipan$sets, kipan$data, folds)
  priors <- c(66/889, 533/889, 290/889)
  kipan_preds <- multiclass_predict(kipan_params$testdata, kipan$sets, kipan_params$zs,
                                    kipan_params$ss, priors, c("Chrom", "Clear", "Papil"))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, conf$byClass[[1]])
  spef <- c(spef, conf$byClass[[2]]) 
  
}
mean(accs)
mean(sens)
mean(spef)


###############################################################
# MESO
###############################################################


biphasic <- read.csv("rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/biphasic_mesothelioma_rsem_processed.csv",
                     sep = "\t")
epitheloid <- read.csv("rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/epithelioid_mesothelioma_rsem_processed.csv",
                       sep = "\t")

meso_sets <- list()
meso_sets[[1]] <- biphasic
meso_sets[[2]] <- epitheloid

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  meso <- multiclass_preprocess(meso_sets, c("Biphasic", "Epitheloid"))
  if (k == 1){
    folds <- do_folds(meso$sets, meso$data)
  }
  meso_params <- multiclass_train(meso$sets, meso$data, folds)
  priors <- c(23/80, 57/80)
  meso_preds <- multiclass_predict(meso_params$testdata, meso$sets, meso_params$zs,
                                   meso_params$ss, priors, c("Biphasic", "Epitheloid"))
  conf <- confusionMatrix(factor(meso_preds$preds[1,]), factor(meso_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, conf$byClass[[1]])
  spef <- c(spef, conf$byClass[[2]])  
  
}
mean(accs)
mean(sens)
mean(spef)


###############################################################
# UCEC
###############################################################


endo <- read.csv("rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/endometrioid_endometrial_adenocarcinoma_rsem_processed.csv",
                 sep = "\t")
mixed <- read.csv("rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/mixed_serous_and_endometrioid_rsem_processed.csv",
                  sep = "\t")
serous <- read.csv("rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/serous_endometrial_adenocarcinoma_rsem_processed.csv",
                   sep = "\t")

ucec_sets <- list()
ucec_sets[[1]] <- endo
ucec_sets[[2]] <- mixed
ucec_sets[[3]] <- serous

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){

  ucec <- multiclass_preprocess(ucec_sets, c("Endo", "Mixed", "Serous"))
  if (k == 1){
    folds <- do_folds(ucec$sets, ucec$data)
  }
  ucec_params <- multiclass_train(ucec$sets, ucec$data, folds)
  priors <- c(113/183, 12/183, 58/183)
  ucec_preds <- multiclass_predict(ucec_params$testdata, ucec$sets, ucec_params$zs,
                                   ucec_params$ss, priors, c("Endo", "Mixed", "Serous"))
  conf <- confusionMatrix(factor(ucec_preds$preds[1,]), factor(ucec_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))  
  
}

mean(accs)
mean(sens)
mean(spef)

###############################################################
# COADREAD
###############################################################


colon_ade <- read.csv("rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/colon_adenocarcinoma_rsem_processed.csv",
                      sep = "\t")
colon_muc <- read.csv("rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/colon_mucinous_adenocarcinoma_rsem_processed.csv",
                      sep = "\t")
rec_ade <- read.csv("rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/rectal_adenocarcinoma_rsem_processed.csv",
                    sep = "\t")
rec_muc <- read.csv("rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/rectal_mucinous_adenocarcinoma_rsem_processed.csv",
                    sep = "\t")


coad_sets <- list()
coad_sets[[1]] <- colon_ade
coad_sets[[2]] <- colon_muc
coad_sets[[3]] <- rec_ade
coad_sets[[4]] <- rec_muc

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){

  coad <- multiclass_preprocess(coad_sets, c("Col ade", "Col muc", 
                                             "Rec ade", "Rec muc"))
  if (k == 1){
    folds <- do_folds(coad$sets, coad$data)
  }
  coad_params <- multiclass_train(coad$sets, coad$data, folds)
  priors <- c(257/390, 38/390, 89/390, 6/390)
  coad_preds <- multiclass_predict(coad_params$testdata, coad$sets, coad_params$zs,
                                   coad_params$ss, priors, c("Col ade", "Col muc", 
                                                             "Rec ade", "Rec muc"))
  conf <- confusionMatrix(factor(coad_preds$preds[1,]), factor(coad_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))
  
}

mean(accs)
mean(sens)
mean(spef)

###############################################################
# THCA
###############################################################

classical <- read.csv("rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_classical_rsem_processed.csv",
                      sep = "\t")
follicular <- read.csv("rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_follicular_rsem_processed.csv",
                       sep = "\t")
tall <- read.csv("rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_tall_cell_rsem_processed.csv",
                 sep = "\t")

thca_sets <- list()
thca_sets[[1]] <- classical
thca_sets[[2]] <- follicular
thca_sets[[3]] <- tall

accs <- c()
sens <- c()
spef <- c()

for(k in 1:5){

  thca <- multiclass_preprocess(thca_sets, c("Classical", "Follicular", "Tall"))
  if (k == 1){
    folds <- do_folds(thca$sets, thca$data)
  }
  thca_params <- multiclass_train(thca$sets, thca$data, folds)
  priors <- c(357/512, 102/512, 35/512)
  thca_preds <- multiclass_predict(thca_params$testdata, thca$sets, thca_params$zs,
                                   thca_params$ss, priors, c("Classical", "Follicular", "Tall"))
  conf <- confusionMatrix(factor(thca_preds$preds[1,]), factor(thca_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))
  
}

mean(accs)
mean(sens)
mean(spef)

###############################################################
# LGG
###############################################################

astro <- read.csv("rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/astrocytoma_rsem_processed.csv",
                  sep = "\t")
oliastro <- read.csv("rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/oligoastrocytoma_rsem_processed.csv",
                     sep = "\t")
olidendro <- read.csv("rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/oligodendroglioma_rsem_processed.csv",
                      sep = "\t")

lgg_sets <- list()
lgg_sets[[1]] <- astro
lgg_sets[[2]] <- oliastro
lgg_sets[[3]] <- olidendro

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){

  lgg <- multiclass_preprocess(lgg_sets, c("Astro", "Oliastro", "Olidendro"))
  if (k == 1){
    folds <- do_folds(lgg$sets, lgg$data)
  }
  lgg_params <- multiclass_train(lgg$sets, lgg$data, folds)
  priors <- c(194/515, 130/515, 191/515)
  lgg_preds <- multiclass_predict(lgg_params$testdata, lgg$sets, lgg_params$zs,
                                  lgg_params$ss, priors, c("Astro", "Oliastro", "Olidendro"))
  conf <- confusionMatrix(factor(lgg_preds$preds[1,]), factor(lgg_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))
  
}
mean(accs)
mean(sens)
mean(spef)


###############################################################
# KIPAN ALL CLASSES
###############################################################

chrom <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_chromophobe_rsem_processed.csv",
                  sep = "\t")
clear <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_clear_cell_renal_carcinoma_rsem_processed.csv",
                  sep = "\t")
papil <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_papillary_renal_cell_carcinoma_rsem_processed.csv",
                  sep = "\t")

kipan_sets <- list()
kipan_sets[[1]] <- chrom
kipan_sets[[2]] <- clear
kipan_sets[[3]] <- papil

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){

  kipan <- multiclass_preprocess(kipan_sets, c("Chrom", "Clear", "Papil"))
  if (k == 1){
    folds <- do_folds(kipan$sets, kipan$data)
  }
  kipan_params <- multiclass_train(kipan$sets, kipan$data, folds)
  priors <- c(66/889, 533/889, 290/889)
  kipan_preds <- multiclass_predict(kipan_params$testdata, kipan$sets, kipan_params$zs,
                                    kipan_params$ss, priors, c("Chrom", "Clear", "Papil"))
  conf <- confusionMatrix(factor(kipan_preds$preds[1,]), factor(kipan_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))
  
}
mean(accs)
mean(sens)
mean(spef)

###############################################################
# SARCOMA
###############################################################

dedi <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/dedifferentiated_liposarcoma_rsem_processed.csv",
                 sep = "\t")
leio <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/leiomyosarcoma_rsem_processed.csv",
                 sep = "\t")
malig <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/malignant_peripheral_nerve_sheath_tumors_rsem_processed.csv",
                  sep = "\t")
myxo <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/myxofibrosarcoma_rsem_processed.csv",
                 sep = "\t")
pleo <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/pleomorphic_mfh_rsem_processed.csv",
                 sep = "\t")
undif <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/undifferentiated_pleomorphic_sarcoma_rsem_processed.csv",
                  sep = "\t")

sarc_sets <- list()
sarc_sets[[1]] <- dedi
sarc_sets[[2]] <- leio
sarc_sets[[3]] <- malig
sarc_sets[[4]] <- myxo
sarc_sets[[5]] <- pleo
sarc_sets[[6]] <- undif

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){

  sarc <- multiclass_preprocess(sarc_sets, c("Dedi", "Leio", "Malig",
                                             "Myxo", "Pleo", "Undif"))
  if (k == 1){
    folds <- do_folds(sarc$sets, sarc$data)
  }
  sarc_params <- multiclass_train(sarc$sets, sarc$data, folds)
  priors <- c(58/246, 104/246, 9/206, 25/206, 29/206, 21/246)
  sarc_preds <- multiclass_predict(sarc_params$testdata, sarc$sets, sarc_params$zs,
                                   sarc_params$ss, priors, c("Dedi", "Leio", "Malig",
                                                             "Myxo", "Pleo", "Undif"))
  conf <- confusionMatrix(factor(sarc_preds$preds[1,]), factor(sarc_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2])) 
  
}
mean(accs)
mean(sens)
mean(spef)

###############################################################
# UCS
###############################################################

malig <- read.csv("rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/malignant_mixed_mullerian_tumor_rsem_processed.csv",
                  sep = "\t")
mheter <- read.csv("rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/mmmt_heterologous_type_rsem_processed.csv",
                   sep = "\t")
mhom <- read.csv("rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/mmmt_homologous_type_rsem_processed.csv",
                 sep = "\t")


ucs_sets <- list()
ucs_sets[[1]] <- malig
ucs_sets[[2]] <- mheter
ucs_sets[[3]] <- mhom

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){

  ucs <- multiclass_preprocess(ucs_sets, c("Malig", "Mheter", "Mhom"))
  if (k == 1){
    folds <- do_folds(ucs$sets, ucs$data)
  }
  ucs_params <- multiclass_train(ucs$sets, ucs$data, folds)
  priors <- c(24/57, 20/57, 13/57)
  ucs_preds <- multiclass_predict(ucs_params$testdata, ucs$sets, ucs_params$zs,
                                  ucs_params$ss, priors, c("Malig", "Mheter", "Mhom"))
  conf <- confusionMatrix(factor(ucs_preds$preds[1,]), factor(ucs_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2])) 
  
}
mean(accs)
mean(sens)
mean(spef)

###############################################################
# THYMOMA
###############################################################

typea <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_a_rsem_processed.csv",
                  sep = "\t")
typeab <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_ab_rsem_processed.csv",
                   sep = "\t")
typeb1 <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b1_rsem_processed.csv",
                   sep = "\t")
typeb2 <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b2_rsem_processed.csv",
                   sep = "\t")
typeb3 <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b3_rsem_processed.csv",
                   sep = "\t")
typec <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_c_rsem_processed.csv",
                  sep = "\t")

thym_sets <- list()
thym_sets[[1]] <- typea
thym_sets[[2]] <- typeab
thym_sets[[3]] <- typeb1
thym_sets[[4]] <- typeb2
thym_sets[[5]] <- typeb3
thym_sets[[6]] <- typec

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){

  thym <- multiclass_preprocess(thym_sets, c("Typea", "Typeab", "Typeb1",
                                             "Typeb2", "Typeb3", "Typec"))
  if (k == 1){
    folds <- do_folds(thym$sets, thym$data)
  }
  thym_params <- multiclass_train(thym$sets, thym$data, folds)
  priors <- c(17/120, 35/120, 15/120, 31/120, 11/120, 11/120)
  thym_preds <- multiclass_predict(thym_params$testdata, thym$sets, thym_params$zs,
                                   thym_params$ss, priors, c("Typea", "Typeab", "Typeb1",
                                                             "Typeb2", "Typeb3", "Typec"))
  conf <- confusionMatrix(factor(thym_preds$preds[1,]), factor(thym_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2])) 
  
}
mean(accs)
mean(sens)
mean(spef)

###############################################################
# BRCA
###############################################################

ductal <- read.csv("rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/brca_ductal_rsem_processed.csv",
                   sep = "\t")
lobular <- read.csv("rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/brca_lobular_rsem_processed.csv",
                    sep = "\t")

brca_sets <- list()
brca_sets[[1]] <- ductal
brca_sets[[2]] <- lobular

brca_data <- cbind(ductal, lobular)
nbrca <- list()
nbrca$x <- as.matrix(brca_data)
nbrca$y <- brca$classes[,1]
brca_nsc <- pamr.train(nbrca)
pamr.geneplot(brca_nsc, nbrca, threshold = 12.718)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  brca <- multiclass_preprocess(brca_sets, c("Ductal", "Lobular"), c(2955, 1459, 3249))
  if (k == 1){
    folds <- do_folds(brca$sets, brca$data)
  }
  brca_params <- multiclass_train(brca$sets, brca$data, folds)
  priors <- c((782/985), (205/985))
  brca_preds <- multiclass_predict(brca_params$testdata, brca$sets, brca_params$zs,
                                   brca_params$ss, priors, c("Ductal", "Lobular"))
  conf <- confusionMatrix(factor(brca_preds$preds[1,]), factor(brca_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, conf$byClass[[1]])
  spef <- c(spef, conf$byClass[[2]])
  
}
mean(accs)
mean(sens)
mean(spef)

###############################################################
# KIPAN
###############################################################

chrom <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_chromophobe_rsem_processed.csv",
                  sep = "\t")
clear <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_clear_cell_renal_carcinoma_rsem_processed.csv",
                  sep = "\t")

kipan_sets <- list()
kipan_sets[[1]] <- chrom
kipan_sets[[2]] <- clear

kipan_data <- cbind(chrom, clear)
nkipan <- list()
nkipan$x <- as.matrix(kipan_data)
nkipan$y <- kipan$classes[,1]
kipan_nsc <- pamr.train(nkipan)
pamr.geneplot(kipan_nsc, nkipan, threshold = 19.921)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  kipan <- multiclass_preprocess(kipan_sets, c("Chrom", "Clear", "Papil"))
  if (k == 1){
    folds <- do_folds(kipan$sets, kipan$data)
  }
  kipan_params <- multiclass_train(kipan$sets, kipan$data, folds)
  priors <- c(66/889, 533/889, 290/889)
  kipan_preds <- multiclass_predict(kipan_params$testdata, kipan$sets, kipan_params$zs,
                                    kipan_params$ss, priors, c("Chrom", "Clear", "Papil"))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, conf$byClass[[1]])
  spef <- c(spef, conf$byClass[[2]]) 
  
}
mean(accs)
mean(sens)
mean(spef)


###############################################################
# MESO
###############################################################


biphasic <- read.csv("rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/biphasic_mesothelioma_rsem_processed.csv",
                     sep = "\t")
epitheloid <- read.csv("rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/epithelioid_mesothelioma_rsem_processed.csv",
                       sep = "\t")

meso_sets <- list()
meso_sets[[1]] <- biphasic
meso_sets[[2]] <- epitheloid

meso_data <- cbind(biphasic, epitheloid)
nmeso <- list()
nmeso$x <- as.matrix(meso_data)
nmeso$y <- meso$classes[,1]
meso_nsc <- pamr.train(nmeso)
pamr.geneplot(meso_nsc, nmeso, threshold = 4.309)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  meso <- multiclass_preprocess(meso_sets, c("Biphasic", "Epitheloid"),
                                c(3106, 2898, 13434, 3243, 5988, 3096, 8210))
  if (k == 1){
    folds <- do_folds(meso$sets, meso$data)
  }
  meso_params <- multiclass_train(meso$sets, meso$data, folds)
  priors <- c(23/80, 57/80)
  meso_preds <- multiclass_predict(meso_params$testdata, meso$sets, meso_params$zs,
                                   meso_params$ss, priors, c("Biphasic", "Epitheloid"))
  conf <- confusionMatrix(factor(meso_preds$preds[1,]), factor(meso_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, conf$byClass[[1]])
  spef <- c(spef, conf$byClass[[2]])  
  
}
mean(accs)
mean(sens)
mean(spef)



###############################################################
# UCEC
###############################################################


endo <- read.csv("rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/endometrioid_endometrial_adenocarcinoma_rsem_processed.csv",
                 sep = "\t")
mixed <- read.csv("rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/mixed_serous_and_endometrioid_rsem_processed.csv",
                  sep = "\t")
serous <- read.csv("rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/serous_endometrial_adenocarcinoma_rsem_processed.csv",
                   sep = "\t")

ucec_sets <- list()
ucec_sets[[1]] <- endo
ucec_sets[[2]] <- mixed
ucec_sets[[3]] <- serous

ucec_data <- cbind(endo, mixed, serous)
nucec <- list()
nucec$x <- as.matrix(ucec_data)
nucec$y <- ucec$classes[,1]
ucec_nsc <- pamr.train(nucec)
pamr.geneplot(ucec_nsc, nucec, threshold = 6.060)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  ucec <- multiclass_preprocess(ucec_sets, c("Endo", "Mixed", "Serous"),
                                c(2712, 7466, 3121, 10945, 15285))
  if (k == 1){
    folds <- do_folds(ucec$sets, ucec$data)
  }
  ucec_params <- multiclass_train(ucec$sets, ucec$data, folds)
  priors <- c(113/183, 12/183, 58/183)
  ucec_preds <- multiclass_predict(ucec_params$testdata, ucec$sets, ucec_params$zs,
                                   ucec_params$ss, priors, c("Endo", "Mixed", "Serous"))
  conf <- confusionMatrix(factor(ucec_preds$preds[1,]), factor(ucec_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))  
  
}

mean(accs)
mean(sens)
mean(spef)

###############################################################
# COADREAD
###############################################################


colon_ade <- read.csv("rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/colon_adenocarcinoma_rsem_processed.csv",
                      sep = "\t")
colon_muc <- read.csv("rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/colon_mucinous_adenocarcinoma_rsem_processed.csv",
                      sep = "\t")
rec_ade <- read.csv("rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/rectal_adenocarcinoma_rsem_processed.csv",
                    sep = "\t")
rec_muc <- read.csv("rnaseqV2data/processed_rnaseqV2/coadread/coadread_rsem_processed/rectal_mucinous_adenocarcinoma_rsem_processed.csv",
                    sep = "\t")

coad_sets <- list()
coad_sets[[1]] <- colon_ade
coad_sets[[2]] <- colon_muc
coad_sets[[3]] <- rec_ade
coad_sets[[4]] <- rec_muc

coad_data <- cbind(colon_ade, colon_muc, rec_ade, rec_muc)
ncoad <- list()
ncoad$x <- as.matrix(coad_data)
ncoad$y <- coad$classes[,1]
coad_nsc <- pamr.train(ncoad)
pamr.geneplot(coad_nsc, ncoad, threshold = 6.178)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  coad <- multiclass_preprocess(coad_sets, c("Col ade", "Col muc", 
                                             "Rec ade", "Rec muc"),
                                c(13084, 8861, 3012))
  if (k == 1){
    folds <- do_folds(coad$sets, coad$data)
  }
  coad_params <- multiclass_train(coad$sets, coad$data, folds)
  priors <- c(257/390, 38/390, 89/390, 6/390)
  coad_preds <- multiclass_predict(coad_params$testdata, coad$sets, coad_params$zs,
                                   coad_params$ss, priors, c("Col ade", "Col muc", 
                                                             "Rec ade", "Rec muc"))
  conf <- confusionMatrix(factor(coad_preds$preds[1,]), factor(coad_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))
  
}

mean(accs)
mean(sens)
mean(spef)

###############################################################
# THCA
###############################################################

classical <- read.csv("rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_classical_rsem_processed.csv",
                      sep = "\t")
follicular <- read.csv("rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_follicular_rsem_processed.csv",
                       sep = "\t")
tall <- read.csv("rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_tall_cell_rsem_processed.csv",
                 sep = "\t")

thca_sets <- list()
thca_sets[[1]] <- classical
thca_sets[[2]] <- follicular
thca_sets[[3]] <- tall

thca_data <- cbind(classical, follicular, tall)
nthca <- list()
nthca$x <- as.matrix(thca_data)
nthca$y <- thca$classes[,1]
thca_nsc <- pamr.train(nthca)
pamr.geneplot(thca_nsc, nthca, threshold = 9.280)

accs <- c()
sens <- c()
spef <- c()

for(k in 1:5){
  
  thca <- multiclass_preprocess(thca_sets, c("Classical", "Follicular", "Tall"),
                                c(8491, 610, 1606, 608, 7635, 6736, 10940, 4551,
                                  5654))
  if (k == 1){
    folds <- do_folds(thca$sets, thca$data)
  }
  thca_params <- multiclass_train(thca$sets, thca$data, folds)
  priors <- c(357/512, 102/512, 35/512)
  thca_preds <- multiclass_predict(thca_params$testdata, thca$sets, thca_params$zs,
                                   thca_params$ss, priors, c("Classical", "Follicular", "Tall"))
  conf <- confusionMatrix(factor(thca_preds$preds[1,]), factor(thca_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))
  
}

mean(accs)
mean(sens)
mean(spef)

###############################################################
# LGG
###############################################################

astro <- read.csv("rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/astrocytoma_rsem_processed.csv",
                  sep = "\t")
oliastro <- read.csv("rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/oligoastrocytoma_rsem_processed.csv",
                     sep = "\t")
olidendro <- read.csv("rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/oligodendroglioma_rsem_processed.csv",
                      sep = "\t")

lgg_sets <- list()
lgg_sets[[1]] <- astro
lgg_sets[[2]] <- oliastro
lgg_sets[[3]] <- olidendro

lgg_data <- cbind(astro, oliastro, olidendro)
nlgg <- list()
nlgg$x <- as.matrix(lgg_data)
nlgg$y <- lgg$classes[,1]
lgg_nsc <- pamr.train(nlgg)
pamr.geneplot(lgg_nsc, nlgg, threshold = 9.473)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  lgg <- multiclass_preprocess(lgg_sets, c("Astro", "Oliastro", "Olidendro"),
                               c(2961, 936, 15339, 2360, 15555, 3466, 12511,
                                 11753, 15223, 9418, 6394, 2339, 7474, 6791, 5495,
                                 3183, 12322, 4840, 14010, 12549, 7258))
  if (k == 1){
    folds <- do_folds(lgg$sets, lgg$data)
  }
  lgg_params <- multiclass_train(lgg$sets, lgg$data, folds)
  priors <- c(194/515, 130/515, 191/515)
  lgg_preds <- multiclass_predict(lgg_params$testdata, lgg$sets, lgg_params$zs,
                                  lgg_params$ss, priors, c("Astro", "Oliastro", "Olidendro"))
  conf <- confusionMatrix(factor(lgg_preds$preds[1,]), factor(lgg_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))
  
}
mean(accs)
mean(sens)
mean(spef)


###############################################################
# KIPAN ALL CLASSES
###############################################################

chrom <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_chromophobe_rsem_processed.csv",
                  sep = "\t")
clear <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_clear_cell_renal_carcinoma_rsem_processed.csv",
                  sep = "\t")
papil <- read.csv("rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_papillary_renal_cell_carcinoma_rsem_processed.csv",
                  sep = "\t")

kipan_sets <- list()
kipan_sets[[1]] <- chrom
kipan_sets[[2]] <- clear
kipan_sets[[3]] <- papil

kipan_data <- cbind(chrom, clear, papil)
kipan_data <- kipan_data[subset,]
nkipan <- list()
nkipan$x <- as.matrix(kipan_data)
nkipan$y <- kipan$classes[,1]
kipan_nsc <- pamr.train(nkipan)
pamr.geneplot(kipan_nsc, nkipan, threshold = 23.499)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  kipan <- multiclass_preprocess(kipan_sets, c("Chrom", "Clear", "Papil"),
                                 c(3159, 3971, 4660, 8199, 5487, 4673, 688, 6707, 5848, 2259))
  if (k == 1){
    folds <- do_folds(kipan$sets, kipan$data)
  }
  kipan_params <- multiclass_train(kipan$sets, kipan$data, folds)
  priors <- c(66/889, 533/889, 290/889)
  kipan_preds <- multiclass_predict(kipan_params$testdata, kipan$sets, kipan_params$zs,
                                    kipan_params$ss, priors, c("Chrom", "Clear", "Papil"))
  conf <- confusionMatrix(factor(kipan_preds$preds[1,]), factor(kipan_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2]))
  
}
conf
mean(accs)
mean(sens)
mean(spef)

###############################################################
# SARCOMA
###############################################################

dedi <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/dedifferentiated_liposarcoma_rsem_processed.csv",
                 sep = "\t")
leio <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/leiomyosarcoma_rsem_processed.csv",
                 sep = "\t")
malig <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/malignant_peripheral_nerve_sheath_tumors_rsem_processed.csv",
                  sep = "\t")
myxo <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/myxofibrosarcoma_rsem_processed.csv",
                 sep = "\t")
pleo <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/pleomorphic_mfh_rsem_processed.csv",
                 sep = "\t")
undif <- read.csv("rnaseqV2data/processed_rnaseqV2/sarc/sarc_rsem_processed/undifferentiated_pleomorphic_sarcoma_rsem_processed.csv",
                  sep = "\t")

sarc_sets <- list()
sarc_sets[[1]] <- dedi
sarc_sets[[2]] <- leio
sarc_sets[[3]] <- malig
sarc_sets[[4]] <- myxo
sarc_sets[[5]] <- pleo
sarc_sets[[6]] <- undif

sarc_data <- cbind(dedi, leio, malig, myxo, pleo, undif)
nsarc <- list()
nsarc$x <- as.matrix(sarc_data)
nsarc$y <- sarc$classes[,1]
sarc_nsc <- pamr.train(nsarc)
pamr.geneplot(sarc_nsc, nsarc, threshold = 12.667)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  sarc <- multiclass_preprocess(sarc_sets, c("Dedi", "Leio", "Malig",
                                             "Myxo", "Pleo", "Undif"),
                                c(8441, 14609, 6961, 5261))
  if (k == 1){
    folds <- do_folds(sarc$sets, sarc$data)
  }
  sarc_params <- multiclass_train(sarc$sets, sarc$data, folds)
  priors <- c(58/246, 104/246, 9/206, 25/206, 29/206, 21/246)
  sarc_preds <- multiclass_predict(sarc_params$testdata, sarc$sets, sarc_params$zs,
                                   sarc_params$ss, priors, c("Dedi", "Leio", "Malig",
                                                             "Myxo", "Pleo", "Undif"))
  conf <- confusionMatrix(factor(sarc_preds$preds[1,]), factor(sarc_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2])) 
  
}
mean(accs)
mean(sens)
mean(spef)

###############################################################
# UCS
###############################################################

malig <- read.csv("rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/malignant_mixed_mullerian_tumor_rsem_processed.csv",
                  sep = "\t")
mheter <- read.csv("rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/mmmt_heterologous_type_rsem_processed.csv",
                   sep = "\t")
mhom <- read.csv("rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/mmmt_homologous_type_rsem_processed.csv",
                 sep = "\t")


ucs_sets <- list()
ucs_sets[[1]] <- malig
ucs_sets[[2]] <- mheter
ucs_sets[[3]] <- mhom

ucs_data <- cbind(malig, mheter, mhom)
nucs <- list()
nucs$x <- as.matrix(ucs_data)
nucs$y <- ucs$classes[,1]
ucs_nsc <- pamr.train(nucs)
pamr.geneplot(ucs_nsc, nucs, threshold = 2.618)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  ucs <- multiclass_preprocess(ucs_sets, c("Malig", "Mheter", "Mhom"),
                               c(10044, 5690, 9439, 6133, 7747, 11271, 1030, 3483,
                                 1669, 9485, 4991, 10338, 1156, 10415))
  if (k == 1){
    folds <- do_folds(ucs$sets, ucs$data)
  }
  ucs_params <- multiclass_train(ucs$sets, ucs$data, folds)
  priors <- c(24/57, 20/57, 13/57)
  ucs_preds <- multiclass_predict(ucs_params$testdata, ucs$sets, ucs_params$zs,
                                  ucs_params$ss, priors, c("Malig", "Mheter", "Mhom"))
  conf <- confusionMatrix(factor(ucs_preds$preds[1,]), factor(ucs_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2])) 
  
}
mean(accs)
mean(sens)
mean(spef)

###############################################################
# THYMOMA
###############################################################

typea <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_a_rsem_processed.csv",
                  sep = "\t")
typeab <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_ab_rsem_processed.csv",
                   sep = "\t")
typeb1 <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b1_rsem_processed.csv",
                   sep = "\t")
typeb2 <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b2_rsem_processed.csv",
                   sep = "\t")
typeb3 <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_b3_rsem_processed.csv",
                   sep = "\t")
typec <- read.csv("rnaseqV2data/processed_rnaseqV2/thym/thym_rsem_processed/thymoma_type_c_rsem_processed.csv",
                  sep = "\t")

thym_sets <- list()
thym_sets[[1]] <- typea
thym_sets[[2]] <- typeab
thym_sets[[3]] <- typeb1
thym_sets[[4]] <- typeb2
thym_sets[[5]] <- typeb3
thym_sets[[6]] <- typec

thym_data <- cbind(typea, typeab, typeb1, typeb2, typeb3, typec)
nthym <- list()
nthym$x <- as.matrix(thym_data)
nthym$y <- thym$classes[,1]
thym_nsc <- pamr.train(nthym)
pamr.geneplot(thym_nsc, nthym, threshold = 8.734)

accs <- c()
sens <- c()
spef <- c()
for(k in 1:5){
  
  thym <- multiclass_preprocess(thym_sets, c("Typea", "Typeab", "Typeb1",
                                             "Typeb2", "Typeb3", "Typec"),
                                c(8561, 2603, 13567, 12897, 6921, 12034, 12780))
  if (k == 1){
    folds <- do_folds(thym$sets, thym$data)
  }
  thym_params <- multiclass_train(thym$sets, thym$data, folds)
  priors <- c(17/120, 35/120, 15/120, 31/120, 11/120, 11/120)
  thym_preds <- multiclass_predict(thym_params$testdata, thym$sets, thym_params$zs,
                                   thym_params$ss, priors, c("Typea", "Typeab", "Typeb1",
                                                             "Typeb2", "Typeb3", "Typec"))
  conf <- confusionMatrix(factor(thym_preds$preds[1,]), factor(thym_preds$preds[2,]))
  accs <- c(accs, conf$overall[[1]]) 
  sens <- c(sens, mean(conf$byClass[,1]))
  spef <- c(spef, mean(conf$byClass[,2])) 
  
}
mean(accs)
mean(sens)
mean(spef)
