library(dplyr)
library(sSeq)
library(caret)

setwd("~/Documents/Kurssimateriaalit/LST_project")
source("RNA_seq_classification/NBLDA.R")

kfolds <- TRUE
set.seed(123)

# Preprocessing data
# load("BRCA_with_clinical.RData")

create_datasets <- function(df, classes){

idx <- which(colnames(df$clinical) == "histologicaltype")
ny <- nrow(df$merged.dat)
samples <- df$merged.dat$bcr

y <- as.integer(factor(df$clinical[samples,idx], levels = classes)) - 1
del <- which(is.na(y))
y <- y[-del]
samples <- samples[-del]
return(list(y = y, samples = samples))
}

accuracy <- function(target, prediction){
  correct <- 0
  n <- length(target)
  for (i in 1:n){
    if (target[i] == prediction[i]) {correct <- correct + 1}
  }
  acc <- correct / n
  return(acc)
}

NBLDA_predict <- function(df, dataset, fraction, method) {
  samples <- dataset$samples
  y <- dataset$y
  samples <- na.exclude(samples)
  indices <- which(df$merged.dat$bcr %in% samples) + 4
  counts <- df$merged.dat[indices,4:length(colnames(df$merged.dat))]
  counts <- na.omit(counts)
  
  train <- sample_frac(counts, fraction)
  sid <- as.numeric(rownames(train))
  y_train <- y[sid]
  test <- counts[-sid,]
  y_test <- y[-sid]
  
  disp <- estimate_dispersion(train)
  
  pred <- NBLDA(train, test, y_train, disp, method = method)
  #acc <- accuracy(y_test, pred)
  return(list(pred = pred, target = y_test))
}

simulate_dataset <- function(ngenes, nsamples) {
  dat <- rnbinom(n = ngenes*nsamples, size = 1000, prob = 0.9)
  labs <- rbinom(nsamples, size = 1, prob = 0.2)
  df <- as.data.frame(matrix(dat, nrow = nsamples, byrow = T))
  return(list(dat = df, lab = labs))
}

transpose_df <- function(df) {
  n <- rownames(df)
  
  # transpose all but the first column (name)
  df <- as.data.frame(t(df))
  colnames(df) <- n
  return(df)
}

split_data <- function(combined){
  train <- sample_frac(combined, 0.7)
  sid <- as.numeric(rownames(train))
  test <- combined[-sid,]
  # Labels
  y_test <- test[,ncol(test)]
  y_train <- train[,ncol(train)]
  train <- train[,-ncol(train)]
  test <- test[,-ncol(test)]
  return(list(test = test, train = train, y_test = y_test, y_train = y_train))
}

kfold <- function(splitlist, df) {
  accuracies <- numeric()
  conf.matrices <- list()
  for (i in 1:length(splitlist)){
    gc()
    train.idx <- c(unname(unlist(splitlist[-i])))
    test.idx <- unname(unlist(splitlist[i]))
    train <- df[train.idx,]
    y_train <- train[,ncol(train)]
    test <- df[test.idx,]
    y_test <- test[,ncol(test)]
    train <- train[,-ncol(train)]
    test <- test[,-ncol(test)]
    
    disp <- estimate_dispersion(train)
    pred <- NBLDA(train, test, y_train, disp, method = "mle")
    conf <- confusionMatrix(factor(y_test, levels = levels(as.factor(pred$pred))), as.factor(pred$pred))
    conf.matrices[[paste0(i)]] <- conf
    print(conf)
    accuracies[i] <- conf$overall["Accuracy"]
  }
  avg <- mean(accuracies)
  print(avg)
  return(list(accuracies=accuracies, conf.matrices = conf.matrices))
}

###################
## BRCA-analysis ##
###################


brca.duct <- read.csv("data/rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/brca_ductal_rsem_processed.csv", header = TRUE, sep = '\t')
brca.lob <- read.csv("data/rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/brca_lobular_rsem_processed.csv", header = TRUE, sep = '\t')

brca.duct <- rbind(brca.duct, rep(0))
brca.lob <- rbind(brca.lob, rep(1))
combined <- cbind(brca.duct, brca.lob)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)
rm(brca.duct, brca.lob)

if (kfolds) {
  print("Starting BRCA analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/BRCA_kfold.RData")
} else {
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.brca <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}
#orig <- NBLDA_orig(final$train, final$y_train, final$test, disp, type = "mle")
#conf.o.brca <- confusionMatrix(as.factor(final$y_test), as.factor(orig$ytehat))


######################
## COAREAD-analysis ##
######################

col.ade <- read.csv("data/rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed/colon_adenocarcinoma_raw.csv", header = TRUE, sep = '\t')
col.muc <- read.csv("data/rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed/colon_mucinous_adenocarcinoma_raw.csv", header = TRUE, sep = '\t')
rec.ade <- read.csv("data/rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed/rectal_adenocarcinoma_raw.csv", header = TRUE, sep = '\t') 
rec.muc <- read.csv("data/rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed/rectal_mucinous_adenocarcinoma_raw.csv", header = TRUE, sep = '\t')

col.ade <- rbind(col.ade, rep(0))
col.muc <- rbind(col.muc, rep(1))
rec.ade <- rbind(rec.ade, rep(2))
rec.muc <- rbind(rec.muc, rep(3))
combined <- cbind(col.ade, col.muc, rec.ade, rec.muc)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)

rm(col.ade, col.muc, rec.ade, rec.muc)

if (kfolds) {
  print("Starting COAREAD analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/COAREAD_kfold.RData")
} else {
  combined <- combined[sample(nrow(combined)),]
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.coaread <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}

#orig <- NBLDA_orig(final$train, final$y_train, final$test, disp, type = "mle")
#conf.o.brca <- confusionMatrix(as.factor(final$y_test), as.factor(orig$pred))

####################
## KIPAN-analysis ##
####################

kipan.chromo <- read.csv("data/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_chromophobe_rsem_processed.csv", header = TRUE, sep = '\t')
kipan.clear <- read.csv("data/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_clear_cell_renal_carcinoma_rsem_processed.csv", header = TRUE, sep = '\t')
kipan.pap <- read.csv("data/rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/kidney_papillary_renal_cell_carcinoma_rsem_processed.csv", header = TRUE, sep = '\t')

kipan.chromo <- rbind(kipan.chromo, rep(0))
kipan.clear <- rbind(kipan.clear, rep(1))
kipan.pap <- rbind(kipan.pap, rep(2))
combined <- cbind(kipan.chromo, kipan.clear, kipan.pap)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)

rm(kipan.chromo, kipan.clear, kipan.pap)

if (kfolds) {
  print("Starting KIPAN analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/KIPAN_kfold.RData")
} else {
  combined <- combined[sample(nrow(combined)),]
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.kipan <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}

##################
## LGG-analysis ##
##################

lgg.astro <- read.csv("data/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/astrocytoma_rsem_processed.csv", header = TRUE, sep = '\t')
lgg.oligo <- read.csv("data/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/oligoastrocytoma_rsem_processed.csv", header = TRUE, sep = '\t')
lgg.oligden <- read.csv("data/rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/oligodendroglioma_rsem_processed.csv", header = TRUE, sep = '\t')

lgg.astro <- rbind(lgg.astro, rep(0))
lgg.oligo <- rbind(lgg.oligo, rep(1))
lgg.oligden <- rbind(lgg.oligden, rep(2))
combined <- cbind(lgg.astro, lgg.oligo, lgg.oligden)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)

rm(lgg.astro, lgg.oligo, lgg.oligden)

if (kfolds) {
  print("Starting LGG analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/LGG_kfold.RData")
} else {
  combined <- combined[sample(nrow(combined)),]
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.lgg <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}

###################
## MESO-analysis ##
###################

meso.bipha <- read.csv("data/rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/biphasic_mesothelioma_rsem_processed.csv", header = TRUE, sep = '\t')
meso.diff <- read.csv("data/rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/diffuse_malignant_mesothelioma_rsem_processed.csv", header = TRUE, sep = '\t')
meso.epit <- read.csv("data/rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/epithelioid_mesothelioma_rsem_processed.csv", header = TRUE, sep = '\t')
#meso.sarco <- read.csv("data/rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/sarcomatoid_mesotheliomarsem_processed.csv", header = TRUE, sep = '\t')

meso.bipha <- rbind(meso.bipha, rep(0))
meso.diff <- rbind(meso.diff, rep(1))
meso.epit <- rbind(meso.epit, rep(2))
#meso.sarco <- rbind(meso.sarco, rep(3))
combined <- cbind(meso.bipha, meso.diff, meso.epit)
#combined <- cbind(meso.bipha, meso.diff, meso.epit, meso.sarco)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)

rm(meso.bipha, meso.diff, meso.epit)

if (kfolds) {
  print("Starting MESO analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/MESO_kfold.RData")
} else {
  combined <- combined[sample(nrow(combined)),]
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.meso <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}

######################
## SARCOMA-analysis ##
######################

sarcoma.ded <- read.csv("data/rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/dedifferentiated_liposarcoma_raw.csv", header = TRUE, sep = '\t')
sarcoma.leio <- read.csv("data/rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/leiomyosarcoma_raw.csv", header = TRUE, sep = '\t')
sarcoma.mal <- read.csv("data/rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/malignant_peripheral_nerve_sheath_tumors_raw.csv", header = TRUE, sep = '\t')
sarcoma.myxo <- read.csv("data/rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/myxofibrosarcoma_raw.csv", header = TRUE, sep = '\t')
sarcoma.pleo <- read.csv("data/rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/pleomorphic_mfh_raw.csv", header = TRUE, sep = '\t')
sarcoma.undif <- read.csv("data/rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/undifferentiated_pleomorphic_sarcoma_raw.csv", header = TRUE, sep = '\t')

sarcoma.ded <- rbind(sarcoma.ded, rep(0))
sarcoma.leio <- rbind(sarcoma.leio, rep(1))
sarcoma.mal <- rbind(sarcoma.mal, rep(2))
sarcoma.myxo <- rbind(sarcoma.myxo, rep(3))
sarcoma.pleo <- rbind(sarcoma.pleo, rep(4))
sarcoma.undif <- rbind(sarcoma.undif, rep(5))
combined <- cbind(sarcoma.ded, sarcoma.leio, sarcoma.mal, sarcoma.myxo, sarcoma.pleo, sarcoma.undif)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)

rm(sarcoma.ded, sarcoma.leio, sarcoma.mal, sarcoma.myxo, sarcoma.pleo, sarcoma.undif)

if (kfolds) {
  print("Starting SARCOMA analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/SARCOMA_kfold.RData")
} else {
  combined <- combined[sample(nrow(combined)),]
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.sarcoma <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}

###################
## THCA-analysis ##
###################

thca.clas <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_classical_rsem_processed.csv", header = TRUE, sep = '\t')
thca.fol <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_follicular_rsem_processed.csv", header = TRUE, sep = '\t')
thca.tall <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/thyroid_papillary_carcinoma_tall_cell_rsem_processed.csv", header = TRUE, sep = '\t')

thca.clas <- rbind(thca.clas, rep(0))
thca.fol <- rbind(thca.fol, rep(1))
thca.tall <- rbind(thca.tall, rep(2))
combined <- cbind(thca.clas, thca.fol, thca.tall)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)

rm(thca.clas, thca.fol, thca.tall)

if (kfolds) {
  print("Starting THCA analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/THCA_kfold.RData")
} else {
  combined <- combined[sample(nrow(combined)),]
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.thca <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}


######################
## THYMOMA-analysis ##
######################

thymoma.ab <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/thymoma_type_ab_rsem_processed.csv", header = TRUE, sep = '\t')
thymoma.a <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/thymoma_type_a_rsem_processed.csv", header = TRUE, sep = '\t')
thymoma.b1 <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/thymoma_type_b1_rsem_processed.csv", header = TRUE, sep = '\t')
thymoma.b2 <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/thymoma_type_b2_rsem_processed.csv", header = TRUE, sep = '\t')
thymoma.b3 <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/thymoma_type_b3_rsem_processed.csv", header = TRUE, sep = '\t')
thymoma.c <- read.csv("data/rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/thymoma_type_c_raw.csv", header = TRUE, sep = '\t')

thymoma.ab <- rbind(thymoma.ab, rep(0))
thymoma.a <- rbind(thymoma.a, rep(1))
thymoma.b1 <- rbind(thymoma.b1, rep(2))
thymoma.b2 <- rbind(thymoma.b2, rep(3))
thymoma.b3 <- rbind(thymoma.b3, rep(4))
thymoma.c <- rbind(thymoma.c, rep(5))
combined <- cbind(thymoma.ab, thymoma.a, thymoma.b1, thymoma.b2, thymoma.b3, thymoma.c)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)

rm(thymoma.ab, thymoma.a, thymoma.b1, thymoma.b2, thymoma.b3, thymoma.c)

if (kfolds) {
  print("Starting THYMOMA analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/THYMOMA_kfold.RData")
} else {
  combined <- combined[sample(nrow(combined)),]
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.thymoma <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}

###################
## UCEC-analysis ##
###################

ucec.endo <- read.csv("data/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/endometrioid_endometrial_adenocarcinoma_rsem_processed.csv", header = TRUE, sep = '\t')
ucec.mixed <- read.csv("data/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/mixed_serous_and_endometrioid_rsem_processed.csv", header = TRUE, sep = '\t')
ucec.serous <- read.csv("data/rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/serous_endometrial_adenocarcinoma_rsem_processed.csv", header = TRUE, sep = '\t')

ucec.endo <- rbind(ucec.endo, rep(0))
ucec.mixed <- rbind(ucec.mixed, rep(1))
ucec.serous <- rbind(ucec.serous, rep(2))
combined <- cbind(ucec.endo, ucec.mixed, ucec.serous)
split <- createFolds(combined, k = 5)
combined <- transpose_df(combined)

rm(ucec.endo, ucec.mixed, ucec.serous)

if (kfolds) {
  print("Starting UCEC analysis")
  start <- Sys.time()
  res <- kfold(split, combined)
  end <- Sys.time()
  save(res, file = "results/UCEC_kfold.RData")
} else {
  combined <- combined[sample(nrow(combined)),]
  final <- split_data(combined)
  disp <- estimate_dispersion(final$train)
  pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
  conf.ucec <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
}

load("results/BRCA_kfold.RData")
mean(res$accuracies)

load("results/COAREAD_kfold.RData")
mean(res$accuracies)

load("results/KIPAN_kfold.RData")
mean(res$accuracies)

load("results/LGG_kfold.RData")
mean(res$accuracies)

load("results/MESO_kfold.RData")
mean(res$accuracies)

load("results/SARCOMA_kfold.RData")
mean(res$accuracies)

load("results/THCA_kfold.RData")
mean(res$accuracies)

load("results/THYMOMA_kfold.RData")
mean(res$accuracies)

load("results/UCEC_kfold.RData")
mean(res$accuracies)

##################
## UCS-analysis ##
##################

#ucs.mal <- read.csv("data/rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/malignant_mixed_mullerian_tumor_rsem_processed.csv", header = TRUE, sep = '\t')
# ucs.mmmt.hetero <- read.csv("data/rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/mmmt_heterologous_type_rsem_processed.csv", header = TRUE, sep = '\t')
# ucs.mmmt.homo <- read.csv("data/rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/mmmt_homologous_type_rsem_processed.csv", header = TRUE, sep = '\t')
# 
# ucs.mal <- rbind(ucs.mal, rep(0))
# ucs.mmmt.hetero <- rbind(ucs.mmmt.hetero, rep(1))
# ucs.mmmt.homo <- rbind(ucs.mmmt.homo, rep(2))
# combined <- cbind(ucs.mal, ucs.mmmt.hetero, ucs.mmmt.homo)
# split <- createFolds(combined, k = 5)
# combined <- transpose_df(combined)
# 
# rm(ucs.mal, ucs.mmmt.hetero, ucs.mmmt.homo)
# 
# if (kfolds) {
#   start <- Sys.time()
#   res <- kfold(split, combined)
#   end <- Sys.time()
#   save(res, file = "results/USC_kfold.RData")
# } else {
#   combined <- combined[sample(nrow(combined)),]
#   final <- split_data(combined)
#   disp <- estimate_dispersion(final$train)
#   pred <- NBLDA(final$train, final$test, final$y_train, disp, method = "mle")
#   conf.ucs <- confusionMatrix(as.factor(final$y_test), as.factor(pred$pred))
# }

if (kfolds = FALSE) {
confusion.matrices <- list("brca" = conf.brca, "coadread" = conf.coaread, "kipan" = conf.kipan, "lgg" = conf.lgg,
                           "meso" = conf.meso, "sarcoma" = conf.sarcoma, "thca" = conf.thca, "thymoma" = conf.thymoma,
                           "ucec" = conf.ucec)
save(confusion.matrices, file = "confusion_matrices_15_04.RData")
}

# Load the datasets
# 
# load("BRCA_with_clinical.RData")
# load("LAML_clinical_RNASeq.RData")
# load("LUAD_clinical_RNASeq.RData")
# load("LUSC_clinical_RNASeq.RData")
# load("OV_clinical_RNASeq.RData")
# load("THCA_clinical_RNASeq.RData")
# 
# # Binary prediction
# brca <- create_datasets(rnaseq.brca, c("infiltrating ductal carcinoma", "infiltrating lobular carcinoma"))
# 
# # Binary prediction
# luad <- create_datasets(rnaseq.luad, c("lung adenocarcinoma- not otherwise specified (nos)", "lung adenocarcinoma mixed subtype"))
# 
# # Binary prediction
# lusc <- create_datasets(rnaseq.lusc, c("lung squamous cell carcinoma- not otherwise specified (nos)", "lung basaloid squamous cell carcinoma"))
# 
# # Three-way classification
# thca <- create_datasets(rnaseq.thca, c("thyroid papillary carcinoma - classical/usual","thyroid papillary carcinoma - follicular (>= 99% follicular patterned)","thyroid papillary carcinoma - tall cell (>= 50% tall cell features)"))
# 
# # brca.pred <- NBLDA_predict(rnaseq.brca, brca, 0.7, "mle")
# # confusionMatrix(as.factor(brca.pred$y), as.factor(brca.pred$target))
# # 
# # luad.pred <- NBLDA_predict(rnaseq.luad, luad, 0.7, "mle")
# # lusc.pred <- NBLDA_predict(rnaseq.lusc, lusc, 0.7, "mle")
# # thca.pred <- NBLDA_predict(rnaseq.thca, thca, 0.7, "mle")
# # 
# # 
# # print(brca.pred$acc)
# # print(luad.pred$acc)
# # print(lusc.pred$acc)
# # print(thca.pred$acc)
# # 
# # confusionMatrix(as.factor(brca.pred$y), as.factor(brca.pred$target))
# # confusionMatrix(luad.pred$y, luad.pred$target)
# # confusionMatrix(lusc.pred$y, lusc.pred$target)
# # confusionMatrix(thca.pred$y, thca.pred$target)
# 
# df <- simulate_dataset(5, 20)
# print(df)
# train <- sample_frac(df$dat, 0.5)
# sid <- as.numeric(rownames(train))
# y_train <- df$lab[sid]
# test <- df$dat[-sid,]
# y_test <- df$lab[-sid]
# 
# disp <- estimate_dispersion(train)
# mine <- NBLDA(train, test, y_train, disp, method = "quantile")
# acc <- accuracy(y_test, pred)
# 
# orig <- NBLDA_orig(train, y_train, test, disp, type = "quantile")
# print(orig$ytehat)
# print(mine$pred)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# y <- brca$y
# samples <- brca$samples
# 
# #y <- na.omit(y)
# #samples <- na.exclude(samples)
# #indices <- which(rnaseq.brca$merged.dat$bcr %in% samples) + 4
# counts <- rnaseq.brca$merged.dat[,4:length(colnames(rnaseq.brca$merged.dat))]
# counts <- na.omit(counts)
# #y <- y[indices-4]
# 
# #genes <- sample(1:ncol(counts), 40, replace = F)
# #counts <- counts[,genes]
# train <- sample_frac(counts, 0.7)
# sid <- as.numeric(rownames(train))
# y_train <- y[sid]
# test <- counts[-sid,]
# y_test <- y[-sid]
# 
# disp <- estimate_dispersion(train)
# 
# pred <- NBLDA(train, test, y_train, disp, method = "mle")
# confusionMatrix(factor(pred$pred, levels = c(0,1)), as.factor(y_test))
# acc <- accuracy(pred$pred, y_test)
# print("Prediction accuracy with mle")
# print(acc)
# # #
# # #
# pred2 <- NBLDA_Toy_Example(train, test, y_train)
# y2 <- pred2$yte
# acc2 <- accuracy(y_test, y2)
# print(acc2)
