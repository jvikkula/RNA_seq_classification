library(dplyr)
library(sSeq)
library(caret)

setwd("~/Documents/Kurssimateriaalit/LST_project")

# Preprocessing data
load("BRCA_with_clinical.RData")

create_datasets <- function(df, classes){

idx <- which(colnames(df$clinical) == "histologicaltype")
ny <- nrow(df$merged.dat)
samples <- df$merged.dat$bcr

y <- as.integer(factor(df$clinical[samples,idx], levels = classes)) - 1

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
  acc <- accuracy(y_test, pred)
  return(list(acc = acc, pred = pred, target = y_test))
}

# Load the datasets

load("BRCA_with_clinical.RData")
load("LAML_clinical_RNASeq.RData")
load("LUAD_clinical_RNASeq.RData")
load("LUSC_clinical_RNASeq.RData")
load("OV_clinical_RNASeq.RData")
load("THCA_clinical_RNASeq.RData")

# Binary prediction
brca <- create_datasets(rnaseq.brca, c("infiltrating ductal carcinoma", "infiltrating lobular carcinoma"))

# Binary prediction
luad <- create_datasets(rnaseq.luad, c("lung adenocarcinoma- not otherwise specified (nos)", "lung adenocarcinoma mixed subtype"))

# Binary prediction
lusc <- create_datasets(rnaseq.lusc, c("lung squamous cell carcinoma- not otherwise specified (nos)", "lung basaloid squamous cell carcinoma"))

# Three-way classification
thca <- create_datasets(rnaseq.thca, c("thyroid papillary carcinoma - classical/usual","thyroid papillary carcinoma - follicular (>= 99% follicular patterned)","thyroid papillary carcinoma - tall cell (>= 50% tall cell features)"))

brca.pred <- NBLDA_predict(rnaseq.brca, brca, 0.7, "mle")
confusionMatrix(as.factor(brca.pred$y), as.factor(brca.pred$target))

luad.pred <- NBLDA_predict(rnaseq.luad, luad, 0.7, "mle")
lusc.pred <- NBLDA_predict(rnaseq.lusc, lusc, 0.7, "mle")
thca.pred <- NBLDA_predict(rnaseq.thca, thca, 0.7, "mle")


print(brca.pred$acc)
print(luad.pred$acc)
print(lusc.pred$acc)
print(thca.pred$acc)

confusionMatrix(as.factor(brca.pred$y), as.factor(brca.pred$target))
confusionMatrix(luad.pred$y, luad.pred$target)
confusionMatrix(lusc.pred$y, lusc.pred$target)
confusionMatrix(thca.pred$y, thca.pred$target)

y <- brca$y
samples <- brca$samples

y <- na.omit(y)
samples <- na.exclude(samples)
indices <- which(rnaseq.brca$merged.dat$bcr %in% samples) + 4
counts <- rnaseq.brca$merged.dat[indices,4:length(colnames(rnaseq.brca$merged.dat))]
counts <- na.omit(counts)
# #
genes <- sample(1:ncol(counts), 40, replace = F)
counts <- counts[,genes]
train <- sample_frac(counts, 0.7)
sid <- as.numeric(rownames(train))
y_train <- y[sid]
test <- counts[-sid,]
y_test <- y[-sid]

disp <- estimate_dispersion(train)

pred <- NBLDA(train, test, y_train, disp, method = "mle")
confusionMatrix(as.factor(pred), as.factor(y_test))
acc <- accuracy(y_test, pred)
print("Prediction accuracy with mle")
print(acc)
# #
# #
pred2 <- NBLDA_Toy_Example(train, test, y_train)
y2 <- pred2$yte
acc2 <- accuracy(y_test, y2)
print(acc2)