# Use this file to analyse performance of PLDA and ZIPLDA functions
require(caret)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('PLDA.R')
source('ziplda.R')

# Define k-fold cross validation functions for PLDA and ZIPLDA to decrease redundancy
get.accuracy <- function(data, y, func = c('PLDA','ZIPLDA'), prior, method){
  set.seed(123)
  flds <- createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)
  acc <- c()
  output <- list()
  start <- Sys.time()
  for (k in 1:5){
    print(k)
    test <- t(data[,flds[[k]]])
    train <- t(data[,-flds[[k]]])
    ytrue <- y[flds[[k]]]
    ytrain <- y[-flds[[k]]]
    if (func == 'PLDA'){
      out <- PLDA(train, test, ytrain, prior = prior, method = method)
    }
    else if (func == 'ZIPLDA'){
      prob0 <- point.mass.at.zero(train, test, ytrain, beta=1, method=method)
      out <- ziplda(train, test, ytrain, prior = prior, method = method, prob0 = prob0)
    }
    #output[[paste('yhat',k,sep='')]] <- out$yhat
    #output[[paste('ytrue',k,sep='')]] <- ytrue
    tb <- table(ytrue, out$yhat)
    acc[k] <- sum(diag(tb))/sum(tb)
  }
  stop <- Sys.time()
  print(stop-start)
  cat('Accuracies are:', acc, '\n')
  cat('Mean prediction accuracy for', method, 'with 5-fold cross validation is', round(mean(acc)*100,2), '% \n')
  #return(output)
}


#########################################
########## RSEM PROCESSED DATA ##########
#########################################

#### BRCA: ductal and lobular / all labels ####
# Read and combine BRCA data
path <- '../rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/'
ductal <- read.csv(paste(path, 'brca_ductal_rsem_processed.csv', sep=''), sep='')
lobular <- read.csv(paste(path, 'brca_lobular_rsem_processed.csv', sep=''), sep='')
data <- cbind(ductal, lobular)
# Define labels
y <- c(rep(1,dim(ductal)[2]), rep(2,dim(lobular)[2]))
# Define prior
prior <- c(dim(ductal)[2], dim(lobular)[2])/length(y)
# Compute prediction accuracy
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)
rm(ductal)
rm(lobular)

#### COADREAD: all labels ####
# colon and rectal carcinomas are very similar: https://www.nature.com/articles/nature11252
# maybe this is why the classifictaion task is so difficult
path <- '../rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed/'
colon <- read.csv(paste(path, 'colon_adenocarcinoma_raw.csv', sep=''), sep='')
colon_muc <- read.csv(paste(path, 'colon_mucinous_adenocarcinoma_raw.csv',sep=''),sep='')
rectal <- read.csv(paste(path, 'rectal_adenocarcinoma_raw.csv',sep=''),sep='')
data <- cbind(colon, colon_muc, rectal)
# Define labels
y <- c(rep(1, dim(colon)[2]), rep(2, dim(colon_muc)[2]), rep(3, dim(rectal)[2]))
# Deifne prior
prior <- c(dim(colon)[2], dim(colon_muc)[2], dim(rectal)[2])/length(y)
#Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variable to save memory
rm(data)
rm(colon)
rm(colon_muc)
rm(rectal)

#### KIPAN: clear cell and papillary renal carcinoma ####
# kidney
path <- '../rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/'
clear <- read.csv(paste(path, 'kidney_clear_cell_renal_carcinoma_rsem_processed.csv', sep=''), sep='')
papillary <- read.csv(paste(path, 'kidney_papillary_renal_cell_carcinomarsem_processed_log2.csv', sep=''),sep='')
data <- cbind(clear, papillary)
# Define labels
y <- c(rep(1, dim(clear)[2]), rep(2, dim(papillary)[2]))
# Define prior
prior <- c(dim(clear)[2], dim(papillary)[2])/length(y)
#Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
putput <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)

#### KIPAN: all labels ####
chromophobe <- read.csv(paste(path, 'kidney_chromophobe_rsem_processed.csv', sep=''), sep='')
data <- cbind(chromophobe, clear, papillary)
# Define labels
y <- c(rep(1, dim(chromophobe)[2]), rep(2, dim(clear)[2]), rep(3, dim(papillary)[2]))
# Define prior
prior <- c(dim(chromophobe)[2], dim(clear)[2], dim(papillary)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(clear)
rm(papillary)
rm(chromophobe)
rm(data)

#### LGG: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/'
astro <- read.csv(paste(path, 'astrocytoma_rsem_processed_log2.csv', sep=''), sep='')
oligo_dendroglioma <- read.csv(paste(path, 'oligodendroglioma_rsem_processed.csv', sep=''), sep='')
oligo_astro <- read.csv(paste(path, 'oligoastrocytoma_rsem_processed.csv', sep=''), sep='')
data <- cbind(astro, oligo_astro, oligo_dendroglioma)
# Define labels
y <- c(rep(1, dim(astro)[2]), rep(2, dim(oligo_astro)[2]), rep(3, dim(oligo_dendroglioma)[2]))
# Define prior
prior <- c(dim(astro)[2], dim(oligo_astro)[2], dim(oligo_dendroglioma)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(astro)
rm(oligo_dendroglioma)
rm(oligo_astro)
rm(data)

#### MESO: epithelioid and biphasic malignant / all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/'
epithelioid <- read.csv(paste(path, 'epithelioid_mesothelioma_rsem_processed.csv', sep=''), sep = '')
biphasic <- read.csv(paste(path, 'biphasic_mesothelioma_rsem_processed.csv', sep=''), sep='')
data <- cbind(biphasic, epithelioid)
# Define labels
y <- c(rep(1, dim(biphasic)[2]), rep(2, dim(epithelioid)[2]))
# Define prior
prior <- c(dim(biphasic)[2], dim(epithelioid)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)
rm(epithelioid)
rm(biphasic)

#### SARCOMA: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/'
leiomyosarcoma <- read.csv(paste(path, 'leiomyosarcoma_raw.csv', sep=''), sep='')
liposarcoma <- read.csv(paste(path, 'dedifferentiated_liposarcoma_raw.csv', sep=''), sep='')
pleomorphic <- read.csv(paste(path, 'pleomorphic_mfh_raw.csv', sep=''), sep='')
myxofibrosarcoma <- read.csv(paste(path, 'myxofibrosarcoma_raw.csv', sep=''), sep='')
undiff_pleomorphic <- read.csv(paste(path, 'undifferentiated_pleomorphic_sarcoma_raw.csv', sep=''), sep='')
nerve <- read.csv(paste(path, 'malignant_peripheral_nerve_sheath_tumors_raw.csv', sep=''), sep='')
#data <- cbind(leiomyosarcoma, liposarcoma, pleomorphic, myxofibrosarcoma, undiff_pleomorphic, nerve)
data <- cbind(liposarcoma, leiomyosarcoma, nerve, myxofibrosarcoma, pleomorphic, undiff_pleomorphic)
# Define labels
y <- c(rep(1, dim(liposarcoma)[2]), rep(2, dim(leiomyosarcoma)[2]), rep(3, dim(nerve)[2]), rep(4, dim(myxofibrosarcoma)[2]), 
       rep(5, dim(pleomorphic)[2]), rep(6, dim(undiff_pleomorphic)[2]))
# Define prior
prior <- c(dim(liposarcoma)[2], dim(leiomyosarcoma)[2], dim(nerve)[2], dim(myxofibrosarcoma)[2], 
           dim(pleomorphic)[2], dim(undiff_pleomorphic)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(leiomyosarcoma)
rm(liposarcoma)
rm(pleomorphic)
rm(myxofibrosarcoma)
rm(undiff_pleomorphic)
rm(nerve)
rm(data)

#### THCA: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/'
classical <- read.csv(paste(path, 'thyroid_papillary_carcinoma_classical_rsem_processed.csv', sep=''), sep='')
follicular <- read.csv(paste(path, 'thyroid_papillary_carcinoma_follicular_rsem_processed.csv', sep=''), sep='')
tall_cell <- read.csv(paste(path, 'thyroid_papillary_carcinoma_tall_cell_rsem_processed.csv', sep=''), sep='')
data <- cbind(classical, follicular, tall_cell)
# Define labels
y <- c(rep(1, dim(classical)[2]), rep(2, dim(follicular)[2]), rep(3, dim(tall_cell)[2]))
# Define prior
prior <- c(dim(classical)[2], dim(follicular)[2], dim(tall_cell)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(classical)
rm(follicular)
rm(tall_cell)
rm(data)

#### THYMOMA: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/'
ab <- read.csv(paste(path, 'thymoma_type_ab_rsem_processed.csv', sep=''), sep='')
b2 <- read.csv(paste(path, 'thymoma_type_b2_rsem_processed_log2.csv', sep=''), sep='')
a <- read.csv(paste(path, 'thymoma_type_a_rsem_processed.csv', sep=''), sep='')
b1 <- read.csv(paste(path, 'thymoma_type_b1_rsem_processed.csv', sep=''), sep='')
b3 <- read.csv(paste(path, 'thymoma_type_b3_rsem_processed_log2.csv', sep=''), sep='')
c <- read.csv(paste(path, 'thymoma_type_c_raw.csv', sep=''), sep='')
data <- cbind(a, ab, b1, b2, b3, c)
# Define labels
y <- c(rep(1, dim(a)[2]), rep(2, dim(ab)[2]), rep(3, dim(b1)[2]), rep(4, dim(b2)[2]), rep(5, dim(b3)[2]), rep(6, dim(c)[2]))
# Define prior
prior <- c(dim(a)[2], dim(ab)[2], dim(b1)[2], dim(b2)[2], dim(b3)[2], dim(c)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(a)
rm(ab)
rm(b1)
rm(b2)
rm(b3)
rm(c)
rm(data)

#### UCEC: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/'
endometrioid <- read.csv(paste(path, 'endometrioid_endometrial_adenocarcinoma_rsem_processed.csv', sep=''),sep='')
serous <- read.csv(paste(path, 'serous_endometrial_adenocarcinoma_rsem_processed.csv', sep=''), sep='')
mixed <- read.csv(paste(path, 'mixed_serous_and_endometrioid_rsem_processed.csv', sep=''), sep='')
data <- cbind(endometrioid, mixed, serous)
# Define labels
y <- c(rep(1, dim(endometrioid)[2]), rep(2, dim(mixed)[2]), rep(3, dim(serous)[2]))
# Define prior
prior <- c(dim(endometrioid)[2], dim(mixed)[2], dim(serous)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(endometrioid)
rm(serous)
rm(mixed)
rm(data)

#### UCS: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/'
mmmt <- read.csv(paste(path, 'malignant_mixed_mullerian_tumor_rsem_processed.csv', sep=''), sep='')
mmmt <- mmmt[1:14975,]
mmmt_heterologous <- read.csv(paste(path, 'mmmt_heterologous_type_rsem_processed.csv', sep=''), sep='')
mmmt_homologous <- read.csv(paste(path, 'mmmt_homologous_type_rsem_processed.csv', sep=''), sep='')
data <- cbind(mmmt, mmmt_heterologous, mmmt_homologous)
# Define labels
y <- c(rep(1, dim(mmmt)[2]), rep(2, dim(mmmt_heterologous)[2]), rep(3, dim(mmmt_homologous)[2]))
# Define prior
prior <- c(dim(mmmt)[2], dim(mmmt_heterologous)[2], dim(mmmt_homologous)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(mmmt)
rm(mmmt_heterologous)
rm(mmmt_homologous)
rm(data)


##############################################
########## RSEM PROCESSED LOG2 DATA ##########
##############################################

#### BRCA: ductal and lobular / all labels ####
# Read and combine BRCA data
path <- '../rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed_log2/'
ductal <- read.csv(paste(path, 'brca_ductal_rsem_processed_log2.csv', sep=''), sep='')
lobular <- read.csv(paste(path, 'brca_lobular_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(ductal, lobular)
# Define labels
y <- c(rep(1,dim(ductal)[2]), rep(2,dim(lobular)[2]))
# Define prior
prior <- c(dim(ductal)[2], dim(lobular)[2])/length(y)
# Compute prediction accuracy
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)
rm(ductal)
rm(lobular)

#### COAREAD: all labels ####
# colon and rectal carcinomas are very similar: https://www.nature.com/articles/nature11252
# maybe this is why the classifictaion task is so difficult
path <- '../rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed_log2/'
colon <- read.csv(paste(path, 'colon_adenocarcinoma_rsem_processed_log2.csv', sep=''), sep='')
colon_muc <- read.csv(paste(path, 'colon_mucinous_adenocarcinoma_rsem_processed_log2.csv',sep=''),sep='')
rectal <- read.csv(paste(path, 'rectal_adenocarcinoma_log2.csv',sep=''),sep='')
data <- cbind(colon, colon_muc, rectal)
# Define labels
y <- c(rep(1, dim(colon)[2]), rep(2, dim(colon_muc)[2]), rep(3, dim(rectal)[2]))
# Deifne prior
prior <- c(dim(colon)[2], dim(colon_muc)[2], dim(rectal)[2])/length(y)
#Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variable to save memory
rm(data)
rm(colon)
rm(colon_muc)
rm(rectal)

#### KIPAN: clear cell and papillary renal carcinoma ####
# kidney
path <- '../rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed_log2/'
clear <- read.csv(paste(path, 'kidney_clear_cell_renal_carcinoma_rsem_processed_log2.csv', sep=''), sep='')
papillary <- read.csv(paste(path, 'kidney_papillary_renal_cell_carcinoma_rsem_processed_log2.csv', sep=''),sep='')
data <- cbind(clear, papillary)
# Define labels
y <- c(rep(1, dim(clear)[2]), rep(2, dim(papillary)[2]))
# Define prior
prior <- c(dim(clear)[2], dim(papillary)[2])/length(y)
#Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
putput <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)

#### KIPAN: all labels ####
chromophobe <- read.csv(paste(path, 'kidney_chromophobe_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(chromophobe, clear, papillary)
# Define labels
y <- c(rep(1, dim(chromophobe)[2]), rep(2, dim(clear)[2]), rep(3, dim(papillary)[2]))
# Define prior
prior <- c(dim(chromophobe)[2], dim(clear)[2], dim(papillary)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(clear)
rm(papillary)
rm(chromophobe)
rm(data)


#### LGG: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed_log2/'
astro <- read.csv(paste(path, 'astrocytoma_rsem_processed_log2.csv', sep=''), sep='')
oligo_dendroglioma <- read.csv(paste(path, 'oligodendroglioma_rsem_processed_log2.csv', sep=''), sep='')
oligo_astro <- read.csv(paste(path, 'oligoastrocytoma_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(astro, oligo_astro, oligo_dendroglioma)
# Define labels
y <- c(rep(1, dim(astro)[2]), rep(2, dim(oligo_astro)[2]), rep(3, dim(oligo_dendroglioma)[2]))
# Define prior
prior <- c(dim(astro)[2], dim(oligo_astro)[2], dim(oligo_dendroglioma)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(astro)
rm(oligo_dendroglioma)
rm(oligo_astro)
rm(data)


#### MESO: epithelioid and biphasic malignant / all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed_log2/'
epithelioid <- read.csv(paste(path, 'epithelioid_mesothelioma_rsem_processed_log2.csv', sep=''), sep = '')
biphasic <- read.csv(paste(path, 'biphasic_mesothelioma_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(biphasic, epithelioid)
# Define labels
y <- c(rep(1, dim(biphasic)[2]), rep(2, dim(epithelioid)[2]))
# Define prior
prior <- c(dim(biphasic)[2], dim(epithelioid)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)
rm(epithelioid)
rm(biphasic)


#### SARCOMA: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed_log2/'
leiomyosarcoma <- read.csv(paste(path, 'leiomyosarcoma_rsem_processed_log2.csv', sep=''), sep='')
liposarcoma <- read.csv(paste(path, 'dedifferentiated_liposarcoma_rsem_processed_log2.csv', sep=''), sep='')
pleomorphic <- read.csv(paste(path, 'pleomorphic_mfh_rsem_processed_log2.csv', sep=''), sep='')
myxofibrosarcoma <- read.csv(paste(path, 'myxofibrosarcoma_rsem_processed_log2.csv', sep=''), sep='')
undiff_pleomorphic <- read.csv(paste(path, 'undifferentiated_pleomorphic_sarcoma_rsem_processed_log2.csv', sep=''), sep='')
nerve <- read.csv(paste(path, 'malignant_peripheral_nerve_sheath_tumors_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(liposarcoma, leiomyosarcoma, nerve, myxofibrosarcoma, pleomorphic, undiff_pleomorphic)
# Define labels
y <- c(rep(1, dim(liposarcoma)[2]), rep(2, dim(leiomyosarcoma)[2]), rep(3, dim(nerve)[2]), rep(4, dim(myxofibrosarcoma)[2]), 
       rep(5, dim(pleomorphic)[2]), rep(6, dim(undiff_pleomorphic)[2]))
# Define prior
prior <- c(dim(liposarcoma)[2], dim(leiomyosarcoma)[2], dim(nerve)[2], dim(myxofibrosarcoma)[2], 
           dim(pleomorphic)[2], dim(undiff_pleomorphic)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(leiomyosarcoma)
rm(liposarcoma)
rm(pleomorphic)
rm(myxofibrosarcoma)
rm(undiff_pleomorphic)
rm(nerve)
rm(data)

#### THCA: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed_log2/'
classical <- read.csv(paste(path, 'thyroid_papillary_carcinoma_classical_rsem_processed_log2.csv', sep=''), sep='')
follicular <- read.csv(paste(path, 'thyroid_papillary_carcinoma_follicular_rsem_processed_log2.csv', sep=''), sep='')
tall_cell <- read.csv(paste(path, 'thyroid_papillary_carcinoma_tall_cell_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(classical, follicular, tall_cell)
# Define labels
y <- c(rep(1, dim(classical)[2]), rep(2, dim(follicular)[2]), rep(3, dim(tall_cell)[2]))
# Define prior
prior <- c(dim(classical)[2], dim(follicular)[2], dim(tall_cell)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(classical)
rm(follicular)
rm(tall_cell)
rm(data)

#### THYMOMA: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed_log2/'
ab <- read.csv(paste(path, 'thymoma_type_ab_rsem_processed_log2.csv', sep=''), sep='')
b2 <- read.csv(paste(path, 'thymoma_type_b2_rsem_processed_log2.csv', sep=''), sep='')
a <- read.csv(paste(path, 'thymoma_type_a_rsem_processed_log2.csv', sep=''), sep='')
b1 <- read.csv(paste(path, 'thymoma_type_b1_rsem_processed_log2.csv', sep=''), sep='')
b3 <- read.csv(paste(path, 'thymoma_type_b3_rsem_processed_log2.csv', sep=''), sep='')
c <- read.csv(paste(path, 'thymoma_type_c_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(a, ab, b1, b2, b3, c)
# Define labels
y <- c(rep(1, dim(a)[2]), rep(2, dim(ab)[2]), rep(3, dim(b1)[2]), rep(4, dim(b2)[2]), rep(5, dim(b3)[2]), rep(6, dim(c)[2]))
# Define prior
prior <- c(dim(a)[2], dim(ab)[2], dim(b1)[2], dim(b2)[2], dim(b3)[2], dim(c)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(a)
rm(ab)
rm(b1)
rm(b2)
rm(b3)
rm(c)
rm(data)

#### UCEC: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed_log2/'
endometrioid <- read.csv(paste(path, 'endometrioid_endometrial_adenocarcinoma_rsem_processed_log2.csv', sep=''),sep='')
serous <- read.csv(paste(path, 'serous_endometrial_adenocarcinoma_rsem_processed_log2.csv', sep=''), sep='')
mixed <- read.csv(paste(path, 'mixed_serous_and_endometrioid_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(endometrioid, mixed, serous)
# Define labels
y <- c(rep(1, dim(endometrioid)[2]), rep(2, dim(mixed)[2]), rep(3, dim(serous)[2]))
# Define prior
prior <- c(dim(endometrioid)[2], dim(mixed)[2], dim(serous)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(endometrioid)
rm(serous)
rm(mixed)
rm(data)

#### UCS: all labels ####
path <- '../rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed_log2/'
mmmt <- read.csv(paste(path, 'malignant_mixed_mullerian_tumor_rsem_processed_log2.csv', sep=''), sep='')
mmmt <- mmmt[1:14975,]
mmmt_heterologous <- read.csv(paste(path, 'mmmt_heterologous_type_rsem_processed_log2.csv', sep=''), sep='')
mmmt_homologous <- read.csv(paste(path, 'mmmt_homologous_type_rsem_processed_log2.csv', sep=''), sep='')
data <- cbind(mmmt, mmmt_heterologous, mmmt_homologous)
# Define labels
y <- c(rep(1, dim(mmmt)[2]), rep(2, dim(mmmt_heterologous)[2]), rep(3, dim(mmmt_homologous)[2]))
# Define prior
prior <- c(dim(mmmt)[2], dim(mmmt_heterologous)[2], dim(mmmt_homologous)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(mmmt)
rm(mmmt_heterologous)
rm(mmmt_homologous)
rm(data)


##############################################
########## RAW DATA ##########################
##############################################

#### BRCA: ductal and lobular / all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/brca/brca_rsem_unprocessed/'
ductal <- read.csv(paste(path, 'brca_ductal_rsem_unprocessed.csv', sep=''), sep='')
lobular <- read.csv(paste(path, 'brca_lobular_rsem_unprocessed.csv', sep=''), sep='')
data <- cbind(ductal, lobular)
# Define labels
y <- c(rep(1,dim(ductal)[2]), rep(2,dim(lobular)[2]))
# Define prior
prior <- c(dim(ductal)[2], dim(lobular)[2])/length(y)
# Compute prediction accuracy
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)
rm(ductal)
rm(lobular)

#### COADREAD: all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/coadread/coadread_rsem_unprocessed/'
colon <- read.csv(paste(path, 'colon_adenocarcinoma_rsem_unprocessed.csv', sep=''), sep='')
colon_muc <- read.csv(paste(path, 'colon_mucinous_adenocarcinoma_rsem_unprocessed.csv',sep=''),sep='')
rectal <- read.csv(paste(path, 'rectal_adenocarcinoma_rsem_unprocessed.csv',sep=''),sep='')
data <- cbind(colon, colon_muc, rectal)
# Define labels
y <- c(rep(1, dim(colon)[2]), rep(2, dim(colon_muc)[2]), rep(3, dim(rectal)[2]))
# Deifne prior
prior <- c(dim(colon)[2], dim(colon_muc)[2], dim(rectal)[2])/length(y)
#Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variable to save memory
rm(data)
rm(colon)
rm(colon_muc)
rm(rectal)

#### KIPAN: clear cell and papillary renal carcinoma ####
# kidney
path <- '../rnaseqV2data/unprocessed_rnaseqV2/kipan/kipan_rsem_unprocessed/'
clear <- read.csv(paste(path, 'kidney_clear_cell_renal_carcinoma_rsem_processed.csv', sep=''), sep='')
papillary <- read.csv(paste(path, 'kidney_papillary_renal_cell_carcinoma_rsem_processed.csv', sep=''),sep='')
data <- cbind(clear, papillary)
# Define labels
y <- c(rep(1, dim(clear)[2]), rep(2, dim(papillary)[2]))
# Define prior
prior <- c(dim(clear)[2], dim(papillary)[2])/length(y)
#Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
putput <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)

#### KIPAN: all labels ####
chromophobe <- read.csv(paste(path, 'kidney_chromophobe_rsem_processed.csv', sep=''), sep='')
data <- cbind(chromophobe, clear, papillary)
# Define labels
y <- c(rep(1, dim(chromophobe)[2]), rep(2, dim(clear)[2]), rep(3, dim(papillary)[2]))
# Define prior
prior <- c(dim(chromophobe)[2], dim(clear)[2], dim(papillary)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(clear)
rm(papillary)
rm(chromophobe)
rm(data)

#### LGG: all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/lgg/lgg_rsem_unprocessed/'
astro <- read.csv(paste(path, 'astrocytoma_rsem_unprocessed.csv', sep=''), sep='')
oligo_dendroglioma <- read.csv(paste(path, 'oligodendroglioma_rsem_unprocessed.csv', sep=''), sep='')
oligo_astro <- read.csv(paste(path, 'oligoastrocytoma_rsem_unprocessed.csv', sep=''), sep='')
data <- cbind(astro, oligo_astro, oligo_dendroglioma)
# Define labels
y <- c(rep(1, dim(astro)[2]), rep(2, dim(oligo_astro)[2]), rep(3, dim(oligo_dendroglioma)[2]))
# Define prior
prior <- c(dim(astro)[2], dim(oligo_astro)[2], dim(oligo_dendroglioma)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(astro)
rm(oligo_dendroglioma)
rm(oligo_astro)
rm(data)

#### MESO: epithelioid and biphasic malignant / all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/meso/meso_rsem_unprocessed/'
epithelioid <- read.csv(paste(path, 'epithelioid_mesothelioma_rsem_unprocessed.csv', sep=''), sep = '')
biphasic <- read.csv(paste(path, 'biphasic_mesothelioma_rsem_unprocessed.csv', sep=''), sep='')
data <- cbind(biphasic, epithelioid)
# Define labels
y <- c(rep(1, dim(biphasic)[2]), rep(2, dim(epithelioid)[2]))
# Define prior
prior <- c(dim(biphasic)[2], dim(epithelioid)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)
rm(epithelioid)
rm(biphasic)

#### SARCOMA: all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/sarc/sarc_rsem_unprocessed/'
leiomyosarcoma <- read.csv(paste(path, 'leiomyosarcoma_rsem_unprocessed.csv', sep=''), sep='')
liposarcoma <- read.csv(paste(path, 'dedifferentiated_liposarcoma_rsem_unprocessed.csv', sep=''), sep='')
pleomorphic <- read.csv(paste(path, 'pleomorphic_mfh_rsem_unprocessed.csv', sep=''), sep='')
myxofibrosarcoma <- read.csv(paste(path, 'myxofibrosarcoma_rsem_unprocessed.csv', sep=''), sep='')
undiff_pleomorphic <- read.csv(paste(path, 'undifferentiated_pleomorphic_sarcoma_rsem_unprocessed.csv', sep=''), sep='')
nerve <- read.csv(paste(path, 'malignant_peripheral_nerve_sheath_tumors_rsem_unprocessed.csv', sep=''), sep='')
#data <- cbind(leiomyosarcoma, liposarcoma, pleomorphic, myxofibrosarcoma, undiff_pleomorphic, nerve)
data <- cbind(liposarcoma, leiomyosarcoma, nerve, myxofibrosarcoma, pleomorphic, undiff_pleomorphic)
# Define labels
y <- c(rep(1, dim(liposarcoma)[2]), rep(2, dim(leiomyosarcoma)[2]), rep(3, dim(nerve)[2]), rep(4, dim(myxofibrosarcoma)[2]), 
       rep(5, dim(pleomorphic)[2]), rep(6, dim(undiff_pleomorphic)[2]))
# Define prior
prior <- c(dim(liposarcoma)[2], dim(leiomyosarcoma)[2], dim(nerve)[2], dim(myxofibrosarcoma)[2], 
           dim(pleomorphic)[2], dim(undiff_pleomorphic)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(leiomyosarcoma)
rm(liposarcoma)
rm(pleomorphic)
rm(myxofibrosarcoma)
rm(undiff_pleomorphic)
rm(nerve)
rm(data)

#### THCA: all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/thca/thca_rsem_unprocessed/'
classical <- read.csv(paste(path, 'thyroid_papillary_carcinoma_classical_rsem_unprocessed.csv', sep=''), sep='')
follicular <- read.csv(paste(path, 'thyroid_papillary_carcinoma_follicular_rsem_unprocessed.csv', sep=''), sep='')
tall_cell <- read.csv(paste(path, 'thyroid_papillary_carcinoma_tall_cell_rsem_unprocessed.csv', sep=''), sep='')
data <- cbind(classical, follicular, tall_cell)
# Define labels
y <- c(rep(1, dim(classical)[2]), rep(2, dim(follicular)[2]), rep(3, dim(tall_cell)[2]))
# Define prior
prior <- c(dim(classical)[2], dim(follicular)[2], dim(tall_cell)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(classical)
rm(follicular)
rm(tall_cell)
rm(data)

#### THYMOMA: all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/thym/thym_rsem_unprocessed/'
ab <- read.csv(paste(path, 'thymoma_type_ab_rsem_unprocessed.csv', sep=''), sep='')
b2 <- read.csv(paste(path, 'thymoma_type_b2_rsem_unprocessed.csv', sep=''), sep='')
a <- read.csv(paste(path, 'thymoma_type_a_rsem_unprocessed.csv', sep=''), sep='')
b1 <- read.csv(paste(path, 'thymoma_type_b1_rsem_unprocessed.csv', sep=''), sep='')
b3 <- read.csv(paste(path, 'thymoma_type_b3_rsem_unprocessed.csv', sep=''), sep='')
c <- read.csv(paste(path, 'thymoma_type_c_rsem_unprocessed.csv', sep=''), sep='')
data <- cbind(a, ab, b1, b2, b3, c)
# Define labels
y <- c(rep(1, dim(a)[2]), rep(2, dim(ab)[2]), rep(3, dim(b1)[2]), rep(4, dim(b2)[2]), rep(5, dim(b3)[2]), rep(6, dim(c)[2]))
# Define prior
prior <- c(dim(a)[2], dim(ab)[2], dim(b1)[2], dim(b2)[2], dim(b3)[2], dim(c)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(a)
rm(ab)
rm(b1)
rm(b2)
rm(b3)
rm(c)
rm(data)

#### UCEC: all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/ucec/ucec_rsem_unprocessed/'
endometrioid <- read.csv(paste(path, 'endometrioid_endometrial_adenocarcinoma_rsem_unprocessed.csv', sep=''),sep='')
serous <- read.csv(paste(path, 'serous_endometrial_adenocarcinoma_rsem_unprocessed.csv', sep=''), sep='')
mixed <- read.csv(paste(path, 'mixed_serous_and_endometrioid_rsem_unprocessed.csv', sep=''), sep='')
data <- cbind(endometrioid, mixed, serous)
# Define labels
y <- c(rep(1, dim(endometrioid)[2]), rep(2, dim(mixed)[2]), rep(3, dim(serous)[2]))
# Define prior
prior <- c(dim(endometrioid)[2], dim(mixed)[2], dim(serous)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(endometrioid)
rm(serous)
rm(mixed)
rm(data)

#### UCS: all labels ####
path <- '../rnaseqV2data/unprocessed_rnaseqV2/usc/ucs_rsem_unprocessed/'
mmmt <- read.csv(paste(path, 'malignant_mixed_mullerian_tumor_rsem_unprocessed.csv', sep=''), sep='')
mmmt_heterologous <- read.csv(paste(path, 'mmmt_heterologous_type_rsem_unprocessed.csv', sep=''), sep='')
mmmt_homologous <- read.csv(paste(path, 'mmmt_homologous_type_rsem_unprocessed.csv', sep=''), sep='')
data1 <- cbind(mmmt, mmmt_heterologous, mmmt_homologous)
# Define labels
y <- c(rep(1, dim(mmmt)[2]), rep(2, dim(mmmt_heterologous)[2]), rep(3, dim(mmmt_homologous)[2]))
# Define prior
prior <- c(dim(mmmt)[2], dim(mmmt_heterologous)[2], dim(mmmt_homologous)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(mmmt)
rm(mmmt_heterologous)
rm(mmmt_homologous)
rm(data)

##############################################
########## DE DATA ##########################
##############################################

#### BRCA: ductal and lobular / all labels ####
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'BRCA112samples_DESeq2_DEgenes.csv', sep=''),sep=',')
genes <- as.character(csv$Gene)
# Read and combine BRCA data
path <- '../rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/'
ductal <- read.csv(paste(path, 'brca_ductal_rsem_processed.csv', sep=''), sep='')
lobular <- read.csv(paste(path, 'brca_lobular_rsem_processed.csv', sep=''), sep='')
data <- cbind(ductal, lobular)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1,dim(ductal)[2]), rep(2,dim(lobular)[2]))
# Define prior
prior <- c(dim(ductal)[2], dim(lobular)[2])/length(y)
# Compute prediction accuracy
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)
rm(csv)
rm(ductal)
rm(lobular)

#### COADREAD: all labels ####
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'coadread_DESeq2_DEgenes.csv',sep=''),sep=',')
genes <- as.character(csv$Gene)
path <- '../rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed/'
colon <- read.csv(paste(path, 'colon_adenocarcinoma_raw.csv', sep=''), sep='')
colon_muc <- read.csv(paste(path, 'colon_mucinous_adenocarcinoma_raw.csv',sep=''),sep='')
rectal <- read.csv(paste(path, 'rectal_adenocarcinoma_raw.csv',sep=''),sep='')
data <- cbind(colon, colon_muc, rectal)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(colon)[2]), rep(2, dim(colon_muc)[2]), rep(3, dim(rectal)[2]))
# Deifne prior
prior <- c(dim(colon)[2], dim(colon_muc)[2], dim(rectal)[2])/length(y)
#Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variable to save memory
rm(data)
rm(csv)
rm(colon)
rm(colon_muc)
rm(rectal)

#### KIPAN: clear cell and papillary renal carcinoma ####
# kidney
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'kipan_DESeq2_DEgenes.csv',sep=''),sep=',')
genes <- as.character(csv$Gene)
path <- '../rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/'
clear <- read.csv(paste(path, 'kidney_clear_cell_renal_carcinoma_rsem_processed.csv', sep=''), sep='')
papillary <- read.csv(paste(path, 'kidney_papillary_renal_cell_carcinomarsem_processed_log2.csv', sep=''),sep='')
data <- cbind(clear, papillary)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(clear)[2]), rep(2, dim(papillary)[2]))
# Define prior
prior <- c(dim(clear)[2], dim(papillary)[2])/length(y)
#Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
putput <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(data)

#### KIPAN: all labels ####
chromophobe <- read.csv(paste(path, 'kidney_chromophobe_rsem_processed.csv', sep=''), sep='')
data <- cbind(chromophobe, clear, papillary)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(chromophobe)[2]), rep(2, dim(clear)[2]), rep(3, dim(papillary)[2]))
# Define prior
prior <- c(dim(chromophobe)[2], dim(clear)[2], dim(papillary)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(clear)
rm(papillary)
rm(chromophobe)
rm(data)

#### LGG: all labels ####
# IS THIS CORRECT FILE??? LUAD??
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'luad_DESeq2_DEgenes.csv',sep=''),sep=',')
genes <- as.character(csv$Gene)
path <- '../rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/'
astro <- read.csv(paste(path, 'astrocytoma_rsem_processed_log2.csv', sep=''), sep='')
oligo_dendroglioma <- read.csv(paste(path, 'oligodendroglioma_rsem_processed.csv', sep=''), sep='')
oligo_astro <- read.csv(paste(path, 'oligoastrocytoma_rsem_processed.csv', sep=''), sep='')
data <- cbind(astro, oligo_astro, oligo_dendroglioma)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(astro)[2]), rep(2, dim(oligo_astro)[2]), rep(3, dim(oligo_dendroglioma)[2]))
# Define prior
prior <- c(dim(astro)[2], dim(oligo_astro)[2], dim(oligo_dendroglioma)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(astro)
rm(oligo_dendroglioma)
rm(oligo_astro)
rm(data)
rm(csv)

#### SARCOMA: all labels ####
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'sarc_DESeq2_DEgenes.csv',sep=''),sep=',')
genes <- as.character(csv$Gene)
path <- '../rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/'
leiomyosarcoma <- read.csv(paste(path, 'leiomyosarcoma_raw.csv', sep=''), sep='')
liposarcoma <- read.csv(paste(path, 'dedifferentiated_liposarcoma_raw.csv', sep=''), sep='')
pleomorphic <- read.csv(paste(path, 'pleomorphic_mfh_raw.csv', sep=''), sep='')
myxofibrosarcoma <- read.csv(paste(path, 'myxofibrosarcoma_raw.csv', sep=''), sep='')
undiff_pleomorphic <- read.csv(paste(path, 'undifferentiated_pleomorphic_sarcoma_raw.csv', sep=''), sep='')
nerve <- read.csv(paste(path, 'malignant_peripheral_nerve_sheath_tumors_raw.csv', sep=''), sep='')
data <- cbind(liposarcoma, leiomyosarcoma, nerve, myxofibrosarcoma, pleomorphic, undiff_pleomorphic)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(liposarcoma)[2]), rep(2, dim(leiomyosarcoma)[2]), rep(3, dim(nerve)[2]), rep(4, dim(myxofibrosarcoma)[2]), 
       rep(5, dim(pleomorphic)[2]), rep(6, dim(undiff_pleomorphic)[2]))
# Define prior
prior <- c(dim(liposarcoma)[2], dim(leiomyosarcoma)[2], dim(nerve)[2], dim(myxofibrosarcoma)[2], 
           dim(pleomorphic)[2], dim(undiff_pleomorphic)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(leiomyosarcoma)
rm(liposarcoma)
rm(pleomorphic)
rm(myxofibrosarcoma)
rm(undiff_pleomorphic)
rm(nerve)
rm(data)
rm(csv)

#### THCA: all labels ####
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'thca_DESeq2_DEgenes.csv',sep=''),sep=',')
genes <- as.character(csv$Gene)
path <- '../rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/'
classical <- read.csv(paste(path, 'thyroid_papillary_carcinoma_classical_rsem_processed.csv', sep=''), sep='')
follicular <- read.csv(paste(path, 'thyroid_papillary_carcinoma_follicular_rsem_processed.csv', sep=''), sep='')
tall_cell <- read.csv(paste(path, 'thyroid_papillary_carcinoma_tall_cell_rsem_processed.csv', sep=''), sep='')
data <- cbind(classical, follicular, tall_cell)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(classical)[2]), rep(2, dim(follicular)[2]), rep(3, dim(tall_cell)[2]))
# Define prior
prior <- c(dim(classical)[2], dim(follicular)[2], dim(tall_cell)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(classical)
rm(follicular)
rm(tall_cell)
rm(data)
rm(csv)

#### THYMOMA: all labels ####
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'thym_DESeq2_DEgenes.csv',sep=''),sep=',')
genes <- as.character(csv$Gene)
path <- '../rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/'
ab <- read.csv(paste(path, 'thymoma_type_ab_rsem_processed.csv', sep=''), sep='')
b2 <- read.csv(paste(path, 'thymoma_type_b2_rsem_processed_log2.csv', sep=''), sep='')
a <- read.csv(paste(path, 'thymoma_type_a_rsem_processed.csv', sep=''), sep='')
b1 <- read.csv(paste(path, 'thymoma_type_b1_rsem_processed.csv', sep=''), sep='')
b3 <- read.csv(paste(path, 'thymoma_type_b3_rsem_processed_log2.csv', sep=''), sep='')
c <- read.csv(paste(path, 'thymoma_type_c_raw.csv', sep=''), sep='')
data <- cbind(a, ab, b1, b2, b3, c)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(a)[2]), rep(2, dim(ab)[2]), rep(3, dim(b1)[2]), rep(4, dim(b2)[2]), rep(5, dim(b3)[2]), rep(6, dim(c)[2]))
# Define prior
prior <- c(dim(a)[2], dim(ab)[2], dim(b1)[2], dim(b2)[2], dim(b3)[2], dim(c)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(a)
rm(ab)
rm(b1)
rm(b2)
rm(b3)
rm(c)
rm(data)
rm(csv)

#### UCEC: all labels ####
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'ucec_DESeq2_DEgenes.csv',sep=''),sep=',')
genes <- as.character(csv$Gene)
path <- '../rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/'
endometrioid <- read.csv(paste(path, 'endometrioid_endometrial_adenocarcinoma_rsem_processed.csv', sep=''),sep='')
serous <- read.csv(paste(path, 'serous_endometrial_adenocarcinoma_rsem_processed.csv', sep=''), sep='')
mixed <- read.csv(paste(path, 'mixed_serous_and_endometrioid_rsem_processed.csv', sep=''), sep='')
data <- cbind(endometrioid, mixed, serous)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(endometrioid)[2]), rep(2, dim(mixed)[2]), rep(3, dim(serous)[2]))
# Define prior
prior <- c(dim(endometrioid)[2], dim(mixed)[2], dim(serous)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(endometrioid)
rm(serous)
rm(mixed)
rm(data)
rm(csv)


#### UCS: all labels ####
path <- '../DE genes csvs DESeq2/'
csv <- read.csv(paste(path, 'lusc_DESeq2_DEgenes.csv',sep=''),sep=',')
genes <- as.character(csv$Gene)
path <- '../rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/'
mmmt <- read.csv(paste(path, 'malignant_mixed_mullerian_tumor_rsem_processed.csv', sep=''), sep='')
mmmt <- mmmt[1:14975,]
mmmt_heterologous <- read.csv(paste(path, 'mmmt_heterologous_type_rsem_processed.csv', sep=''), sep='')
mmmt_homologous <- read.csv(paste(path, 'mmmt_homologous_type_rsem_processed.csv', sep=''), sep='')
data <- cbind(mmmt, mmmt_heterologous, mmmt_homologous)
idx <- which(rownames(data)%in%genes)
data <- data[idx,]
# Define labels
y <- c(rep(1, dim(mmmt)[2]), rep(2, dim(mmmt_heterologous)[2]), rep(3, dim(mmmt_homologous)[2]))
# Define prior
prior <- c(dim(mmmt)[2], dim(mmmt_heterologous)[2], dim(mmmt_homologous)[2])/length(y)
# Compute predictions
output <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
output <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
# Remove some variables to save memory
rm(mmmt)
rm(mmmt_heterologous)
rm(mmmt_homologous)
rm(data)
rm(csv)


