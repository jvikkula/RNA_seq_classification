require(caret)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('PLDA.R')
source('ziplda.R')
library(DESeq2)

get.accuracy <- function(data, y, func = c('PLDA','ZIPLDA'), prior, method){
	set.seed(123)
	flds <- createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)
	acc <- c()
	ypred <- matrix(NA, nrow=1, ncol=length(y))
	ycorrect <- matrix(NA, nrow=1, ncol=length(y))
	start <- 1
	stop <- 0
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
		tb <- table(ytrue, out$yhat)
		acc[k] <- sum(diag(tb))/sum(tb)
		stop <- stop + length(out$yhat)
		ypred[start:stop] <- out$yhat
		ycorrect[start:stop] <- ytrue
		start <- start + length(out$yhat)
	}
	cat('Accuracies are:', acc, '\n')
  cat('Mean accuracy with 5-fold cross validation:', round(mean(acc)*100,2), '% \n')

	confusion <- confusionMatrix(factor(ycorrect), factor(ypred), dnn = c('Reference', 'Prediction'))
	cat('Accuracy from confusion matrix:', confusion$overall[1])

	if (dim(table(ycorrect, ypred))[1]==2){
    sensitivity <- confusion$byClass[1]
    specificity <- confusion$byClass[2]		
	}
	else if (dim(table(ycorrect, ypred))[1]==3){
		sensitivity <- sum(confusion$byClass[1:3])/3
		specificity <- sum(confusion$byClass[4:6])/3
	
	}
	else if (dim(table(ycorrect, ypred))[1]==6){
		sensitivity <- sum(confusion$byClass[1:6])/6
		specificity <- sum(confusion$byClass[7:12])/6
	}

	print('Confusion matrix:')
	print(confusion)
	cat('Sensitivity is:', sensitivity, '\n')
	cat('Specificity is:', specificity, '\n')
	return(list(ytrue=ycorrect, ypred=ypred, acc=acc, confusion=confusion))
}


print('BRCA \n')
path <- 'rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/'
ductal <- read.csv(paste(path, 'brca_ductal_rsem_processed.csv', sep=''), sep='')
lobular <- read.csv(paste(path, 'brca_lobular_rsem_processed.csv', sep=''), sep='')
data <- cbind(ductal,lobular)
y <- c(rep(1,dim(ductal)[2]), rep(2,dim(lobular)[2]))
prior <- c(dim(ductal)[2], dim(lobular)[2])/length(y)
print('Perform PLDA...\n')
brca_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(brca_plda, file='../results/brca_plda.Rdata')
print('Perform ZIPLDA...\n')
brca_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(brca_ziplda, file='../results/brca_ziplda.Rdata')
print('BRCA done!\n')
rm(data)
rm(ductal)
rm(lobular)


print('COADREAD:\n')
path <- 'rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed/'
colon <- read.csv(paste(path, 'colon_adenocarcinoma_raw.csv', sep=''), sep='')
colon_muc <- read.csv(paste(path, 'colon_mucinous_adenocarcinoma_raw.csv',sep=''),sep='')
rectal <- read.csv(paste(path, 'rectal_adenocarcinoma_raw.csv',sep=''),sep='')
data <- cbind(colon, colon_muc, rectal)
y <- c(rep(1, dim(colon)[2]), rep(2, dim(colon_muc)[2]), rep(3, dim(rectal)[2]))
prior <- c(dim(colon)[2], dim(colon_muc)[2], dim(rectal)[2])/length(y)
print('Perform PLDA...\n')
coadread_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(coadread_plda, file='../results/coadread_plda.Rdata')
print('Perform ZIPLDA...\n')
coadread_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(coadread_ziplda, file='../results/coadread_ziplda.Rdata')
print('COADREAD DONE!\n')
rm(colon)
rm(colon_muc)
rm(rectal)
rm(data)


print('KIPAN:\n')
path <- 'rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/'
clear <- read.csv(paste(path, 'kidney_clear_cell_renal_carcinoma_rsem_processed.csv', sep=''), sep='')
papillary <- read.csv(paste(path, 'kidney_papillary_renal_cell_carcinomarsem_processed_log2.csv', sep=''),sep='')
chromophobe <- read.csv(paste(path, 'kidney_chromophobe_rsem_processed.csv', sep=''), sep='')
data <- cbind(chromophobe, clear, papillary)
y <- c(rep(1, dim(chromophobe)[2]), rep(2, dim(clear)[2]), rep(3, dim(papillary)[2]))
prior <- c(dim(chromophobe)[2], dim(clear)[2], dim(papillary)[2])/length(y)
print('Perform PLDA...\n')
kipan_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(kipan_plda, file='../results/kipan_plda.Rdata')
print('Perform ZIPLDA...\n')
kipan_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(kipan_ziplda, file='../results/kipan_ziplda.Rdata')
rm(clear)
rm(papillary)
rm(chromophobe)
rm(data)
print('KIPAN DONE!\n')

print('LGG')
path <- 'rnaseqV2data/processed_rnaseqV2/lgg/lgg_rsem_processed/'
astro <- read.csv(paste(path, 'astrocytoma_rsem_processed_log2.csv', sep=''), sep='')
oligo_dendroglioma <- read.csv(paste(path, 'oligodendroglioma_rsem_processed.csv', sep=''), sep='')
oligo_astro <- read.csv(paste(path, 'oligoastrocytoma_rsem_processed.csv', sep=''), sep='')
data <- cbind(astro, oligo_astro, oligo_dendroglioma)
# Define labels
y <- c(rep(1, dim(astro)[2]), rep(2, dim(oligo_astro)[2]), rep(3, dim(oligo_dendroglioma)[2]))
# Define prior
prior <- c(dim(astro)[2], dim(oligo_astro)[2], dim(oligo_dendroglioma)[2])/length(y)
print('Perform PLDA...\n')
lgg_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(lgg_plda, file='../results/lgg_plda.Rdata')
print('Perform ZIPLDA...\n')
lgg_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(lgg_ziplda, file='../results/lgg_ziplda.Rdata')
rm(astro)
rm(oligo_dendroglioma)
rm(oligo_astro)
rm(data)
print('LGG DONE!\n')

print('MESO')
path <- 'rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/'
epithelioid <- read.csv(paste(path, 'epithelioid_mesothelioma_rsem_processed.csv', sep=''), sep = '')
biphasic <- read.csv(paste(path, 'biphasic_mesothelioma_rsem_processed.csv', sep=''), sep='')
data <- cbind(biphasic, epithelioid)
# Define labels
y <- c(rep(1, dim(biphasic)[2]), rep(2, dim(epithelioid)[2]))
# Define prior
prior <- c(dim(biphasic)[2], dim(epithelioid)[2])/length(y)
print('Perform PLDA...\n')
meso_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(meso_plda, file='../results/meso_plda.Rdata')
print('Perform ZIPLDA...\n')
meso_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(meso_ziplda, file='../results/meso_ziplda.Rdata')
rm(data)
rm(epithelioid)
rm(biphasic)
print('MESO DONE!\n')

print('SARCOMA')
path <- 'rnaseqV2data/processed_rnaseqV2/sarcoma/sarcoma_rsem_processed/'
leiomyosarcoma <- read.csv(paste(path, 'leiomyosarcoma_raw.csv', sep=''), sep='')
liposarcoma <- read.csv(paste(path, 'dedifferentiated_liposarcoma_raw.csv', sep=''), sep='')
pleomorphic <- read.csv(paste(path, 'pleomorphic_mfh_raw.csv', sep=''), sep='')
myxofibrosarcoma <- read.csv(paste(path, 'myxofibrosarcoma_raw.csv', sep=''), sep='')
undiff_pleomorphic <- read.csv(paste(path, 'undifferentiated_pleomorphic_sarcoma_raw.csv', sep=''), sep='')
nerve <- read.csv(paste(path, 'malignant_peripheral_nerve_sheath_tumors_raw.csv', sep=''), sep='')
data <- cbind(liposarcoma, leiomyosarcoma, nerve, myxofibrosarcoma, pleomorphic, undiff_pleomorphic)
# Define labels
y <- c(rep(1, dim(liposarcoma)[2]), rep(2, dim(leiomyosarcoma)[2]), rep(3, dim(nerve)[2]), rep(4, dim(myxofibrosarcoma)[2]), 
              rep(5, dim(pleomorphic)[2]), rep(6, dim(undiff_pleomorphic)[2]))
# Define prior
prior <- c(dim(liposarcoma)[2], dim(leiomyosarcoma)[2], dim(nerve)[2], dim(myxofibrosarcoma)[2], 
                     dim(pleomorphic)[2], dim(undiff_pleomorphic)[2])/length(y)
print('Perform PLDA...\n')
sarcoma_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(sarcoma_plda, file='../results/sarcoma_plda.Rdata')
print('Perform ZIPLDA...\n')
sarcoma_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(sarcoma_ziplda, file='../results/sarcoma_ziplda.Rdata')
rm(leiomyosarcoma)
rm(liposarcoma)
rm(pleomorphic)
rm(myxofibrosarcoma)
rm(undiff_pleomorphic)
rm(nerve)
rm(data)
print('SARCOMA DONE!\n')

print('THCA')
path <- 'rnaseqV2data/processed_rnaseqV2/thca/thca_rsem_processed/'
classical <- read.csv(paste(path, 'thyroid_papillary_carcinoma_classical_rsem_processed.csv', sep=''), sep='')
follicular <- read.csv(paste(path, 'thyroid_papillary_carcinoma_follicular_rsem_processed.csv', sep=''), sep='')
tall_cell <- read.csv(paste(path, 'thyroid_papillary_carcinoma_tall_cell_rsem_processed.csv', sep=''), sep='')
data <- cbind(classical, follicular, tall_cell)
# Define labels
y <- c(rep(1, dim(classical)[2]), rep(2, dim(follicular)[2]), rep(3, dim(tall_cell)[2]))
# Define prior
prior <- c(dim(classical)[2], dim(follicular)[2], dim(tall_cell)[2])/length(y)
print('Perform PLDA...\n')
thca_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(thca_plda, file='../results/thca_plda.Rdata')
print('Perform ZIPLDA...\n')
thca_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(thca_ziplda, file='../results/thca_ziplda.Rdata')
rm(classical)
rm(follicular)
rm(tall_cell)
rm(data)
print('THCA DONE!\n')

print('THYMOMA')
path <- 'rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/'
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
print('Perform PLDA...\n')
thymoma_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(thymoma_plda, file='../results/thymoma_plda.Rdata')
print('Perform ZIPLDA...\n')
thymoma_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(thymoma_ziplda, file='../results/thymoma_ziplda.Rdata')
rm(a)
rm(ab)
rm(b1)
rm(b2)
rm(b3)
rm(c)
rm(data)
print('THYMOMA DONE!\n')

print('UCEC')
path <- 'rnaseqV2data/processed_rnaseqV2/ucec/ucec_rsem_processed/'
endometrioid <- read.csv(paste(path, 'endometrioid_endometrial_adenocarcinoma_rsem_processed.csv', sep=''),sep='')
serous <- read.csv(paste(path, 'serous_endometrial_adenocarcinoma_rsem_processed.csv', sep=''), sep='')
mixed <- read.csv(paste(path, 'mixed_serous_and_endometrioid_rsem_processed.csv', sep=''), sep='')
data <- cbind(endometrioid, mixed, serous)
# Define labels
y <- c(rep(1, dim(endometrioid)[2]), rep(2, dim(mixed)[2]), rep(3, dim(serous)[2]))
# Define prior
prior <- c(dim(endometrioid)[2], dim(mixed)[2], dim(serous)[2])/length(y)
print('Perform PLDA...\n')
ucec_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(ucec_plda, file='../results/ucec_plda.Rdata')
print('Perform ZIPLDA...\n')
ucec_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(ucec_ziplda, file='../results/ucec_ziplda.Rdata')
rm(endometrioid)
rm(serous)
rm(mixed)
rm(data)
print('UCEC DONE!\n')


print('UCS')
path <- 'rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/'
mmmt <- read.csv(paste(path, 'malignant_mixed_mullerian_tumor_rsem_processed.csv', sep=''), sep='')
mmmt <- mmmt[1:14975,]
mmmt_heterologous <- read.csv(paste(path, 'mmmt_heterologous_type_rsem_processed.csv', sep=''), sep='')
mmmt_homologous <- read.csv(paste(path, 'mmmt_homologous_type_rsem_processed.csv', sep=''), sep='')
data <- cbind(mmmt, mmmt_heterologous, mmmt_homologous)
# Define labels
y <- c(rep(1, dim(mmmt)[2]), rep(2, dim(mmmt_heterologous)[2]), rep(3, dim(mmmt_homologous)[2]))
# Define prior
prior <- c(dim(mmmt)[2], dim(mmmt_heterologous)[2], dim(mmmt_homologous)[2])/length(y)
print('Perform PLDA...\n')
ucs_plda <- get.accuracy(data, y, 'PLDA', prior, 'deseq')
save(ucs_plda, file='../results/ucs_plda.Rdata')
print('Perform ZIPLDA...\n')
ucs_ziplda <- get.accuracy(data, y, 'ZIPLDA', prior, 'deseq')
save(ucs_ziplda, file='../results/ucs_ziplda.Rdata')
rm(mmmt)
rm(mmmt_heterologous)
rm(mmmt_homologous)
rm(data)
print('UCS DONE!\n')


#################################
#### PLOT CONFUSION MATRICES ####
#################################
library("dplyr")
library("tidyr")
library('ggplot2')

#### PLDA ####

# BRCA #
load('../data/brca_plda.Rdata')
tb <- table(brca_plda$ytrue, brca_plda$ypred)
df <- data.frame(reference=c('ductal','lobular','ductal','lobular'), 
                 prediction=c('ductal','ductal','lobular','lobular'), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 79.59 %\nSensitivity: 94.76 %\nSpecificity: 50.30 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# COADREAD #
load('../data/coadread_plda.Rdata')
tb <- table(coadread_plda$ytrue, coadread_plda$ypred)
df <- data.frame(reference=rep(c('adeno','mucinous','rectal'),3), 
                 prediction=c(rep('adeno',3), rep('mucinous',3), rep('rectal',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 48.7 %\nSensitivity: 46.53 %\nSpecificity: 72.66 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# KIPAN #
load('../data/kipan_plda.Rdata')
tb <- table(kipan_plda$ytrue, kipan_plda$ypred)
df <- data.frame(reference=rep(c('chromophobe','clear cell','papillary'),3), 
                 prediction=c(rep('chromophobe',3), rep('clear cell',3), rep('papillary',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 84.59 %\nSensitivity: 75.36 %\nSpecificity: 90.56 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# LGG #
load('../data/lgg_plda.Rdata')
tb <- table(lgg_plda$ytrue, lgg_plda$ypred)
df <- data.frame(reference=rep(c('astro','oligoastro','oligodendro'),3), 
                 prediction=c(rep('astro',3), rep('oligoastro',3), rep('oligodendro',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 54.56 %\nSensitivity: 53.45 %\nSpecificity: 77.36 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# MESO 
load('../data/meso_plda.Rdata')
tb <- table(meso_plda$ytrue, meso_plda$ypred)
df <- data.frame(reference=c('biphasic','epithelioid','biphasic','epithelioid'), 
                 prediction=c('biphasic','biphasic','epithelioid','epithelioid'), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 80 %\nSensitivity: 62.96 %\nSpecificity: 88.68 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# SARCOMA 
load('../data/sarcoma_plda.Rdata')
tb <- table(sarcoma_plda$ytrue, sarcoma_plda$ypred)
df <- data.frame(reference=rep(c('delipo','leiomyo','mpnst','myxofibro',
                                 'pleomorphic','ups'),6), 
                 prediction=c(rep('delipo',6), rep('leiomyo',6), rep('mpnst',6), rep('myxo',6),
                              rep('pleom', 6), rep('ups',6)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 64.63 %\nSensitivity: 52.09 %\nSpecificity: 92.53 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=6)

# THCA #
load('../data/thca_plda.Rdata')
tb <- table(thca_plda$ytrue, thca_plda$ypred)
df <- data.frame(reference=rep(c('classical','follicular','tall cell'),3), 
                 prediction=c(rep('classical',3), rep('follicular',3), rep('tall cell',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 51.01 %\nSensitivity: 49.11%\nSpecificity: 74.79%') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# THYMOMA #
load('../data/thymoma_plda.Rdata')
tb <- table(thymoma_plda$ytrue, thymoma_plda$ypred)
df <- data.frame(reference=rep(c('a','ab','b1','b2','b3','c'),6), 
                 prediction=c(rep('a',6), rep('ab',6), rep('b1',6), rep('b2',6), rep('b3', 6), rep('c',6)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 59.17 %\nSensitivity: 62.69%\nSpecificity: 91.64%') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=6)

# UCEC # 
load('../data/ucec_plda.Rdata')
tb <- table(ucec_plda$ytrue, ucec_plda$ypred)
df <- data.frame(reference=rep(c('endometrioid','mixed','serous'),3), 
                 prediction=c(rep('endometrioid',3), rep('mixed',3), rep('serous',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 73.22 %\nSensitivity: 52.35 %\nSpecificity: 83.93 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# UCS #
load('../data/ucs_plda.Rdata')
tb <- table(ucs_plda$ytrue, ucs_plda$ypred)
df <- data.frame(reference=rep(c('mmmt','mmmt hetero','mmmt homo'),3), 
                 prediction=c(rep('mmmt',3), rep('mmmt hetero',3), rep('mmmt homo',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 47.37 %\nSensitivity: 47.67 %\nSpecificity: 72.97 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

####### ZIPLDA ##########

# BRCA #
load('../data/brca_ziplda.Rdata')
tb <- table(brca_ziplda$ytrue, brca_ziplda$ypred)
df <- data.frame(reference=c('ductal','lobular','ductal','lobular'), 
                 prediction=c('ductal','ductal','lobular','lobular'), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 79.7 %\nSensitivity: 94.77 %\nSpecificity: 50.45 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)


# COADREAD #
load('../data/coadread_ziplda.Rdata')
tb <- table(coadread_ziplda$ytrue, coadread_ziplda$ypred)
df <- data.frame(reference=rep(c('adeno','mucinous','rectal'),3), 
                 prediction=c(rep('adeno',3), rep('mucinous',3), rep('rectal',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 48.97 %\nSensitivity: 56.72 %\nSpecificity: 72.72 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# KIPAN #
load('../data/kipan_ziplda.Rdata')
tb <- table(kipan_ziplda$ytrue, kipan_ziplda$ypred)
df <- data.frame(reference=rep(c('chromophobe','clear cell','papillary'),3), 
                 prediction=c(rep('chromophobe',3), rep('clear cell',3), rep('papillary',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 84.59 %\nSensitivity: 75.36 %\nSpecificity: 90.56 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# LGG #
load('../data/lgg_ziplda.Rdata')
tb <- table(lgg_ziplda$ytrue, lgg_ziplda$ypred)
df <- data.frame(reference=rep(c('astro','oligoastro','oligodendro'),3), 
                 prediction=c(rep('astro',3), rep('oligoastro',3), rep('oligodendro',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 54.76 %\nSensitivity: 53.57 %\nSpecificity: 77.46 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# MESO 
load('../data/meso_ziplda.Rdata')
tb <- table(meso_ziplda$ytrue, meso_ziplda$ypred)
df <- data.frame(reference=c('biphasic','epithelioid','biphasic','epithelioid'), 
                 prediction=c('biphasic','biphasic','epithelioid','epithelioid'), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 80 %\nSensitivity: 62.96 %\nSpecificity: 88.68 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# SARCOMA 
load('../data/sarcoma_ziplda.Rdata')
tb <- table(sarcoma_ziplda$ytrue, sarcoma_ziplda$ypred)
df <- data.frame(reference=rep(c('delipo','leiomyo','mpnst','myxofibro',
                                 'pleomorphic','ups'),6), 
                 prediction=c(rep('delipo',6), rep('leiomyo',6), rep('mpnst',6), rep('myxo',6),
                              rep('pleom', 6), rep('ups',6)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 64.23 %\nSensitivity: 51.63 %\nSpecificity: 92.45 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=6)

# THCA #
load('../data/thca_ziplda.Rdata')
tb <- table(thca_ziplda$ytrue, thca_ziplda$ypred)
df <- data.frame(reference=rep(c('classical','follicular','tall cell'),3), 
                 prediction=c(rep('classical',3), rep('follicular',3), rep('tall cell',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 51.01 %\nSensitivity: 49.11%\nSpecificity: 74.79%') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# THYMOMA #
load('../data/thymoma_ziplda.Rdata')
tb <- table(thymoma_ziplda$ytrue, thymoma_ziplda$ypred)
df <- data.frame(reference=rep(c('a','ab','b1','b2','b3','c'),6), 
                 prediction=c(rep('a',6), rep('ab',6), rep('b1',6), rep('b2',6), rep('b3', 6), rep('c',6)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 59.17 %\nSensitivity: 62.69%\nSpecificity: 91.64%') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=6)

# UCEC # 
load('../data/ucec_ziplda.Rdata')
tb <- table(ucec_ziplda$ytrue, ucec_ziplda$ypred)
df <- data.frame(reference=rep(c('endometrioid','mixed','serour'),3), 
                 prediction=c(rep('endometrioid',3), rep('mixed',3), rep('serous',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 73.22 %\nSensitivity: 52.35 %\nSpecificity: 83.93 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)

# UCS #
load('../data/ucs_ziplda.Rdata')
tb <- table(ucs_ziplda$ytrue, ucs_ziplda$ypred)
df <- data.frame(reference=rep(c('mmmt','mmmt hetero','mmmt homo'),3), 
                 prediction=c(rep('mmmt',3), rep('mmmt hetero',3), rep('mmmt homo',3)), 
                 value=c(tb))
df <- df %>% mutate(prediction = factor(prediction), reference = factor(reference, levels=rev(unique(reference))))
ggplot(df, aes(x=prediction, y=reference, fill=value)) +
  geom_tile() + theme_bw() + coord_equal() + 
  scale_fill_distiller(palette="Blues", direction=1) + 
  guides(fill=F) + # removing legend for `fill`
  ggtitle('Accuracy: 47.37 %\nSensitivity: 47.67 %\nSpecificity: 72.97 %') + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  geom_text(aes(label=value), color="black", size=10)


