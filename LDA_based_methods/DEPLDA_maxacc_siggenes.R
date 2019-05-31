## TO run DEPLDA with optimizatoin method 'max accuracy' or 'significant genes', run first deplda_increase_accuracy.R to
## get the significantly differentially expressed genes. Then you can run this file.

require(caret)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('PLDA.R')
source('ziplda.R')
library(DESeq2)


simplify.deplda <- function(full_data, full_y, prior, method, genes, n_top_acc, optimizatin='significant_genes'){
  
  if (optimization=='max_accuracy'){
  # find max accuracy for each fold and save genes until that genes that gives max accuracy
	max_idx <- which(n_top_acc[1,]==max(n_top_acc[1,]))
	if (length(max_idx)>1){
		max_idx <- max_idx[1]
	}
	top_genes1 <- genes[1,1:(max_idx+1)]

  max_idx <- which(n_top_acc[2,]==max(n_top_acc[2,]))
  if (length(max_idx)>1){
    max_idx <- max_idx[1]
  }	
	top_genes2 <- genes[2,1:(max_idx+1)]

  max_idx <- which(n_top_acc[3,]==max(n_top_acc[3,]))
  if (length(max_idx)>1){
    max_idx <- max_idx[1]
  }	
	top_genes3 <- genes[3,1:(max_idx+1)]

  max_idx <- which(n_top_acc[4,]==max(n_top_acc[4,]))
  if (length(max_idx)>1){
    max_idx <- max_idx[1]
  }	
	top_genes4 <- genes[4,1:(max_idx+1)]

  max_idx <- which(n_top_acc[5,]==max(n_top_acc[5,]))
  if (length(max_idx)>1){
    max_idx <- max_idx[1]
  }
	top_genes5 <- genes[5,1:(max_idx+1)]
  }
  
  if (optimization=='significant_genes'){
    top_genes1 <- genes[1,]
    top_genes2 <- genes[2,]
    top_genes3 <- genes[3,]
    top_genes4 <- genes[4,]
    top_genes5 <- genes[5,]
  }

	final_genes <- Reduce(intersect, list(top_genes1, top_genes2, top_genes3, top_genes4, top_genes5))
	if (length(final_genes)<=1){
		final_genes <- unique(c(top_genes1, top_genes2, top_genes3, top_genes4, top_genes5))
		print('Overlap not used!!')
	}
	cat('Number of final genes is:', length(final_genes), '\n')

 	# leave 20% of the data for testing as previously
  set.seed(12345)
  folds <- createFolds(full_y, k=5, list=TRUE, returnTrain=FALSE)
  data <- full_data[,-folds[[1]]]
  test_data <- full_data[,folds[[1]]]
  y <- full_y[-folds[[1]]]
  test_y <- full_y[folds[[1]]]

	all_genes <- rownames(data)
        test <- t(test_data[which(all_genes%in%final_genes),])
        train <- t(data[which(all_genes%in%final_genes),])
        ytrue <- test_y
        ytrain <- y

	out <- PLDA(train, test, ytrain, prior=prior, method=method)

	tb <- table(ytrue, out$yhat)
	accuracy <- sum(diag(tb))/sum(tb)
	cat('Accuracy is', accuracy, '\n')


  confusion <- confusionMatrix(factor(ytrue), factor(out$yhat), dnn=c('Reference', 'Prediction'))
  cat('Accuracy from confusion matrix:', round(confusion$overall[1]*100,2), '\n')


  if (dim(table(ytrue, out$yhat))[1]==2){
    sensitivity <- confusion$byClass[1]
    specificity <- confusion$byClass[2]
  }
  else if (dim(table(ytrue, out$yhat))[1]==3){
          sensitivity <- sum(confusion$byClass[1:3])/3
          specificity <- sum(confusion$byClass[4:6])/3
  }
  else if (dim(table(ytrue, out$yhat))[1]==6){
    sensitivity <- sum(confusion$byClass[1:6])/6
    specificity <- sum(confusion$byClass[7:12])/6
  }	

  print(confusion)
	cat('Sensitivity is:', sensitivity, '\n')
  cat('Specificity is:', specificity, '\n')

	return(list(confusion=confusion, final_genes=final_genes, ytrue=ytrue, ypred=out$yhat))
}


print('BRCA')
path <- 'rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/'
ductal <- read.csv(paste(path, 'brca_ductal_rsem_processed.csv', sep=''), sep='')
lobular <- read.csv(paste(path, 'brca_lobular_rsem_processed.csv', sep=''), sep='')
data <- cbind(ductal,lobular)
y <- c(rep(1,dim(ductal)[2]), rep(2,dim(lobular)[2]))
prior <- c(dim(ductal)[2], dim(lobular)[2])/length(y)
load('results/results_brca_deplda.Rdata')
genes <- results_brca_deplda$genes
n_top_acc <- results_brca_deplda$n_top_acc
brca <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(brca, file='results/simpde_brca.Rdata')
rm(data)
rm(ductal)
rm(lobular)


print('COADREAD:')
path <- 'rnaseqV2data/processed_rnaseqV2/coaread/coadread_rsem_processed/'
colon <- read.csv(paste(path, 'colon_adenocarcinoma_raw.csv', sep=''), sep='')
colon_muc <- read.csv(paste(path, 'colon_mucinous_adenocarcinoma_raw.csv',sep=''),sep='')
rectal <- read.csv(paste(path, 'rectal_adenocarcinoma_raw.csv',sep=''),sep='')
data <- cbind(colon, colon_muc, rectal)
y <- c(rep(1, dim(colon)[2]), rep(2, dim(colon_muc)[2]), rep(3, dim(rectal)[2]))
prior <- c(dim(colon)[2], dim(colon_muc)[2], dim(rectal)[2])/length(y)
load('results/results_coadread_deplda.Rdata')
genes <- results_coadread_deplda$genes
n_top_acc <- results_coadread_deplda$n_top_acc
coadread <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(coadread, file='results/simpde_coadread.Rdata')
rm(colon)
rm(colon_muc)
rm(rectal)
rm(data)

print('KIPAN:')
path <- 'rnaseqV2data/processed_rnaseqV2/kipan/kipan_rsem_processed/'
clear <- read.csv(paste(path, 'kidney_clear_cell_renal_carcinoma_rsem_processed.csv', sep=''), sep='')
papillary <- read.csv(paste(path, 'kidney_papillary_renal_cell_carcinomarsem_processed_log2.csv', sep=''),sep='')
chromophobe <- read.csv(paste(path, 'kidney_chromophobe_rsem_processed.csv', sep=''), sep='')
data <- cbind(chromophobe, clear, papillary)
y <- c(rep(1, dim(chromophobe)[2]), rep(2, dim(clear)[2]), rep(3, dim(papillary)[2]))
prior <- c(dim(chromophobe)[2], dim(clear)[2], dim(papillary)[2])/length(y)
load('results/results_kipan_deplda.Rdata')
genes <- results_kipan_deplda$genes
n_top_acc <- results_kipan_deplda$n_top_acc
kipan <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(kipan, file='results/simpde_kipan.Rdata')
rm(clear)
rm(papillary)
rm(chromophobe)
rm(data)

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
load('results/results_lgg_deplda.Rdata')
genes <- results_lgg_deplda$genes
n_top_acc <- results_lgg_deplda$n_top_acc
lgg <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(lgg, file='results/simpde_lgg.Rdata')
rm(astro)
rm(oligo_dendroglioma)
rm(oligo_astro)
rm(data)


print('MESO')
path <- 'rnaseqV2data/processed_rnaseqV2/meso/meso_rsem_processed/'
epithelioid <- read.csv(paste(path, 'epithelioid_mesothelioma_rsem_processed.csv', sep=''), sep = '')
biphasic <- read.csv(paste(path, 'biphasic_mesothelioma_rsem_processed.csv', sep=''), sep='')
data <- cbind(biphasic, epithelioid)
## Define labels
y <- c(rep(1, dim(biphasic)[2]), rep(2, dim(epithelioid)[2]))
## Define prior
prior <- c(dim(biphasic)[2], dim(epithelioid)[2])/length(y)
load('results/results_meso_deplda.Rdata')
genes <- results_meso_deplda$genes
n_top_acc <- results_meso_deplda$n_top_acc
meso <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(meso, file='results/simpde_meso.Rdata')
rm(data)
rm(epithelioid)
rm(biphasic)

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
load('results/results_sarcoma_deplda.Rdata')
genes <- results_sarcoma_deplda$genes
n_top_acc <- results_sarcoma_deplda$n_top_acc
sarcoma <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(sarcoma, file='results/simpde_sarcoma.Rdata')
rm(leiomyosarcoma)
rm(liposarcoma)
rm(pleomorphic)
rm(myxofibrosarcoma)
rm(undiff_pleomorphic)
rm(nerve)
rm(data)

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
load('results/results_thca_deplda.Rdata')
genes <- results_thca_deplda$genes
n_top_acc <- results_thca_deplda$n_top_acc
thca <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(thca, file='results/simpde_thca.Rdata')
rm(classical)
rm(follicular)
rm(tall_cell)
rm(data)


print('THYMOMA')
path <- 'rnaseqV2data/processed_rnaseqV2/thymoma/thymoma_rsem_processed/'
ab <- read.csv(paste(path, 'thymoma_type_ab_rsem_processed.csv', sep=''), sep='')
b2 <- read.csv(paste(path, 'thymoma_type_b2_rsem_processed_log2.csv', sep=''), sep='')
a <- read.csv(paste(path, 'thymoma_type_a_rsem_processed.csv', sep=''), sep='')
b1 <- read.csv(paste(path, 'thymoma_type_b1_rsem_processed.csv', sep=''), sep='')
b3 <- read.csv(paste(path, 'thymoma_type_b3_rsem_processed_log2.csv', sep=''), sep='')
c <- read.csv(paste(path, 'thymoma_type_c_raw.csv', sep=''), sep='')
data <- cbind(a, ab, b1, b2, b3, c)
## Define labels
y <- c(rep(1, dim(a)[2]), rep(2, dim(ab)[2]), rep(3, dim(b1)[2]), rep(4, dim(b2)[2]), rep(5, dim(b3)[2]), rep(6, dim(c)[2]))
# Define prior
prior <- c(dim(a)[2], dim(ab)[2], dim(b1)[2], dim(b2)[2], dim(b3)[2], dim(c)[2])/length(y)
load('results/results_thymoma_deplda.Rdata')
genes <- results_thymoma_deplda$genes
n_top_acc <- results_thymoma_deplda$n_top_acc
thymoma <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(thymoma, file='results/simpde_thymoma.Rdata')
rm(a)
rm(ab)
rm(b1)
rm(b2)
rm(b3)
rm(c)
rm(data)

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
load('results/results_ucec_deplda.Rdata')
genes <- results_ucec_deplda$genes
n_top_acc <- results_ucec_deplda$n_top_acc
ucec <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(ucec, file='results/simpde_ucec.Rdata')
rm(endometrioid)
rm(serous)
rm(mixed)
rm(data)

print('UCS')
path <- 'rnaseqV2data/processed_rnaseqV2/ucs/ucs_rsem_processed/'
mmmt <- read.csv(paste(path, 'malignant_mixed_mullerian_tumor_rsem_processed.csv', sep=''), sep='')
mmmt <- mmmt[1:14975,]
mmmt_heterologous <- read.csv(paste(path, 'mmmt_heterologous_type_rsem_processed.csv', sep=''), sep='')
mmmt_homologous <- read.csv(paste(path, 'mmmt_homologous_type_rsem_processed.csv', sep=''), sep='')
data <- cbind(mmmt, mmmt_heterologous, mmmt_homologous)
## Define labels
y <- c(rep(1, dim(mmmt)[2]), rep(2, dim(mmmt_heterologous)[2]), rep(3, dim(mmmt_homologous)[2]))
## Define prior
prior <- c(dim(mmmt)[2], dim(mmmt_heterologous)[2], dim(mmmt_homologous)[2])/length(y)
load('results/results_ucs_deplda.Rdata')
genes <- results_ucs_deplda$genes
n_top_acc <- results_ucs_deplda$n_top_acc
ucs <- simplify.deplda(data, y, prior, 'mle', genes, n_top_acc)
save(ucs, file='results/simpde_ucs.Rdata')
rm(mmmt)
rm(mmmt_heterologous)
rm(mmmt_homologous)
rm(data)


print('DONE!!!')



