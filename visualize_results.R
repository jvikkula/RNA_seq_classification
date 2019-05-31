# VISUALISATION OF THE DATA
tumor <- c('BRCA','BRCA',
           'COADREAD','COADREAD','COADREAD',
           'KIPAN','KIPAN','KIPAN',
           'LGG','LGG','LGG',
           'MESO','MESO',
           'SARCOMA','SARCOMA','SARCOMA','SARCOMA','SARCOMA','SARCOMA',
           'THCA','THCA','THCA',
           'THYM','THYM','THYM','THYM','THYM','THYM',
           'UCEC','UCEC','UCEC',
           'UCS','UCS','UCS')
#subtype <- c('ductal','lobular',
#             'colon','rectal','colon mucinous',
#             'clear cell','papillary','chromophobe',
#             'astrocytoma','oligodendroglioma', 'oligoastrocytoma',
#             'epithelioid','biphasic',
#             'leiomyosarcoma', 'dedifferentiated','pleomorphic', 'myxofibrosarcoma', 'undifferentiated','peripheral',
#             'classical','follicular','tall cell',
#             'ab','b2','a','b1','b3','c',
#             'endometrioid','serous','mixed',
#             'mmmt', 'heterologous','homologous')
subtype <- c('type1','type2',
             'type1','type2','type3',
             'type1','type2','type3',
             'type1','type2','type3',
             'type1','type2',
             'type1','type2','type3','type4','type5','type6',
             'type1','type2','type3',
             'type1','type2','type3','type4','type5','type6',
             'type1','type2','type3',
             'type1','type2','type3')
n_samples <- c(782, 203, 
               257, 89, 38,
               533, 290, 66, 
               194, 191, 130, 
               57, 23,
               104, 58, 29, 25, 21, 9, 
               358, 102, 36, 
               35, 31, 17, 15, 11, 11,
               113, 58, 12, 
               24, 20, 13)
data <- data.frame(tumor,subtype, n_samples)

ggplot(data, aes(x=as.factor(subtype), y=n_samples, fill=subtype)) + 
  geom_bar(stat='identity') + 
  facet_wrap(~as.factor(tumor),nrow=1) + 
  #theme(panel.background = element_rect(fill = "white")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


########################################

tumor <- c('BRCA','COADREAD','KIPAN','LGG','MESO','SARCOMA','THCA','THYM','UCEC','UCS')
hist_type <- c('type1','type2','type3','type4','type5','type6')
n_samples <- c(782, 203, 0, 0, 0, 0,
               257, 89, 38, 0, 0, 0,
               533, 290, 66, 0, 0, 0,
               194, 191, 130, 0, 0, 0,
               57, 23, 0, 0, 0, 0, 
               104, 58, 29, 25, 21, 9, 
               358, 102, 36, 0, 0, 0,
               35, 31, 17, 15, 11, 11,
               113, 58, 12, 0, 0, 0,
               24, 20, 13, 0, 0, 0)
data <- matrix(data = n_samples, nrow = 10, ncol = 6, byrow=TRUE)
rownames(data) <- tumor
colnames(data) <- hist_type
data
barplot(t(data), col=c('royalblue3','palegreen3','orchid3','purple3','coral2','darkblue'), 
        border="white", space=0.04, font.axis=2, xlab="Data set", ylab='Number of samples', ylim=c(0,1000))
#barplot(t(data), col=c('royalblue1','palegreen') , border="white", font.axis=2, beside=T, legend=colnames(data), xlab="group", font.lab=2)

#############################################
# VISUALIZATION OF THR DISTRIBUTION OF THE DATA
path <- '../rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/'
ductal <- read.csv(paste(path, 'brca_ductal_rsem_processed.csv', sep=''), sep='')
lobular <- read.csv(paste(path, 'brca_lobular_rsem_processed.csv', sep=''), sep='')
avg_ductal <- rowMeans(ductal)
avg_lobular <- rowMeans(lobular)
rm(ductal) 
rm(lobular)
df <- data.frame(subtype=factor(rep(c('ductal','lobular'),c(length(avg_ductal), length(avg_lobular)))),
                 counts = c(avg_ductal, avg_lobular))
head(df)

# Change colors
ggplot(df, aes(x=counts, color=subtype)) + 
  geom_histogram(fill="white", bins = 200) +
  xlim(c(0,20000)) + 
  ylim(c(0, 6000))


# Histogram with density plot
ggplot(df, aes(x=counts)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 100)+
  geom_density(alpha=.2, fill="#FF6666") +
  xlim(c(0,10000))


x <- rpois(1000, lambda = 0.005)
ggplot(x, aes(x=x)) +
  geom_histogram(fill='white')

########################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(readxl)
rows <- c('BRCA','COADREAD','KIPAN (binary)','KIPAN (all classes)','LGG','MESO','SARCOMA','THCA','THYM','UCEC','UCS')
processed <- read_excel('../data/accuracies.xlsx', sheet = 'processed data')
processed <- processed[-1]
processed <- as.matrix(processed)
rownames(processed) <- rows
keep <- c(1,2,4,5,6,7,8,9,10,11)
processed <- processed[keep,]
keep_cols <- c(1,3,5,6,7,8,9,10)
processed <- processed[,keep_cols]
rownames(processed) <- c('BRCA','COADREAD','KIPAN','LGG','MESO','SARCOMA','THCA','THYM','UCEC','UCS')

log2 <- read_excel('../data/accuracies.xlsx', sheet = 'log2 data')
log2 <- log2[-1]
log2 <- as.matrix(log2)
rownames(log2) <- rows

de <- read_excel('../data/accuracies.xlsx', sheet = 'DE data')
de <- de[-1]
de <- as.matrix(de)
rownames(de) <- rows
de <- de[keep,]
keep_cols <- c(1,2,3,4,7)
de <- de[,keep_cols]
rownames(de) <- c('BRCA','COADREAD','KIPAN','LGG','MESO','SARCOMA','THCA','THYM','UCEC','UCS')

raw <- read_excel('../data/accuracies.xlsx', sheet = 'raw data')
raw <- raw[-1]
raw <- as.matrix(raw)
rownames(raw) <- rows
raw <- raw[keep,]
#keep_cols <- c(1,2,3,4,6,7,8)
#raw <- raw[,keep_cols]
rownames(raw) <- c('BRCA','COADREAD','KIPAN','LGG','MESO','SARCOMA','THCA','THYM','UCEC','UCS')

deplda <- read_excel('../data/accuracies.xlsx', sheet = 'DE PLDA')
deplda <- deplda[-1]
deplda <- as.matrix(deplda)
rownames(deplda) <- c('BRCA','COADREAD','KIPAN','LGG','MESO','SARCOMA','THCA','THYM','UCEC','UCS')

#####
keep <- c(1,2,4,5,6,7,8,9,10,11)
df <- data.frame(PLDA=processed[keep,1], DE_PLDA=deplda[,1])
rownames(df) <- c('BRCA','COADREAD','KIPAN','LGG','MESO','SARCOMA','THCA','THYM','UCEC','UCS')
boxplot(df, ylab = 'Accuracy', xlab = 'Method', main='Processed data', col='skyblue2', font.axis=2)

################################################
###### plot box plot for presenation !!!!! #####
op <- par(mar = c(9,4,4,2) + 0.1)
par(mfrow=c(1,3))
par(cex.lab=1.5) # is for y-axis
par(cex.axis=1.3) # is for x-axis
boxplot(raw, ylab = 'Accuracy', main = 'Raw data', ylim = c(0.2,1),las=2, cex.main=2)
boxplot(processed, ylab = 'Accuracy', main = 'Processed data', ylim = c(0.2,1), las=2, cex.main=2)
boxplot(de, ylab = 'Accuracy', main = 'DE analysed data', ylim = c(0.2,1), las=2, cex.main=2)
par(op)

##### DEPLDA boxplot ####
plda <- c(0.7959, 0.487, 0.8459, 0.5456, 0.8, 0.6101, 0.5102, 0.5893, 0.7321, 0.4758)
deplda_sig_genes <- c(0.7817, 0.4342, 0.8652, 0.5146, 0.6875, 0.6735, 0.4388, 0.625, 0.7222, 0.18)
deplda_inc_acc <- c(0.797, 0.4868, 0.8596, 0.4951, 0.6875, 0.6531, 0.4592, 0.625, 0.75, 0.5)
deplda_max_acc <- c(0.8274, 0.4079, 0.8258, 0.534, 0.6875, 0.5714, 0.4388, 0.625, 0.6111, 0.4545)
df <- data.frame(plda, deplda_sig_genes, deplda_inc_acc, deplda_max_acc)
colnames(df) <- c('PLDA', 'DEPLDA\nsignificant genes', 'DEPLDA\nincrease accuracy', 'DEPLDA\nmax accuracy')
op <- par()
par(cex.lab=1.5) # is for y-axis
boxplot(df, ylab=c('Accuracy'))
par(op)
#########################
library(ggplot2)

# SUBTYPE COUNT
ggplot(df)+
  geom_boxplot(aes(x=dataset, y=accuracy,fill=subtype_count))+
  facet_grid(~subtype_count) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ORDERING BASED ON SAMPLE COUNT
o <- c(1,3,4,7,2,6,9,8,5,10)
processed_ordered <- processed[o,]
df <- data.frame(accuracy=as.vector(processed_ordered), method=method, dataset=dataset[o], subtype_count=subtype_count[o])
#df$dataset <- factor(df$dataset, levels = dataset[o])
ggplot(df) + 
  geom_boxplot(aes(x=dataset, y=accuracy, fill=subtype_count)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle('Processed data ordered based on sample count')

# ORDERIG BASED ON BALANCED DATA
o <- c(4,10,8,6,5,9,3,7,2,1)
processed_ordered <- processed[o,]
df <- data.frame(accuracy=as.vector(processed_ordered), method=method, dataset=dataset[o], subtype_count=subtype_count[o])
df$dataset <- factor(df$dataset, levels = dataset[o])
ggplot(df) + 
  geom_boxplot(aes(x=dataset, y=accuracy, fill=subtype_count)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#boxplot(t(processed_ordered))
###################################################

par(mfrow=c(2,2))
boxplot(processed, ylab = 'Accuracy', xlab = 'Method', main = 'Processed data', ylim = c(0,1))
#boxplot(log2, ylab = 'Accuracy', xlab = 'Method', main = 'Log2 data', ylim = c(0,1))
boxplot(de, ylab = 'Accuracy', xlab = 'Method', main = 'DE analysed data', ylim = c(0,1))
boxplot(raw, ylab = 'Accuracy', xlab = 'Method', main = 'Raw data', ylim = c(0,1))

op <- par(mar = c(9,4,2,2) + 0.1)
boxplot(t(processed), col = 'skyblue3', las = 2, ylab =  'Accuracy', ylim = c(0.3,1))
par(op)

avg <- rowMeans(processed, na.rm = TRUE)
op <- par(mar = c(9,4,4,2) + 0.1)
barplot(avg, ylim=c(0,1), col = 'skyblue3',ylab = 'Average accuracy',
        border="white", space=0.07, font.axis=1, las=2, main = 'Average accuracy per data set')
par(op)

datasets <- c('BRCA','MESO')
class2 <- processed[which(rownames(processed)%in%datasets),]

datasets <- c('COADREAD','KIPAN','LGG','THCA','UCEC','UCS')
class3 <- processed[which(rownames(processed)%in%datasets),]
class3_avg <- rowMeans(class3)
o <- order(class3_avg, decreasing = TRUE)
class3 <- class3[o,]

datasets <- c('SARCOMA','THYM')
class6 <- processed[which(rownames(processed)%in%datasets),]

par(mfrow=c(1,3), mar = c(9,4,4,2) + 0.1, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
boxplot(t(class2), ylim=c(0.2,1), ylim = 'Accuracy', main='2 classes', las=2)
text(x=1, y=782/(782+203), 'X', col='blue',cex=1.5)
text(x=2, y=57/(57+23), 'X', col='blue', cex=1.5)
boxplot(t(class3), ylim=c(0.2,1), main='3 classes', las=2)
text(x=1, y=533/(533+66+290), 'X', col='blue', cex=1.5)
text(x=2, y=113/(113+12+58), 'X', col='blue', cex=1.5)
text(x=3, y=358/(358+102+36), 'X', col='blue', cex=1.5)
text(x=4, y=257/(257+38+89), 'X', col='blue', cex=1.5)
text(x=5, y=194/(194+191+130), 'X', col='blue', cex=1.5)
text(x=6, y=24/(24+20+13), 'X', col='blue', cex=1.5)
boxplot(t(class6), ylim=c(0.2,1), main='6 classes', las=2)
text(x=1, y=104/(104+58+29+25+21+9), 'X', col='blue', cex=1.5)
text(x=2, y=35/(17+35+15+31+11+11), 'X', col='blue', cex=1.5)

#########################
datasets <- c('BRCA','MESO')
class2 <- processed[which(rownames(processed)%in%datasets),]

datasets <- c('COADREAD','KIPAN','LGG','THCA','UCEC','UCS')
class3 <- processed[which(rownames(processed)%in%datasets),]

datasets <- c('SARCOMA','THYM')
class6 <- processed[which(rownames(processed)%in%datasets),]

op <- par(mfrow=c(1,3), mar = c(11,4,4,2) + 0.1, cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
boxplot(class2, ylim=c(0.2,0.95), main='BRCA & MESO', ylab = 'Accuracy', las=2)
boxplot(class3, ylim=c(0.2,0.95), main='COADREAD, KIPAN, LGG\nTHCA, UCEC & UCS', las=2)
boxplot(class6, ylim=c(0.2,0.95), main='SARC & THYM', las=2)
par(op)

#########################

path <- '../rnaseqV2data/processed_rnaseqV2/brca/brca_rsem_processed/'
ductal <- read.csv(paste(path, 'brca_ductal_rsem_processed.csv', sep=''), sep='')
lobular <- read.csv(paste(path, 'brca_lobular_rsem_processed.csv', sep=''), sep='')

genes <- rownames(ductal)
ductal <- as.matrix(ductal)
lobular <- as.matrix(lobular)
df <- data.frame(c(ductal[which(genes=='BRCA1'),], lobular[which(genes=='BRCA1'),]), 
                 c(ductal[which(genes=='BRCA2'),], lobular[which(genes=='BRCA2'),]),
                 c(ductal[which(genes=='CDH1'),], lobular[which(genes=='CDH1'),]),
                 c(ductal[which(genes=='TP53'),], lobular[which(genes=='TP53'),]),
                 c(rep('ductal',length(ductal[which(genes=='BRCA1'),])),rep('lobular',length(lobular[which(genes=='BRCA1'),]))))
colnames(df) <- c('BRCA1','BRCA2','CDH1','TP53','subtype')

par(mfrow=c(1,4))
boxplot(BRCA1 ~ subtype, df, main = 'BRCA1')
boxplot(BRCA2 ~ subtype, df, main = 'BRCA2')
boxplot(CDH1 ~ subtype, df, main = 'CDH1')
boxplot(TP53 ~ subtype, df, main = 'TP53')

rm(ductal)
rm(lobular)
