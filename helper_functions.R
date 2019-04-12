# Compute shat for PLDA and ZIPLDA
compute.shat <- function(x_train,x_test,method = c('mle','deseq','quatile')){
  # Compute size factor based on defined method
  # 'mle' : total count for the ith observation, which is based on the MLE for Nij 
  # 'deseq' : median ratio
  # 'quantile' 
  
  if (method == 'mle'){
    shat = rowSums(x_test)/sum(x_train)
  }
  else if (method == 'quantile'){
    shat = pmax(1, apply(x_test, 1, quantile, 0.75))
    shat = shat/sum(apply(x_train, 1, quantile, 0.75))
  }  
  else if (method == 'deseq'){
    counts_train <- t(x_train)
    counts_test <- t(x_test)
    geomeans <- exp(rowMeans(log(counts_train)))
    rawsizestr <- apply(counts_train, 2, function(cnts) median((cnts/geomeans)[geomeans > 0]))
    shat <- apply(counts_test, 2, function(cnts) median((cnts/geomeans)[geomeans > 0], na.rm = TRUE))/sum(rawsizestr)
  }
  return(shat)
}


compute.nhat <- function(x){
  Nhat = outer(rowSums(x), colSums(x), '*')/sum(x)
  return(Nhat)
}


compute.nhat.test <- function(x_train, x_test, method){
  shat = compute.shat(x_train, x_test, method)
  ghat = colSums(x_train)
  Nhat_test = outer(shat, ghat, '*')
  return(Nhat_test)
}


compute.dhat <- function(x,y,Nhat,beta=1){
  classes = sort(unique(y))
  dhat <- matrix(0,nrow=length(classes),ncol=length(colSums(x)))
  for(i in 1:length(classes)){
    numerator = colSums(x[y==classes[i],])+beta
    denominator = colSums(Nhat[y==classes[i],])+beta
    division = numerator/denominator
    colnames(division) = NULL
    dhat[i,] = division
  }
  return(dhat)
}

# Compute the probability of point mass at zero for ZIPLDA
# largely based on the function estimatep
point.mass.at.zero <- function(x_train, x_test, y, beta = 1, method = c('mle','deseq','quantile')){
  Nhat <- compute.nhat(x_train)
  Nhat_test <- compute.nhat.test(x_train, x_test, method)
  classes <- sort(unique(y))
  signx <- sign(x_test==0)
  dhat <- compute.dhat(x_train, classes, Nhat, beta)
  
  mu <- matrix(NA, nrow=nrow(x_test),ncol=ncol(x_test))
  for (i in 1:nrow(x_test)){
    dstar = dhat[y[i],]
    mu[i,] = Nhat_test[i,]*dstar    
  }
  
  #G = ncol(x_train)
  G = length(x_test[1,])
  #lib <- colSums(x_train)
  lib <- rowSums(x_test)
  #x1 <- t(x_train)
  x1 <- t(x_test)
  
  mu1 <- as.vector(t(mu))
  librep <- rep(lib,rep(G,length(lib)))/(lib[1])
  x2<-as.vector(x1)
  
  y<-x2
  y[y!=0]<-1
  xreg<-cbind(y,librep,mu1)
  glm.out<-glm(y~ librep+ mu1,family=binomial("logit"),data=data.frame(xreg))
  summary(glm.out)
  coef<-as.matrix(glm.out$coefficients)
  inter<-rep(1,G)
  muu<-t(mu)
  
  x_test1=t(x_test)
  p<-x_test1
  for(i in 1:length(x_test1[1,])){
    libsize<-rep(sum(x_test1[,i]),G)/lib[1]
    estx1<-cbind(inter,libsize,muu[,1])
    dd<-estx1%*% coef
    dd[dd>50]<-50
    dd[dd<(-50)]<--50
    p1<-exp(dd)/(1+exp(dd))
    p[,i]<-1-p1
  }
  return(t(p)) 
}
