library(limma)

# This implementation of voom-transformation is heavily influenced by the 
# R package voom

my_voom <- function(counts, plot = FALSE){
  Ri <- colSums(counts)
  z <- t(log2(t(counts + 0.5)/(Ri + 1) * 1e+06))
  
  lin <- lmFit(z)
  sx <- lin$Amean + mean(log2(Ri + 1)) - 6*log2(10)
  sy <- sqrt(lin$sigma)
  l <- lowess(sx, sy, f = 0.5)
  
  if (plot == TRUE){
    plot(sx, sy, xlab = "log_2(counts)", ylab = "sqrt(stdev)", 
         pch = 16, cex = 0.25)
    lines(l, col = "red")
  }
  
  f <- approxfun(l, rule = 2)
  fvalues <- lin$coef %*% t(lin$design)
  fcpm <- 2^fvalues
  fcount <- 1e-06 * t(t(fcpm) * (Ri + 1))
  flogcount <- log2(fcount)
  w <- 1/f(flogcount)^4
  dim(w) <- dim(flogcount)
  
  ret <- list()
  ret$weights <- w
  ret$logvals <- z
  return(ret)
}

multiclass_preprocess <- function(sets, names, subset = c()){
  
  if (length(subset) > 0){
    for (i in 1:length(sets)){
      sets[[i]] <- sets[[i]][subset,]
    }
  }
  
  data <- do.call(cbind, sets)
  cols <- c()
  for (i in 1:length(sets)){
    cols <- c(cols, colnames(sets[1]))
  }
  
  classes <- matrix(nrow = ncol(data), ncol = 1)
  rownames(classes) <- colnames(data)
  colnames(classes) <- "Class"
  for (i in 1:length(names)){
    set <- which(rownames(classes) %in% colnames(as.data.frame(sets[i])))
    classes[set] <- names[i]
  }
  
  ret <- list()
  ret$sets <- sets
  ret$data <- data
  ret$classes <- classes
  
  return(ret)
  
}

do_folds <- function(sets, data){
  
  
  test <- FALSE
  
  n <- ncol(data) %/% 5
  
  while (test == FALSE){
    data <- data[,sample(ncol(data))]
    folds <- list()
    folds$f1 <- data[,1:n]
    folds$f2 <- data[,(n+1):(2*n+1)]
    folds$f3 <- data[,(2*n+2):(3*n+2)]
    folds$f4 <- data[,(3*n+3):(4*n+3)]
    folds$f5 <- data[,(4*n+4):ncol(data)]
    
    smpl <- TRUE
    for (u in 1:5){
      for (i in 1:length(sets)){
        index <- which(colnames(folds[[u]]) %in% colnames(sets[[i]]))
        if (length(index) < 1){
          smpl <- FALSE
        }
      }
    }
    
    if (smpl == TRUE){
      test <- TRUE
    }
    
  }
  
  
  return(folds)
}

multiclass_train <- function(sets, data, folds){ # Training the voomDQDA
  
  if (k == 1){
    testdata <- folds$f1
    traindata <- cbind(folds$f2, folds$f3, folds$f4, folds$f5)
  }
  
  else if (k == 2){
    testdata <- folds$f2
    traindata <- cbind(folds$f1, folds$f5, folds$f3, folds$f4)
  }
  
  else if (k == 3){
    testdata <- folds$f3
    traindata <- cbind(folds$f1, folds$f2, folds$f5, folds$f4)
  }
  
  else if (k == 4){
    testdata <- folds$f4
    traindata <- cbind(folds$f1, folds$f2, folds$f5, folds$f3)
  }
  
  else if (k == 5){
    testdata <- folds$f5
    traindata <- cbind(folds$f1, folds$f2, folds$f4, folds$f3)
  }
  
  trainsets <- list()
  for (i in 1:length(sets)){
    index <- which(colnames(traindata) %in% colnames(sets[[i]]))
    trainsets[[i]] <- traindata[,index]
  }
  
  
  vooms <- list()
  for (i in 1:length(sets)){
    vooms[[i]] <- my_voom(trainsets[[i]])
  }
  
  ss <- list()
  zs <- list()
  for (i in 1:length(sets)){
    s2w <- rowSums(vooms[[i]]$weights) / (rowSums(vooms[[i]]$weights)^2 - rowSums(vooms[[i]]$weights^2))
    zwk <- rowSums(vooms[[i]]$weights * vooms[[i]]$logvals) / rowSums(vooms[[i]]$weights)
    s2w <- s2w * rowSums(vooms[[i]]$weights * (vooms[[i]]$logvals - zwk)^2)
    ss[[i]] <- s2w
    zs[[i]] <- zwk
  }
  
  ret = list()
  ret$ss <- ss
  ret$zs <- zs
  ret$trainsets <- trainsets
  ret$testdata <- testdata
  
  return(ret)
  
}

multiclass_predict <- function(testdata, sets, zs, ss, priors, classes){ 
  # Classifying using the voomDQDA
  testvoom <- my_voom(testdata)
  
  preds <- matrix(nrow = 2, ncol = ncol(testvoom$logvals))
  colnames(preds) <- colnames(testdata)
  
  for (n in 1:ncol(testvoom$logvals)){
    is <- c()
    for (i in 1:length(sets)){
      is <- c(is, (-sum(((testvoom$logvals[,n] - zs[[i]])^2)/ss[[i]]) + 2 * log(priors[i])))
    }
    for (i in 1:length(is)){
      
      if (is[i] == max(is)){
        preds[1,n] <- classes[i]
      }
      
    }
    
  }
  
  for (n in 1:ncol(testvoom$logvals)){
    for (i in 1:length(sets)){
      if (colnames(preds)[n] %in% colnames(sets[[i]])){
        preds[2,n] <- classes[i]
      }
    }
  }
  
  correct <- 0
  
  for (i in 1:ncol(testvoom$logvals)){
    if (preds[1,i] == preds[2,i]){
      correct <- correct + 1
    }
  }
  
  ret <- list()
  ret$cor <- correct
  ret$preds <- preds
  
  return(ret)
  
}
