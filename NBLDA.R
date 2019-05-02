NBLDA <- function(x, x_test, y, disperhat, beta = 1, prior = NULL, method = c('mle', 'deseq', 'quantile')) 
{
  # x: Matrix of size n x p where n is amount of observations and p is amount of features
  # x_test: Matrix of size m x p for predicting class labels 
  # y: Labels for each observation in x
  # disperhat: Dispersion parameter for NBLDA model
  # beta: Prior for Gamma(beta, beta) distribution
  # method: method used to compute the size factor
  
  # Set prior for all classes if not set before
  if (is.null(prior)){
    prior <- rep(1/length(unique(y)), length(unique(y))) 
  }
  # Colsums, amount of reads per gene
  lambda_g <- colSums(x)
  # Rowsums, amount of reads per sample
  x_j <- rowSums(x)
  
  # Number of classes in the data
  classes <- sort(unique(y))
  n_classes <- length(classes)
  
  shats <- size_estimation(x, x_test, method)
  
  # Estimate dkgs
  # DKGs are fine
  dkg <- matrix(1, nrow = n_classes, ncol = ncol(x))
  for (i in 1:n_classes){
    num <- colSums(x[y == classes[i],]) + beta
    denom <- sum(shats$train[y == classes[i]]) * lambda_g + beta
    div <- num/denom
    colnames(div) = NULL
    dkg[i,] = div
  }
  
  disc <- matrix(NA, nrow = nrow(x_test), ncol = n_classes)
  
  # Loop through number of classes
  for (i in 1:n_classes) {
    # Loop through samples in test data
    for (j in 1:nrow(x_test))
    {
      dstar <- dkg[i,]
      p2 <- 1 + shats$test[j] * lambda_g * dstar * disperhat 
      p1 <- dstar / p2
      disc[j,i] <- sum(x_test[j,] * log(p1)) - sum( (1/disperhat) * log(p2) + log(prior[i]))
    }
  }
  y_test <- classes[apply(disc, 1, which.max)]
  return(list(pred = y_test, ns = shats, ds = dkg, disc = disc))
}

size_estimation <- function(x_train, x_test, method) {
  if (method == 'mle'){
    shat_train <- rowSums(x_train) / sum(x_train)
    shat_test <- rowSums(x_test) / sum(x_train)
  }
  else if (method == 'deseq') {
    train_counts <- t(x_train)
    test_counts <- t(x_test)
    geom_train <- exp(rowMeans(log(train_counts)))
    geom_test <- exp(rowMeans(log(test_counts)))
    rawsize_train <- apply(train_counts, 2, function(cnts) median((cnts/geom_train)[geom_train > 0]))
    rawsize_test <- apply(test_counts, 2, function(cnts) median((cnts/geom_test)[geom_test > 0]))
    shat_train <- apply(train_counts, 2, function(cnts) median((cnts/geom_train)[geom_train > 0], 
                                                               na.rm = TRUE))/sum(rawsize_train)
    shat_test <- apply(test_counts, 2, function(cnts) median((cnts/geom_train)[geom_test > 0], 
                                                             na.rm = TRUE))/sum(rawsize_test)
  }
  else if (method == 'quantile') {
    train_topq <- pmax(1, apply(x_train, 1, quantile, 0.75))
    test_topq <- pmax(1, apply(x_test, 1, quantile, 0.75))
    shat_train <- train_topq / sum(apply(x_train, 1, quantile, 0.75))
    shat_test <- test_topq / sum(apply(x_train, 1, quantile, 0.75))
  }
  return(list("train" = shat_train, "test" = shat_test))
}

estimate_dispersion <- function(X) {
  X = t(X)
  s_means <- rowMeans(X)
  s_vars <- rowVars(X)
  tt = getT(X, sizeFactors=rep(1,ncol(X)), plotASD = T)$target  
  print(tt)
  
  disp <- (s_vars - s_means) / s_means^2
  disp0 <- numeric()
  for (i in 1:length(disp)){
    d <- disp[i]
    disp0[i] <- max(0,d)
  }
  
  adj_disp <- getAdjustDisp(disp0, shrinkTarget = tt)$adj
  
  return(adj_disp)
  #delta <- calculate_weights(X,disp0)
  #adj_disp <- delta * epsilon + (1 - delta) * disp0
}

calculate_weights <- function(d0) {
  G <- length(d0)
  upbound <- quantile(d0, na.rm = TRUE)
  epsilon <- NULL
  
  delta <- (sum(d0 - (1/G) * sum(d0))^2 / (G-1)) / (sum(d0-epsilon)^2 / (G-2))
  return(delta)
}
