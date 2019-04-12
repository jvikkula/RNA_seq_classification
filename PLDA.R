# Poisson Linear Discriminant Analysis classifier
# This implementation is largely based on the PoiClaClu package
# Own implementation for the purpose of understanding it better

source('helper_functions.R')
PLDA <- function(x,x_test, y, beta = 1, prior=NULL, method = c('mle','deseq','quantile')){
  # x : n x p matrix, where n indicates observations and p indicates features
  # x_test : m x p matrix for testing, where m indicates observations and p indicates features
  # y : class labels of each observation in x
  # beta : prior for Gamma(beta, beta) distribution
  # prior : prior probabilities for classes y. Vector of lenght unique(y)
  # method: method to compute size factor
  # return predicted classes and the log probabbilities for each observation
  
  Xi_dot = rowSums(x)
  Xdot_j = colSums(x)
  Xdotdot = sum(x)
  
  Nhat <- compute.nhat(x)
  Nhat_test <- compute.nhat.test(x, x_test, method)
  classes = sort(unique(y))
  dhat <- compute.dhat(x, y, Nhat, beta)

  p = matrix(0,ncol = length(classes),nrow = nrow(x_test))
  for (i in 1:length(classes)){
    p[,i] = rowSums(scale(x_test, center=FALSE, scale=(1/log(dhat[i,])))) - rowSums(scale(Nhat_test, center = FALSE, scale = (1/dhat[i,]))) + log(prior[i])
  }
  yhat = classes[apply(p,1,which.max)]
  return(list(p = p, yhat = yhat)) 
}
