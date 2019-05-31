# Zero-inflated Poisson Linear Discriminant Classifier
# This implementation is largely based on the ZIPLDA function made by Zhou Yan
# Own implementation for the purpose of understanding it better

source('helper_functions.R')
ziplda <- function(x,x_test, y, beta = 1, prior, method = c('mle','deseq','quantile'), prob0){
  # x : n x p matrix, where n indicates observations and p indicates features
  # x_test : m x p matrix for testing, where m indicates observations and p indicates features
  # y : class labels of each observation in x
  # beta : prior for Gamma(beta, beta) distribution
  # prior : prior probabilities for classes y. Vector of lenght unique(y)
  # method: method to compute size factor
  # return predicted classes and the log probabbilities for each observation
  
  Nhat <- compute.nhat(x)
  Nhat_test <- compute.nhat.test(x, x_test, method)
  classes = sort(unique(y))
  dhat <- compute.dhat(x, y, Nhat, beta)
  signx = sign(x_test==0)
  p = matrix(0,ncol = length(classes),nrow = nrow(x_test))
  for (i in 1:length(classes)){
    for (j in 1:nrow(x_test)){
      dstar = dhat[i,]
      part2=Nhat_test[j,]*dstar
      part1=prob0[j,]+(1-prob0[j,])*exp(-part2)
      part1[part1==0]=1
      p[j,i] <-sum(signx[j,]*log(part1))+sum(x_test[j,]*(1-signx[j,])*log(dstar))-sum((1-signx[j,])*part2)+log(prior[i])
    }
  }
  yhat = classes[apply(p,1,which.max)]
  return(list(p = p, yhat = yhat))
}
