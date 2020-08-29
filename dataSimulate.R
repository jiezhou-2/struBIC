###generate the network matrix
library(MASS)
simulate=function(n,p,m1,m2){
  real_stru=matrix(0, nrow = p, ncol = p)
  real_stru[lower.tri(real_stru,diag = T)]=1
  index=which(real_stru==0,arr.ind = T)
  a=sample(1:nrow(index),m1, replace = F)
  real_stru[index[a,]]=1
  real_stru[lower.tri(real_stru,diag = T)]=0
  real_stru=real_stru+t(real_stru)+diag(p)
  distrubance=real_stru
  distrubance[lower.tri(distrubance,diag = T)]=0
  index1=which(distrubance==1,arr.ind = T)
  distrubance=(distrubance+1)%%2
  distrubance[lower.tri(distrubance,diag = T)]=0
  index2=which(distrubance==1,arr.ind = T)
  a=sample(1:nrow(index1),m2, replace = F)
  b=sample(1:nrow(index2),m2, replace = F)
  distrubance=matrix(0,nrow = p,ncol = p)
  distrubance[index1[a,]]=1
  distrubance[index2[b,]]=1
  distrubance=distrubance+t(distrubance)
  prior_stru=(real_stru+distrubance)%%2
  diag(prior_stru)=1
  
  
  # generate the symmetrical precision matrix
  theta = matrix(rnorm(p^2), p)
  theta[lower.tri(theta, diag = TRUE)] = 0
  theta = theta + t(theta) + diag(p)
  
  # apply the reqired sparsity
  theta = theta * real_stru
  
  # force it to be positive definite
  eig=eigen(theta)$values
  theta = theta - (min(eig)-0.1) * diag(p)
  Sigma=solve(theta)
  mu=rep(0,p)
  data=MASS::mvrnorm(n,mu=mu,Sigma = Sigma)
  return(list(data=data, realnetwork=real_stru, priornetwork=prior_stru))
}
