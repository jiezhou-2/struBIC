##For a given network matrix, compute
##the mle of precision matrix 
mle=function(data,priori){
  priori=priori+t(priori)
  priori=ifelse(priori==0,0,1)
  diag(priori)=0
  p=dim(data)[2]
  n=dim(data)[1]
  precision=matrix(0,nrow = p,ncol = p)
  ##for the first row
  for (j in 1:p) {
    y=data[,j]
    index=which(priori[j,]==1)
    if (length(index)==0) {
      precision[j,]=0
      precision[j,j]=1/var(y)
      next
    }
    x=data[,index]
    result=lm(y~0+x)
    alpha=result$coefficients
    sigma=t(result$residuals)%*%(result$residuals)/result$df.residual
    precision[j,index]=-alpha/sigma[1]
    precision[j,j]=1/sigma
  }
  precision=(precision+t(precision))/2
  return(precision)
}
