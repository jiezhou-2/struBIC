comparison=function(real, estimate){
  ## M1 is the true net matrix
  ##Me is the estimated net matrix
  real=real+t(real)
  diag(real)=1
  estimate=estimate+t(estimate)
  diag(estimate)=1
  N1=ifelse(real==0,0,1)
  N2=ifelse(estimate==0,0,1)
   if (any(dim(N1)!=dim(N2))[1])
  stop("Two matrixes should have the same dimension")
  p=dim(real)[1]
  real=(sum(N1)-p)/2
  select=(sum(N2)-p)/2
  real_select=(sum(N2[N1==1])-p)/2
  fause_select=sum(N2[N1==0])/2
  if (real==0){return( list(PSR=1, FDR=0))}
  if (select==0) {return(list(PSR=0, FDR=0))}
  PSR=real_select/real
  FDR=fause_select/select
    return(list(PSR=PSR, FDR=FDR)) 
}
