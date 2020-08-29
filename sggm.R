
sggm=function(data,lambda,M,prob){
p=dim(data)[2]
n=dim(data)[1]
m=length(lambda)
t=log(p)/((log(1/prob-1)))
web=modelset(data = data, lambda = lambda, P=M)
web=web[(sapply(web, sum)-p)/2<n]
sbic=seq(0,length=m)
for (i in 1:m) {
  ggbic=gbic(data=data, theta= web[[i]], prob = prob, P=M)
  sbic[i]=ggbic
}
mm=which.min(sbic)
wweb=web[[mm]]

return(list(networkhat=wweb, candidates=web))
}


