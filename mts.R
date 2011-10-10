## normal: A data frame or matrix for the normal data set
## abnormal: A data frame or matrix for the abnormal data set
## method: "im" instead for inverse matrix method and "gs" instead for Gram-Schmidt method

mts <- function(normal,abnormal,method="im"){
  if(ncol(normal)!=ncol(abnormal))
    stop("Normal data and Abnormal data should hhave the same number of columns")
  if(!is.matrix(normal))
    normal=as.matrix(normal)
  if(!is.matrix(abnormal))
    abnormal=as.matrix(abnormal)
  dimnames(normal) <- NULL
  dimnames(abnormal) <- NULL
  nvar = ncol(normal)
  mean.normal = colMeans(normal)
  cov.normal = cov(normal)
  if(method == "im"){                        #inverse matrix method
    md.normal = mahalanobis(normal,mean.normal,cov.normal)/nvar #mahalanobis distance for normal data
    md.abnormal = mahalanobis(abnormal,mean.normal,cov.normal)/nvar #mahalanobis distance for abnormal data
  }else
  if(method == "gs"){                        #gram-schmidt method
  ##  U = normal.center = sweep(normal,2,mean.normal)
  ## for(j in 2:nvar)
  ##   for(i in 1:(j-1))            #gram-schmidt general method
  ##     U[,j] = U[,j] - sum(normal.center[,j]*U[,i])/sum(U[,i]*U[,i])*U[,i]
  ##   s = apply(U,2,var)
  ## md.normal = colSums(t(U*U)/s)/nvar
  
    U = normal.center = sweep(normal,2,mean.normal)
    for(j in 2:nvar)
      for(i in 1:(j-1))             #gram-schmidt numerical stability method
        U[,j] = U[,j] - sum(U[,i]*U[,j])/sum(U[,i]*U[,i])*U[,i]
    s = apply(U,2,var)
    md.normal = colSums(t(U*U)/s)/nvar

    Uab = sweep(abnormal,2,mean.normal)
    for(j in 2:nvar)
      for(i in 1:(j-1))              #gram-schmidt general method
        Uab[,j] = Uab[,j] - sum(normal.center[,j]*U[,i])/sum(U[,i]*U[,i])*Uab[,i]
    md.abnormal = colSums(t(Uab*Uab)/s)/nvar
  }else
  cat("There is no such method!")
  md = list(md.normal=md.normal,md.abnormal=md.abnormal)
  class(md) <- "mts"
  return(md)
}



###### liver disease example
data=read.table(file.choose(),header=T,sep=",")
data1=data[,-1]
data_normal=data1[data1$diagnosis==1,]
data_abnormal=data1[data1$diagnosis!=1,]
dnormal=data_normal[,-18]
dabnormal=data_abnormal[,-18]


score = mts(dnormal,dabnormal)
score
#score = mts(dnormal,dabnormal,"gs")
plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))))
title("MD score before optimization")


mdoa <- function(dnormal,dabnormal,nvar=17,OA="L32.txt"){
  OA = read.table(paste("http://dreamhunter.me/guest/oa/",OA,sep=""))
  oa=OA[,1:nvar]
  oa.nrow=nrow(oa)
  md.oa = list()
  for(i in 1:oa.nrow)
    md.oa[[i]] = mts(dnormal[,oa[i,]==1],dabnormal[,oa[i,]==1])
  class(md.oa) <- "mdoa"
  md.oa
}
md.oa=mdoa(dnormal,dabnormal)

## OA = read.table("http://dreamhunter.me/guest/oa/L32.txt")
## nvar=17
## oa=OA[,1:nvar]
## oa.nrow=nrow(oa)
## md.oa = list()
## for(i in 1:oa.nrow)
##   md.oa[[i]] = mts(dnormal[,oa[i,]==1],dabnormal[,oa[i,]==1])


## make the contingency table
ctable <- function(x,...) UseMethod("ctable")
ctable.mdoa <- function(x,L=5){
  ## x: a list containing the returned results of mts function
  ## L: tolerance limit of the mahalanobis distance
  n = length(x)
  judge=array(NA,dim=c(2,2,n))
  dimnames(judge) <- list(c("Normal","Abnormal"),c("Normal.test","Abnormal.test"),paste("OA",1:n,sep=""))
  for(i in 1:n){
    nn = sum(x[[i]][[1]]<=5)
    an = sum(x[[i]][[2]]<=5)
    judge[,,i] = matrix(c(nn,an,length(x[[i]][[1]])-nn,length(x[[i]][[2]])-an),2)
  }
  return(judge)
}
ctable(md.oa)

ctable.mts <- function(x,L=5){
  nn = sum(x[[1]]<=5)
  an = sum(x[[2]]<=5)
  judge = matrix(c(nn,an,length(x[[1]])-nn,length(x[[2]])-an),2)
  dimnames(judge) <- list(c("Normal","Abnormal"),c("Normal.test","Abnormal.test"))
  return(judge)
}

####scoreselect
for(i in 1:oa.nrow){
score = md.oa[[i]]
plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))))
Sys.sleep(1)
}


# calculate larger the better type S/N ratio for abnormal data
sn=c()
for (i in 1:oa.nrow){
sn[i]=-10*log10(mean(1/(md.oa[[i]][[2]])^2))
}

snmean=matrix(NA,nvar,2)
for(i in 1:nvar){
snmean[i,]=aggregate(sn,by=list(oa[,i]),mean)[,2]
}
gain=snmean[,1]-snmean[,2]
plot(gain)
abline(h=0)

score = mts(dnormal[,gain>0],dabnormal[,gain>0])
plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))),main="MTS score after optimization")
ctable(score)
