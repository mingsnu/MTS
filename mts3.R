## normal: A data frame or matrix for the normal data set
## abnormal: A data frame or matrix for the abnormal data set
## method: "im" instead for inverse matrix method and "gs" instead for Gram-Schmidt method

mts <- function(normal,abnormal,group=NULL,method="im"){
  if(NCOL(normal)!= NCOL(abnormal))
    stop("Normal data and Abnormal data should hhave the same number of columns")
  if(!is.matrix(normal))
    normal=as.matrix(normal)
  if(!is.matrix(abnormal))
    abnormal=as.matrix(abnormal)
  dimnames(normal) <- NULL
  dimnames(abnormal) <- NULL
  nvar = ncol(normal)
  mean.normal = colMeans(normal)
  sd.normal = apply(normal,2,sd)
  normal = scale(normal)
  abnormal = scale(abnormal,mean.normal,sd.normal)
  cov.normal = cov(normal)
  mean.normal = colMeans(normal)
  sd.normal = apply(normal,2,sd)
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
  if(is.null(group))
    md = list(md.normal=md.normal,md.abnormal=md.abnormal) else{
      level = levels(group <- as.factor(as.character(group)))
      n = length(level)
      mdabnormal = list()
      for(i in 1:n)
        mdabnormal[[i]] = md.abnormal[group==level[i]]
      names(mdabnormal) = level
      md = list(md.normal=md.normal,md.abnormal=mdabnormal)
    }
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

plot.mts <- function(score){
  plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))),ylab="mahalanobis distance")
}

mdoa <- function(dnormal,dabnormal,nvar=17,OA="L32.txt",...){
  Args <- match.call()
  OA <- read.table(paste("http://dreamhunter.me/guest/oa/",OA,sep=""))
  oa <- OA[,1:nvar]
  oa.nrow <- nrow(oa)
  md.oa <- list()
  for(i in 1:oa.nrow)
    md.oa[[i]] <- mts(dnormal[,oa[i,]==1],dabnormal[,oa[i,]==1],...)
  md.oa$call <- Args
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
ctable.default <- function(x,y,L=5){
  nn = sum(x <= L)
  an = sum(y <= L)
  judge = matrix(c(nn,an,length(x)-nn,length(y)-an),2)
  dimnames(judge) <- list(c("Normal","Abnormal"),c("Normal.test","Abnormal.test"))
  return(judge)
}

ctable.mdoa <- function(x,L=5){
  ## x: a list containing the returned results of mts function
  ## L: tolerance limit of the mahalanobis distance
  n = length(x)-1
  judge=array(NA,dim=c(2,2,n))
  dimnames(judge) <- list(c("Normal","Abnormal"),c("Normal.test","Abnormal.test"),paste("OA",1:n,sep=""))
  for(i in 1:n){
    nn = sum(x[[i]][[1]]<=L)
    an = sum(x[[i]][[2]]<=L)
    judge[,,i] = matrix(c(nn,an,length(x[[i]][[1]])-nn,length(x[[i]][[2]])-an),2)
  }
  return(judge)
}
ctable(md.oa)

ctable.mts <- function(x,L=5){
  nn = sum(x[[1]]<=L)
  an = sum(x[[2]]<=L)
  judge = matrix(c(nn,an,length(x[[1]])-nn,length(x[[2]])-an),2)
  dimnames(judge) <- list(c("Normal","Abnormal"),c("Normal.test","Abnormal.test"))
  return(judge)
}

####scoreselect
plot.mdoa <- function(md.oa,sleep=1){
  Args <- as.list(md.oa$call)
  aa=unlist(strsplit(Args$OA,split=NULL))
  bb=which(aa==".")-1
  cc=aa[2:bb]
  if(length(cc)==1)
    oa.nrow = as.numeric(cc) else{
      if(length(cc)>1)
        dd=cc[1]
      for(i in 2:length(cc))
        oa.nrow <- paste(dd,cc[i],sep="")
      oa.nrow <- as.numeric(oa.nrow)
    }
  for(i in 1:oa.nrow){
    score = md.oa[[i]]
    plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))),ylab="mahalanobis distance")
    Sys.sleep(sleep)
  }
}

as.list(md.oa$call)$OA


# calculate larger the better type S/N ratio for abnormal data
snchart <- function(md.oa){
  Args <- as.list(md.oa$call)
  nvar <- Args$nvar
  OA <- read.table(paste("http://dreamhunter.me/guest/oa/",Args$OA,sep=""))
  oa <- OA[,1:nvar]
  oa.nrow <- nrow(oa)
  sn <- c()
  for (i in 1:oa.nrow){
    sn[i] <- -10*log10(mean(1/(md.oa[[i]][[2]])^2))
  }
  
  snmean <- matrix(NA,nvar,2)
  for(i in 1:nvar){
    snmean[i,] <- aggregate(sn,by=list(oa[,i]),mean)[,2]
  }
  x0 <- seq(1,(2*nvar),by=2)
  x1 <- seq(2,(2*nvar),by=2)
  plot(0,xlim=c(1,2*nvar),ylim=range(snmean),type="n",xaxt="n",ylab="SN Ratio",xlab="Variables")
  axis(1,at=1:(2*nvar),labels=rep(1:2,nvar))
  segments(x0,snmean[,1],x1,snmean[,2])
  gain=snmean[,1]-snmean[,2]
  return(gain)
}

gain <- snchart(md.oa)
score = mts(dnormal[,gain>0],dabnormal[,gain>0])
plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))),main="MTS score after optimization")
ctable(score)


############ Optimum threshold
threshold <- function(x,...) UseMethod("threshold")
threshold.default <- function(x,y){
  m <- length(x)
  n <- length(y)
  x.max <- max(x)
  x.min <- min(x)
  y.max <- max(y)
  y.min <- min(y)
  if(x.max <= y.min)
    threshold <- (x.max + y.min)/2 else
  if(x.max < y.max && x.min < y.min){
    cross <- sort(unique(c(x[x >= y.min],y[y <= x.max],max(x[x < y.min]),min(y[y > x.max]))))
    l.cross <- length(cross)
    power <- rep(NA,l.cross)
    for(i in 1:l.cross){
      power[i] <- sum(x <= cross[i])/m + sum(y > cross[i])/n
    }
    index <- which(power==max(power))
    threshold <- (cross[index] + cross[index+1])/2
  }else
  if(x.max <= y.max &&  x.min >= y.min){
    if(sum(y < x.min) == 0 & sum(y > x.max) > 0)
      cross <- sort(unique(c(x,y[y <= x.max & y >= x.min],min(y[y > x.max])))) else
    if(sum(y < x.min) > 0 & sum(y > x.max) == 0)
      cross <- sort(unique(c(x,y[y <= x.max & y >= x.min],max(y[y < x.min]))))else
    if(sum(y < x.min) == 0 & sum(y > x.max) == 0)
      cross <- sort(unique(c(x,y[y <= x.max & y >= x.min]))) else
    cross <- sort(unique(c(x,y[y <= x.max & y >= x.min],max(y[y < x.min]),min(y[y > x.max]))))
    y.lo <- sum(y < x.min)
    y.hi <- sum(y > x.max)
    l.cross <- length(cross)
    power <- rep(NA,l.cross)
    if(y.lo <= y.hi){ 
      for(i in 1:l.cross){
        power[i] <- sum(x <= cross[i])/m + sum(y > cross[i])/n
      }
      index <- which(power==max(power))
      threshold <- (cross[index] + cross[index+1])/2 #index+1 must exist, because the highest power points can't appear at the end point
    }else{
      for(i in 1:l.cross){
        power[i] <- sum(x >= cross[i])/m + sum(y < cross[i])/n
      }
      index <- which(power==max(power))
      threshold <- (cross[index] + cross[index+1])/2
    }
  }else
  if(x.max > y.max &&  x.min < y.min){
    cross <- sort(unique(c(y,x[x <= y.min & x >= y.max],max(x[x < y.max]),min(x[x > y.max]))))
    x.lo <- sum(x < y.min)
    x.hi <- sum(x > y.max)
    l.cross <- length(cross)
    power <- rep(NA,l.cross)
    if(x.lo <= x.hi){ 
      for(i in 1:l.cross){
        power[i] <- sum(x >= cross[i])/m + sum(y < cross[i])/n
      }
      index <- which(power==max(power))
      threshold <- (cross[index] + cross[index+1])/2 #index+1 must exist, because the highest power points can't appear at the end point
    }else{
      for(i in 1:l.cross){
        power[i] <- sum(x <= cross[i])/m + sum(y > cross[i])/n
      }
      index <- which(power==max(power))
      threshold <- (cross[index] + cross[index+1])/2
    }
  } else
  if(x.min < y.max && x.min > y.min && x.max > y.max){
    cross <- sort(unique(c(x[x <= y.max],y[y >= x.min],min(x[x > y.max]),max(y[y < x.min]))))
    l.cross <- length(cross)
    power <- rep(NA,l.cross)
    for(i in 1:l.cross){
      power[i] <- sum(x >= cross[i])/m + sum(y < cross[i])/n
    }
    index <- which(power==max(power))
    threshold <- (cross[index] + cross[index+1])/2
  } else
  if(x.min >= y.max)
    threshold <- (x.min + y.max)/2
  return(list(threshold=threshold,power=max(power)))
}


  
threshold.mts <- function(x){
    normal <- x$md.normal
    abnormal <- x$md.abnormal
    thresh <- threshold.default(normal,abnormal)
  return(thresh)
}

print.mts <- function(score){
  class(score) <- NULL
  print(score)
}
print.mdoa <- function(md.oa){
  md.oa$call <- NULL
  class(md.oa) <- NULL
  print(md.oa)
}

  



###### liver disease example
data=read.table(file.choose(),header=T,sep=",")
data1=data[,c(-14,-16)]
data_normal=data1[data1[,14]=="Y",]
data_abnormal=data1[data1[,14]=="N",]
dnormal=data_normal[,-14]
dabnormal=data_abnormal[,-14]


score = mts(dnormal,dabnormal,"gs")
score
#score = mts(dnormal,dabnormal,"gs")
plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))))
title("MD score before optimization")

md.oa=mdoa(dnormal,dabnormal,13,"L16.txt","gs")
ctable(md.oa,L=2)


data1=data[,c(-13,-14,-16)]
data_normal=data1[data1[,13]=="Y",]
data_abnormal=data1[data1[,13]=="N",]
dnormal=data_normal[,-13]
dabnormal=data_abnormal[,-13]
score = mts(dnormal,dabnormal,"gs")
sort(score$md.abnormal)[round(length(score$md.abnormal)*0.05)]
md.oa=mdoa(dnormal,dabnormal,nvar=12,"L16.txt")
ctable(md.oa,1.119777)

OA = read.table("http://dreamhunter.me/guest/oa/L16.txt")
oa=oa=OA[,1:12]
score=mts(dnormal[,oa[9,]==1],dabnormal[,oa[9,]==1])
plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))))





data1=data[,c(-13,-15,-16)]
data_normal=data1[data1[,13]=="Y",]
data_abnormal=data1[data1[,13]=="N",]
dnormal=data_normal[,-13]
dabnormal=data_abnormal[,-13]
score = mts(dnormal,dabnormal,"gs")
sort(score$md.abnormal)[round(length(score$md.abnormal)*0.05)]
md.oa=mdoa(dnormal,dabnormal,nvar=12,"L16.txt")
ctable(md.oa,1.119777)

OA = read.table("http://dreamhunter.me/guest/oa/L16.txt")
oa=oa=OA[,1:12]
score=mts(dnormal[,oa[9,]==1],dabnormal[,oa[9,]==1])
plot(log(unlist(score)),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))))






#### 2.5~4 include 13column
data1=data[,-16]
data_normal=data1[data1[,13]=="Y"&data1[,15]=="Y",]
data_abnormal=data1[data1[,13]=="N"&data1[,14]=="Y",]
dnormal=data_normal[,c(-13,-14)]
dabnormal=data_abnormal[,c(-13,-14)]
score = mts(dnormal,dabnormal)
sort(score$md.abnormal)[round(length(score$md.abnormal)*0.05)]
md.oa=mdoa(dnormal,dabnormal,nvar=12,"L16.txt")
ctable(md.oa,1.119777)

OA = read.table("http://dreamhunter.me/guest/oa/L16.txt")
oa=oa=OA[,1:12]
score=mts(dnormal[,oa[9,]==1],dabnormal[,oa[9,]==1])
plot(log(unlist(score)),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))))



pnorm(-0.1)+pnorm(0.1)





#iris data
Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
head(Iris)
Iris1 <- Iris[,1:4]
iris.normal <- Iris1[Iris[,5]=="s",]
iris.abnormal1 <- Iris1[Iris[,5]=="c",]
iris.abnormal2 <- Iris1[Iris[,5]=="v",]
score1 <- mts(iris.normal,iris.abnormal1)
score2 <- mts(iris.normal,iris.abnormal2)
score1
score2
plot(c(score1$md.normal,score1$md.abnormal,score2$md.abnormal),pch=rep(1:3,each=50))
threshold(score1)
threshold(score2)
plot(score2)
threshold(score1$md.abnormal,score2$md.abnormal)
threshold.default(score1$md.abnormal,score2$md.abnormal)
abline(h=122)
ctable(score1$md.abnormal,score2$md.abnormal,L=122)
md.oa=mdoa(iris.normal,iris.abnormal1,nvar=4,OA="L8.txt")
md.oa
plot(md.oa,sleep=3)
ctable(md.oa,10)
gain <- snchart(md.oa)
score1 <- mts(iris.normal[,gain>0],iris.abnormal1[,gain>0])


iris.normal <- Iris1[Iris[,5]=="c",]
iris.abnormal1 <- Iris1[Iris[,5]=="s",]
iris.abnormal2 <- Iris1[Iris[,5]=="v",]
score1 <- mts(iris.normal,iris.abnormal1)
score2 <- mts(iris.normal,iris.abnormal2)
score1
score2
plot(c(score1$md.normal,score1$md.abnormal,score2$md.abnormal),pch=rep(1:3,each=50))

iris.normal <- Iris1[Iris[,5]=="v",]
iris.abnormal1 <- Iris1[Iris[,5]=="s",]
iris.abnormal2 <- Iris1[Iris[,5]=="c",]
score1 <- mts(iris.normal,iris.abnormal1)
score2 <- mts(iris.normal,iris.abnormal2)
score1
score2
plot(c(score1$md.normal,score1$md.abnormal,score2$md.abnormal),pch=rep(1:3,each=50))


md.oa=mdoa(iris.normal,iris.abnormal1,nvar=4,OA="L8.txt")
md.oa

ctable(md.oa,10)

for(i in 1:oa.nrow){
score = md.oa[[i]]
plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))))
Sys.sleep(1)
}

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


score11 <- mts(iris.normal[,gain>0],iris.abnormal1[,gain>0])
score22 <- mts(iris.normal[,gain>0],iris.abnormal2[,gain>0])
plot(c(score11$md.normal,score11$md.abnormal,score22$md.abnormal),pch=rep(1:3,each=50))
x11()
plot(c(score1$md.normal,score1$md.abnormal,score2$md.abnormal),pch=rep(1:3,each=50))

##
Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
head(Iris)
iris.normal <- Iris[Iris[,5]=="s",][,1:4]
iris.abnormal <- Iris[Iris[,5]!="s",][,1:4]
group <- Iris[Iris[,5]!="s",][,5]
score <- mts(iris.normal,iris.abnormal,group)
#L <- threshod()
plot.mts(score$md.abnormal)
abline(h=122)
ctable(score$md.abnormal$c,score$md.abnormal$v,L=122)




### spect.test
spect.t <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.test",sep=",")
head(spect.t)
normal <- spect.t[spect.t[,1]==1,][,-1]
abnormal <- spect.t[spect.t[,1]==0,][,-1]
score <- mts(normal,abnormal)
score <- mts(normal,abnormal,"gs")
score
plot(score)
ctable(score,L=1)

md.oa <- mdoa(normal,abnormal,nvar=22,OA="L32.txt")
md.oa
plot.mdoa(md.oa)
#L <- threshod(md.oa)
ctable(md.oa,L=0.5)
gain <- snchart(md.oa)
score1 = mts(normal[,gain>0],abnormal[,gain>0])
plot(score1)
ctable(score1,L=1)



spectf.t <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECTF.test",sep=",")
head(spectf.t)
normal <- spectf.t[spectf.t[,1]==1,][,-1]
abnormal <- spectf.t[spectf.t[,1]==0,][,-1]
nrow(normal)
nrow(abnormal)
score <- mts(normal,abnormal)
score
plot(score)

score


?lda

Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
z <- lda(Sp ~ ., Iris)
table(Iris$Sp, predict(z)$class)

plot(predict(z)$x, type="n",
     xlab="LD I", ylab="LD II",
     main="Untrained Iris LDA (n=150)")
text(predict(z)$x,
     levels(predict(z)$class)[predict(z)$class],
     col=unclass(Iris$Sp), cex=1.5)
abline(h=0);abline(v=0)
