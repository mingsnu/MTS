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
