

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

# using orthogonal array
md.oa=mdoa(dnormal,dabnormal)

# make the contingency table
ctable(md.oa)

# calculate larger the better type S/N ratio for abnormal data

gain <- snchart(md.oa)
score = mts(dnormal[,gain>0],dabnormal[,gain>0])
plot(unlist(score),pch=c(rep(21,length(score[[1]])),rep(24,length(score[[2]]))),main="MTS score after optimization")
ctable(score)


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
