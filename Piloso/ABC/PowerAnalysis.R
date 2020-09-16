library(abc)
#setwd("/Users/manolo/Documents/Novas_spDelPiloso2020/ABC")
sust<-read.table("SuSt.txt")
numsim=30000
numloc=15
sust <- data.frame(s1=tapply(sust[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s2=tapply(sust[,2], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s3=tapply(sust[,3], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s4=tapply(sust[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s5=tapply(sust[,5], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s6=tapply(sust[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s7=tapply(sust[,7], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s8=tapply(sust[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s9=tapply(sust[,9], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s10=tapply(sust[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s11=tapply(sust[,11], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s12=tapply(sust[,12], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s13=tapply(sust[,13], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s14=tapply(sust[,14], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s15=tapply(sust[,15], factor(rep(1:numsim, each=numloc)), mean,na.rm=T))

models<-scan("models.txt")
emp<-read.table("Emp.txt")
numsim=1
numloc=15

emp <- data.frame(s1=tapply(emp[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s2=tapply(emp[,2], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s3=tapply(emp[,3], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s4=tapply(emp[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s5=tapply(emp[,5], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s6=tapply(emp[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s7=tapply(emp[,7], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s8=tapply(emp[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s9=tapply(emp[,9], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s10=tapply(emp[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s11=tapply(emp[,11], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s12=tapply(emp[,12], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s13=tapply(emp[,13], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s14=tapply(emp[,14], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                  s15=tapply(emp[,15], factor(rep(1:numsim, each=numloc)), mean,na.rm=T))

pcasust <- prcomp(sust, scale = TRUE)
pcaemp<-predict(pcasust, emp, scale=TRUE)

summary(pcasust)
pcasust<-pcasust$x[,1:4]
pcaemp<-pcaemp[1:4]

cv.modsel1K <- cv4postpr(models[c(1:1000,10001:11000,20001:21000)], sust[c(1:1000,10001:11000,20001:21000),], nval=100, tol=.01, method="neuralnet")
pca.cv.modsel1K <- cv4postpr(models[c(1:1000,10001:11000,20001:21000)], pcasust[c(1:1000,10001:11000,20001:21000),], nval=100, tol=.01, method="neuralnet")
summary(cv.modsel1K)
summary(pca.cv.modsel1K)

cv.modsel2.5K <- cv4postpr(models[c(1:2500,10001:12500,20001:22500)], sust[c(1:2500,10001:12500,20001:22500),], nval=100, tol=.01, method="neuralnet")
pca.cv.modsel2.5K <- cv4postpr(models[c(1:2500,10001:12500,20001:22500)], pcasust[c(1:2500,10001:12500,20001:22500),], nval=100, tol=.01, method="neuralnet")
summary(cv.modsel2.5K)
summary(pca.cv.modsel2.5K)

cv.modsel5K <- cv4postpr(models[c(1:5000,10001:15000,20001:25000)], sust[c(1:5000,10001:15000,20001:25000),], nval=100, tol=.01, method="neuralnet")
pca.cv.modsel5K <- cv4postpr(models[c(1:5000,10001:15000,20001:25000)], pcasust[c(1:5000,10001:15000,20001:25000),], nval=100, tol=.01, method="neuralnet")
summary(cv.modsel5K)
summary(pca.cv.modsel5K)

cv.modsel10K <- cv4postpr(models, sust, nval=100, tol=.01, method="neuralnet")
pca.cv.modsel10K <- cv4postpr(models, pcasust, nval=100, tol=.01, method="neuralnet")
summary(cv.modsel10K)
summary(pca.cv.modsel10K)


PCANN1<-postpr(pcaemp, models, pcasust, tol = 0.01, method = "neuralnet")
NN1<-postpr(emp, models, sust, tol = 0.01, method = "neuralnet")

summary(PCANN1)
summary(NN1)

