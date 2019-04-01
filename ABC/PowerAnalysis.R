library(abc)
setwd("/Volumes/HD2/manolo/OneDrive/LaGEVol/Artigos_DoutoradoBelManolo/ArtigoDelimitacao/ABC/ABC_CO/IndMarkers/ABC")
sust<-read.table("SuSt.txt")
#numsim=100000
#numloc=26
#sust <- data.frame(tapply(sust[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,2], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,3], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,5], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,7], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,9], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,10], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,11], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,12], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,13], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,14], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
#                   tapply(sust[,15], factor(rep(1:numsim, each=numloc)), mean,na.rm=T))
#write.table(sust, "sust3.txt", sep="\t", row.names=FALSE, col.names = FALSE)

models<-scan("models.txt")

emp<-read.table("Emp.txt")
numsim=1
numloc=20

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

cv.modsel <- cv4postpr(models, sust, nval=100, tol=.005, method="neuralnet")
pca.cv.modsel <- cv4postpr(models, pcasust, nval=100, tol=.005, method="neuralnet")
summary(cv.modsel)
summary(pca.cv.modsel)

PCANN01<-postpr(pcaemp, models, pcasust, tol = 0.001, method = "neuralnet")

PCANN05<-postpr(pcaemp, models, pcasust, tol = 0.005, method = "neuralnet")
NN05<-postpr(emp, models, sust, tol = 0.005, method = "neuralnet")

summary(PCANN05)
summary(NN05)
summary(PCANN01)

