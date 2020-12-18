##load required libraries
library(abc)

##change working directory if necessary
#setwd("/manolo/GitHub/Drosophila/ABC")

##load models for each simulation
models<-scan("models_melano_simulans.txt")

##load simulated summary statistics
sust<-read.table("SuSt_melano_simulans.txt")

##inform total number of simulations
numsim=20000

##inform total number of loci
numloc=3

##average loci summary statistic values
sust <- data.frame(s1=tapply(sust[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s3=tapply(sust[,3], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s4=tapply(sust[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s5=tapply(sust[,5], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s6=tapply(sust[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s7=tapply(sust[,7], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s8=tapply(sust[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T))

##load empirical data
emp<-read.table("Emp_melano_simulans.txt")

##set number os datasets to one
numsim=1

##average loci summary statistic values
emp <- data.frame(s1=tapply(emp[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s3=tapply(emp[,3], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s4=tapply(emp[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s5=tapply(emp[,5], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s6=tapply(emp[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s7=tapply(emp[,7], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s8=tapply(emp[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T))

##Perform a PCA on the simulated datasets
pcasust <- prcomp(sust, scale = TRUE)

##use the calculated eigenvalues to predict PCA values on the empirical dataset
pcaemp<-predict(pcasust, emp, scale=TRUE)

##visualize the summary of the PCA on the simulated dataset
summary(pcasust)
#Output:
#Importance of components:
#                          PC1    PC2    PC3    PC4       PC5       PC6       PC7
#Standard deviation     1.6827 1.3705 1.1114 1.0272 2.606e-07 9.219e-08 1.264e-14
#Proportion of Variance 0.4045 0.2683 0.1764 0.1507 0.000e+00 0.000e+00 0.000e+00
#Cumulative Proportion  0.4045 0.6728 0.8493 1.0000 1.000e+00 1.000e+00 1.000e+00

##use the first four axes, as they contained at least 99% of the variance
pcasust<-pcasust$x[,1:4]
##use the same axes for the empirical data
pcaemp<-pcaemp[1:4]

##run cross-validation 
pca.cv.modsel10K <- cv4postpr(models, pcasust, nval=100, tol=.2, method="neuralnet")

##visualize the results of the cross validation
summary(pca.cv.modsel10K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.2
#1  2
#1 82 18
#2 16 84
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.2
#1      2
#1 0.7072 0.2928
#2 0.2069 0.7931

##run the rejection step of ABC with the empirical data, estimating the running time
start_time <- Sys.time()
PCANN.2<-postpr(pcaemp, models, pcasust, tol = 0.2, method = "neuralnet")
end_time <- Sys.time()
PCANN.2_time = end_time - start_time

##visualize the results of the rejection step
summary(PCANN.2)
#Output:
#Call: 
#  postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
#         method = "neuralnet")
#Data:
#  postpr.out$values (4000 posterior samples)
#Models a priori:
#  1, 2
#Models a posteriori:
#  1, 2
#
#Proportion of accepted simulations (rejection):
#  1      2 
#0.2145 0.7855 
#
#Bayes factors:
#  1      2
#1 1.0000 0.2731
#2 3.6620 1.0000
#
#
#Posterior model probabilities (neuralnet):
#  1      2 
#0.0175 0.9825 
#
#Bayes factors:
#  1       2
#1  1.0000  0.0178
#2 56.2607  1.0000

##visualize running time
PCANN.2_time
#Output:
#Time difference of 5.037071 secs

######################################################################
#Now we repeat the same steps for the D.melanogaster-D.sechellia pair#
######################################################################

##load models for each simulation
models<-scan("models_melano_sechellia.txt")

##load simulated summary statistics
sust<-read.table("SuSt_melano_sechellia.txt")

##inform total number of simulations
numsim=20000

##inform total number of loci
numloc=3

##average loci summary statistic values
sust <- data.frame(s1=tapply(sust[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s3=tapply(sust[,3], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s4=tapply(sust[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s5=tapply(sust[,5], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s6=tapply(sust[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s7=tapply(sust[,7], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s8=tapply(sust[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T))

##load empirical data
emp<-read.table("Emp_melano_sechellia.txt")

##set number os datasets to one
numsim=1

##average loci summary statistic values
emp <- data.frame(s1=tapply(emp[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s3=tapply(emp[,3], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s4=tapply(emp[,4], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s5=tapply(emp[,5], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s6=tapply(emp[,6], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s7=tapply(emp[,7], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
                   s8=tapply(emp[,8], factor(rep(1:numsim, each=numloc)), mean,na.rm=T))

##Perform a PCA on the simulated datasets
pcasust <- prcomp(sust, scale = TRUE)

##use the calculated eigenvalues to predict PCA values on the empirical dataset
pcaemp<-predict(pcasust, emp, scale=TRUE)

##visualize the summary of the PCA on the simulated dataset
summary(pcasust)
#Output:
#Importance of components:
#                          PC1    PC2    PC3    PC4     PC5     PC6       PC7
#Standard deviation     1.5883 1.3163 1.1698 1.0219 0.57613 1.7e-14 2.343e-15
#Proportion of Variance 0.3604 0.2475 0.1955 0.1492 0.04742 0.0e+00 0.000e+00
#Cumulative Proportion  0.3604 0.6079 0.8034 0.9526 1.00000 1.0e+00 1.000e+00

##use the first five axes, as they contained at least 99% of the variance
pcasust<-pcasust$x[,1:5]
##use the same axes for the empirical data
pcaemp<-pcaemp[1:5]

##run cross-validation 
pca.cv.modsel10K <- cv4postpr(models, pcasust, nval=100, tol=.2, method="neuralnet")

##visualize the results of the cross validation
summary(pca.cv.modsel10K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.2
#1  2
#1 81 19
#2 27 73
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.2
#1      2
#1 0.7221 0.2779
#2 0.3488 0.6512

##run the rejection step of ABC with the empirical data, estimating the running time
start_time <- Sys.time()
PCANN.2<-postpr(pcaemp, models, pcasust, tol = 0.2, method = "neuralnet")
end_time <- Sys.time()
PCANN.2_time = end_time - start_time

##visualize the results of the rejection step
summary(PCANN.2)
#Output:
#Call: 
#  postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
#         method = "neuralnet")
#Data:
#  postpr.out$values (4000 posterior samples)
#Models a priori:
#  1, 2
#Models a posteriori:
#  1, 2
#
#Proportion of accepted simulations (rejection):
#  1      2 
#0.5312 0.4688 
#
#Bayes factors:
#  1      2
#1 1.0000 1.1333
#2 0.8824 1.0000
#
#
#Posterior model probabilities (neuralnet):
#  1      2 
#0.8939 0.1061 
#
#Bayes factors:
#  1      2
#1 1.0000 8.4261
#2 0.1187 1.0000

##visualize running time
PCANN.2_time
#Output:
#Time difference of 5.045474 secs