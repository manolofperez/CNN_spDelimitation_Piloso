##load required libraries
library(abc)

##change working directory if necessary
#setwd("/Users/manolo/Documents/Novas_spDelPiloso2020/Drosophila/melano_simulans/ABC")

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
#                          PC1    PC2    PC3    PC4     PC5       PC6
#Standard deviation     1.6249 1.2608 1.1164 1.0495 0.64982 1.284e-08
#Proportion of Variance 0.3772 0.2271 0.1780 0.1574 0.06032 0.000e+00
#Cumulative Proportion  0.3772 0.6043 0.7823 0.9397 1.00000 1.000e+00
#                             PC7
#Standard deviation     4.626e-15
#Proportion of Variance 0.000e+00
#Cumulative Proportion  1.000e+00

##use the first five axes, as they contained at least 99% of the variance
pcasust<-pcasust$x[,1:5]
##use the same axes for the empirical data
pcaemp<-pcaemp[1:5]

##run cross-validation 
pca.cv.modsel10K <- cv4postpr(models, pcasust, nval=100, tol=.05, method="neuralnet")

##visualize the results of the cross validation
summary(pca.cv.modsel10K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.05
#   1  2
#1 71 29
#2 27 73
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.05
#       1      2
#1 0.6558 0.3442
#2 0.3323 0.6677

##run the rejection step of ABC with the empirical data, estimating the running time
start_time <- Sys.time()
PCANN.05<-postpr(pcaemp, models, pcasust, tol = 0.05, method = "neuralnet")
end_time <- Sys.time()
PCANN.05_time = end_time - start_time

##visualize the results of the rejection step
summary(PCANN.05)
#Output:
#Call: 
#postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.05, 
#    method = "neuralnet")
#Data:
# postpr.out$values (1000 posterior samples)
#Models a priori:
# 1, 2
#Models a posteriori:
# 1, 2
#
#Proportion of accepted simulations (rejection):
#    1     2 
#0.372 0.628 
#
#Bayes factors:
#       1      2
#1 1.0000 0.5924
#2 1.6882 1.0000
#
#
#Posterior model probabilities (neuralnet):
#     1      2 
#0.3265 0.6735 
#
#Bayes factors:
#       1      2
#1 1.0000 0.4849
#2 2.0624 1.0000

##visualize running time
PCANN.05_time
#Output:
#Time difference of 1.896603 secs

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
#                          PC1    PC2    PC3    PC4     PC5       PC6
#Standard deviation     1.6249 1.2608 1.1164 1.0495 0.64982 1.284e-08
#Proportion of Variance 0.3772 0.2271 0.1780 0.1574 0.06032 0.000e+00
#Cumulative Proportion  0.3772 0.6043 0.7823 0.9397 1.00000 1.000e+00
#                             PC7
#Standard deviation     4.626e-15
#Proportion of Variance 0.000e+00
#Cumulative Proportion  1.000e+00

##use the first five axes, as they contained at least 99% of the variance
pcasust<-pcasust$x[,1:5]
##use the same axes for the empirical data
pcaemp<-pcaemp[1:5]

##run cross-validation 
pca.cv.modsel10K <- cv4postpr(models, pcasust, nval=100, tol=.05, method="neuralnet")

##visualize the results of the cross validation
summary(pca.cv.modsel10K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.05
#   1  2
#1 63 37
#2 32 68
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.05
#       1      2
#1 0.6286 0.3714
#2 0.3335 0.6665

##run the rejection step of ABC with the empirical data, estimating the running time
start_time <- Sys.time()
PCANN.05<-postpr(pcaemp, models, pcasust, tol = 0.05, method = "neuralnet")
end_time <- Sys.time()
PCANN.05_time = end_time - start_time

##visualize the results of the rejection step
summary(PCANN.05)
#Output:
#Call: 
#postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.05, 
#    method = "neuralnet")
#Data:
# postpr.out$values (1000 posterior samples)
#Models a priori:
# 1, 2
#Models a posteriori:
# 1, 2
#
#Proportion of accepted simulations (rejection):
#    1     2 
#0.372 0.628 
#
#Bayes factors:
#       1      2
#1 1.0000 0.5924
#2 1.6882 1.0000
#
#
#Posterior model probabilities (neuralnet):
#     1      2 
#0.0748 0.9252 
#
#Bayes factors:
#        1       2
#1  1.0000  0.0809
#2 12.3684  1.0000

##visualize running time
PCANN.05_time
#Output:
#Time difference of 2.025175 secs