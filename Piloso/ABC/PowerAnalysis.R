##load required libraries
library(abc)

##change working directory if necessary
#setwd("/Users/manolo/Documents/Novas_spDelPiloso2020/ABC")

##load models for each simulation
models<-scan("models.txt")

##load simulated summary statistics
sust<-read.table("SuSt.txt")

##inform total number of simulations
numsim=50000

##inform total number of loci
numloc=16

##average loci summary statistic values
sust <- data.frame(s1=tapply(sust[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
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

##load empirical data
emp<-read.table("Emp.txt")

##set number os datasets to one
numsim=1

##average loci summary statistic values
emp <- data.frame(s1=tapply(emp[,1], factor(rep(1:numsim, each=numloc)), mean,na.rm=T),
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

##Perform a PCA on the simulated datasets
pcasust <- prcomp(sust, scale = TRUE)

##use the calculated eigenvalues to predict PCA values on the empirical dataset
pcaemp<-predict(pcasust, emp, scale=TRUE)

##visualize the summary of the PCA on the simulated dataset
summary(pcasust)
#Output:
#Importance of components:
#                          PC1    PC2    PC3     PC4     PC5    PC6     PC7
#Standard deviation     3.0047 1.6324 1.3263 0.41871 0.36864 0.3107 0.22157
#Proportion of Variance 0.6449 0.1903 0.1256 0.01252 0.00971 0.0069 0.00351
#Cumulative Proportion  0.6449 0.8352 0.9608 0.97336 0.98307 0.9900 0.99347
#                           PC8     PC9    PC10    PC11    PC12      PC13
#Standard deviation     0.20274 0.18439 0.10384 0.06453 0.03661 1.536e-07
#Proportion of Variance 0.00294 0.00243 0.00077 0.00030 0.00010 0.000e+00
#Cumulative Proportion  0.99641 0.99884 0.99961 0.99990 1.00000 1.000e+00
#                           PC14
#Standard deviation     4.99e-08
#Proportion of Variance 0.00e+00
#Cumulative Proportion  1.00e+00

##use the first six axes, as they contained at least 99% of the variance
pcasust<-pcasust$x[,1:6]
##use the same axes for the empirical data
pcaemp<-pcaemp[1:6]

##run cross-validation with varying number of simulations per model
pca.cv.modsel0.5K <- cv4postpr(models[c(1:500,10001:10500,20001:20500,30001:30500,40001:40500)], pcasust[c(1:500,10001:10500,20001:20500,30001:30500,40001:40500),], nval=100, tol=.05, method="neuralnet")
pca.cv.modsel1K <- cv4postpr(models[c(1:1000,10001:11000,20001:21000,30001:31000,40001:41000)], pcasust[c(1:1000,10001:11000,20001:21000,30001:31000,40001:41000),], nval=100, tol=.05, method="neuralnet")
pca.cv.modsel2.5K <- cv4postpr(models[c(1:2500,10001:12500,20001:22500,30001:32500,40001:42500)], pcasust[c(1:2500,10001:12500,20001:22500,30001:32500,40001:42500),], nval=100, tol=.05, method="neuralnet")
pca.cv.modsel10K <- cv4postpr(models, pcasust, nval=100, tol=.05, method="neuralnet")

##visualize the results of the cross validation with the different number of simulations per model used
summary(pca.cv.modsel0.5K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.05
#   1  2  3  4  5
#1 83 11  1  5  0
#2  0 89  0  0 11
#3  0  0 56 44  0
#4  0  0 48 52  0
#5  0  1  0  0 99
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.05
#       1      2      3      4      5
#1 0.8031 0.1119 0.0402 0.0437 0.0010
#2 0.0134 0.8570 0.0064 0.0037 0.1196
#3 0.0291 0.0079 0.4827 0.4782 0.0020
#4 0.0172 0.0077 0.4567 0.5136 0.0048
#5 0.0001 0.1314 0.0111 0.0107 0.8467

summary(pca.cv.modsel1K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.05
#   1  2  3  4  5
#1 88  8  2  2  0
#2  0 90  0  0 10
#3  1  0 45 54  0
#4  0  0 48 52  0
#5  0  1  0  0 99
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.05
#       1      2      3      4      5
#1 0.8497 0.0968 0.0262 0.0273 0.0001
#2 0.0043 0.8826 0.0069 0.0055 0.1008
#3 0.0305 0.0035 0.4620 0.5037 0.0003
#4 0.0144 0.0067 0.4841 0.4901 0.0047
#5 0.0001 0.1095 0.0058 0.0079 0.8766

summary(pca.cv.modsel2.5K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.05
#    1   2   3   4   5
#1  83   7   3   7   0
#2   0  89   0   0  11
#3   0   0  48  52   0
#4   1   0  41  58   0
#5   0   0   0   0 100
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.05
#       1      2      3      4      5
#1 0.8102 0.0944 0.0477 0.0476 0.0000
#2 0.0031 0.8915 0.0008 0.0013 0.1032
#3 0.0177 0.0009 0.4856 0.4944 0.0014
#4 0.0173 0.0005 0.4975 0.4844 0.0004
#5 0.0000 0.1054 0.0061 0.0069 0.8816

summary(pca.cv.modsel10K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.05
#   1  2  3  4  5
#1 92  3  0  5  0
#2  1 92  0  1  6
#3  0  0 55 45  0
#4  0  0 52 48  0
#5  0  1  0  0 99
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.05
#       1      2      3      4      5
#1 0.8615 0.0772 0.0302 0.0312 0.0000
#2 0.0101 0.9288 0.0032 0.0049 0.0531
#3 0.0151 0.0011 0.4902 0.4936 0.0000
#4 0.0136 0.0004 0.4944 0.4915 0.0000
#5 0.0000 0.1134 0.0029 0.0035 0.8803

##run the rejection step of ABC with the empirical data, estimating the running time
start_time <- Sys.time()
PCANN.05<-postpr(pcaemp, models, pcasust, tol = 0.05, method = "neuralnet")
end_time <- Sys.time()
PCANN.05_time = end_time - start_time

##visualize the results of the rejection step
summary(PCANN.05)
#Output:
#postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.05, 
#    method = "neuralnet")
#Data:
# postpr.out$values (2500 posterior samples)
#Models a priori:
# 1, 2, 3, 4, 5
#Models a posteriori:
# 1, 2, 3, 4, 5
#
#Proportion of accepted simulations (rejection):
#     1      2      3      4      5 
#0.1728 0.5476 0.0324 0.0304 0.2168 
#
#Bayes factors:
#        1       2       3       4       5
#1  1.0000  0.3156  5.3333  5.6842  0.7970
#2  3.1690  1.0000 16.9012 18.0132  2.5258
#3  0.1875  0.0592  1.0000  1.0658  0.1494
#4  0.1759  0.0555  0.9383  1.0000  0.1402
#5  1.2546  0.3959  6.6914  7.1316  1.0000
#
#
#Posterior model probabilities (neuralnet):
#     1      2      3      4      5 
#0.9931 0.0058 0.0008 0.0004 0.0000 
#
#Bayes factors:
#            1           2           3           4           5
#1      1.0000    172.2953   1232.0596   2717.1738 589615.8565
#2      0.0058      1.0000      7.1509     15.7704   3422.1241
#3      0.0008      0.1398      1.0000      2.2054    478.5612
#4      0.0004      0.0634      0.4534      1.0000    216.9960
#5      0.0000      0.0003      0.0021      0.0046      1.0000

##visualize running time
PCANN.05_time
#Output:
#Time difference of 7.964716 secs