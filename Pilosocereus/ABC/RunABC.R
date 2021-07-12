##load required libraries
library(abc)

##change working directory if necessary
#setwd("/manolo/GitHub/Pilosocereus/ABC")

##load models for each simulation
models<-scan("models.txt")

##load simulated summary statistics
sust<-read.table("SuSt.txt")

##inform total number of simulations
numsim=700000

##inform total number of loci
numloc=14

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
#                          PC1    PC2    PC3     PC4     PC5
#Standard deviation     2.9567 1.8505 1.1587 0.38087 0.31719
#Proportion of Variance 0.6244 0.2446 0.0959 0.01036 0.00719
#Cumulative Proportion  0.6244 0.8690 0.9649 0.97529 0.98248
#                           PC6     PC7     PC8     PC9    
#Standard deviation     0.28164 0.23561 0.19058 0.16953
#Proportion of Variance 0.00567 0.00397 0.00259 0.00205
#Cumulative Proportion  0.98814 0.99211 0.99470 0.99675
#                          PC10    PC11    PC12    PC13
#Standard deviation     0.16686 0.11620 0.04617 0.04433
#Proportion of Variance 0.00199 0.00096 0.00015 0.00014
#Cumulative Proportion  0.99874 0.99971 0.99986 1.00000
#                            PC14
#Standard deviation     7.005e-08
#Proportion of Variance 0.000e+00
#Cumulative Proportion  1.000e+00

##use the first seven axes, as they contained at least 99% of the variance
pcasust<-pcasust$x[,1:7]
##use the same axes for the empirical data
pcaemp<-pcaemp[1:7]


##run cross-validation with varying threshold levels, resgistering the running time
pca.cv.modsel100K <- cv4postpr(models, pcasust, nval=100, tols =c(.05,.1,.2), method="neuralnet")

proc.time()
#Output:
#user      system     elapsed 
#1517884.631    1811.109 1519726.666 

#Visualize the results of cross validation with different thresholds
summary(pca.cv.modsel100K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.05
#1  2  3  4  5  6  7
#1 66  0  9  0  0 25  0
#2  0 53  0  4  6  0 37
#3  0  0 98  1  0  0  1
#4  0  1  1 75 19  0  4
#5  0  1  0  2 97  0  0
#6 32  0 18  0  0 50  0
#7  0 39  0  6  5  0 50
#
#$tol0.1
#1  2  3  4  5  6  7
#1 59  0  9  0  0 32  0
#2  0 54  0  3  2  0 41
#3  0  0 98  1  0  0  1
#4  0  1  1 85  8  0  5
#5  0  1  0  2 97  0  0
#6 21  0 19  0  0 60  0
#7  0 38  0  7  0  0 55
#
#$tol0.2
#1  2  3  4  5  6  7
#1 56  0  3  0  0 41  0
#2  0 55  0  4  1  0 40
#3  0  0 97  2  0  0  1
#4  0  1  1 90  3  0  5
#5  0  1  0  2 96  0  1
#6 14  0 10  0  0 76  0
#7  0 39  0  7  0  0 54
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.05
#1      2      3      4      5      6      7
#1 0.5166 0.0017 0.0623 0.0022 0.0000 0.4155 0.0017
#2 0.0001 0.4947 0.0002 0.0524 0.0280 0.0001 0.4246
#3 0.0668 0.0022 0.8129 0.0374 0.0000 0.0759 0.0047
#4 0.0007 0.0745 0.0130 0.7104 0.0934 0.0009 0.1071
#5 0.0000 0.0753 0.0001 0.1686 0.6590 0.0000 0.0970
#6 0.4330 0.0010 0.1137 0.0037 0.0000 0.4470 0.0015
#7 0.0000 0.4269 0.0001 0.0787 0.0258 0.0000 0.4685
#
#$tol0.1
#1      2      3      4      5      6      7
#1 0.5142 0.0010 0.0636 0.0027 0.0000 0.4174 0.0012
#2 0.0003 0.4979 0.0003 0.0451 0.0226 0.0003 0.4335
#3 0.0881 0.0023 0.7661 0.0422 0.0000 0.0968 0.0044
#4 0.0011 0.0634 0.0127 0.7747 0.0502 0.0014 0.0965
#5 0.0000 0.0422 0.0001 0.0947 0.8104 0.0000 0.0526
#6 0.4175 0.0007 0.1016 0.0037 0.0000 0.4754 0.0012
#7 0.0000 0.4292 0.0000 0.0729 0.0096 0.0000 0.4883
#
#$tol0.2
#1      2      3      4      5      6      7
#1 0.6101 0.0001 0.0184 0.0007 0.0000 0.3705 0.0002
#2 0.0004 0.4963 0.0003 0.0505 0.0178 0.0004 0.4344
#3 0.0907 0.0025 0.7538 0.0518 0.0000 0.0963 0.0048
#4 0.0002 0.0571 0.0158 0.8030 0.0323 0.0004 0.0911
#5 0.0000 0.0135 0.0000 0.0417 0.9300 0.0000 0.0148
#6 0.3704 0.0005 0.0526 0.0029 0.0007 0.5719 0.0010
#7 0.0000 0.4297 0.0000 0.0722 0.0096 0.0000 0.4884

##From the above results, we selected a 0.2 threshold

##run cross-validation with varying number of simulations per model
pca.cv.modsel0.5K <- cv4postpr(models[c(1:500,100001:100500,200001:200500,300001:300500,400001:400500,500001:500500,600001:600500)], pcasust[c(1:500,100001:100500,200001:200500,300001:300500,400001:400500,500001:500500,600001:600500),], nval=100, tol=.2, method="neuralnet")
pca.cv.modsel1K <- cv4postpr(models[c(1:1000,100001:101000,200001:201000,300001:301000,400001:401000,500001:501000,600001:601000)], pcasust[c(1:1000,100001:101000,200001:201000,300001:301000,400001:401000,500001:501000,600001:601000),], nval=100, tol=.2, method="neuralnet")
pca.cv.modsel2.5K <- cv4postpr(models[c(1:2500,100001:102500,200001:202500,300001:302500,400001:402500,500001:502500,600001:602500)], pcasust[c(1:2500,100001:102500,200001:202500,300001:302500,400001:402500,500001:502500,600001:602500),], nval=100, tol=.2, method="neuralnet")
pca.cv.modsel10K <- cv4postpr(models[c(1:10000,100001:110000,200001:210000,300001:310000,400001:410000,500001:510000,600001:610000)], pcasust[c(1:10000,100001:110000,200001:210000,300001:310000,400001:410000,500001:510000,600001:610000),], nval=100, tols =.2, method="neuralnet")


##visualize the results of cross validation with  different number of simulations per model used
summary(pca.cv.modsel0.5K)
summary(pca.cv.modsel1K)
summary(pca.cv.modsel2.5K)
summary(pca.cv.modsel10K)

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
#  postpr.out$values (140000 posterior samples)
#Models a priori:
#  1, 2, 3, 4, 5, 6, 7
#Models a posteriori:
#  1, 2, 3, 4, 5, 6, 7
#
#Proportion of accepted simulations (rejection):
#  1      2      3      4      5      6      7 
#0.1728 0.2601 0.0361 0.0714 0.1273 0.1303 0.2019 
#
#Bayes factors:
#  1      2      3      4      5      6      7
#1 1.0000 0.6644 4.7865 2.4198 1.3577 1.3260 0.8561
#2 1.5052 1.0000 7.2045 3.6423 2.0436 1.9958 1.2886
#3 0.2089 0.1388 1.0000 0.5056 0.2837 0.2770 0.1789
#4 0.4133 0.2746 1.9780 1.0000 0.5611 0.5480 0.3538
#5 0.7365 0.4893 3.5254 1.7823 1.0000 0.9766 0.6306
#6 0.7542 0.5011 3.6099 1.8250 1.0240 1.0000 0.6457
#7 1.1680 0.7760 5.5909 2.8265 1.5859 1.5488 1.0000
#
#
#Posterior model probabilities (neuralnet):
#  1      2      3      4      5      6      7 
#0.9626 0.0017 0.0000 0.0000 0.0000 0.0328 0.0029 
#
#Bayes factors:
#  1            2            3            4            5            6
#1 1.000000e+00 5.825617e+02 2.306470e+07 4.176520e+08 2.130273e+19 2.931870e+01
#2 1.700000e-03 1.000000e+00 3.959185e+04 7.169232e+05 3.656734e+16 5.030000e-02
#3 0.000000e+00 0.000000e+00 1.000000e+00 1.810780e+01 9.236077e+11 0.000000e+00
#4 0.000000e+00 0.000000e+00 5.520000e-02 1.000000e+00 5.100594e+10 0.000000e+00
#5 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00
#6 3.410000e-02 1.986990e+01 7.866877e+05 1.424522e+07 7.265908e+17 1.000000e+00
#7 3.100000e-03 1.779600e+00 7.045819e+04 1.275846e+06 6.507573e+16 8.960000e-02
#7
#1 3.273529e+02
#2 5.619000e-01
#3 0.000000e+00
#4 0.000000e+00
#5 0.000000e+00
#6 1.116530e+01
#7 1.000000e+00
##visualize running time
PCANN.2_time
#Output:
#Time difference of 19.83608 mins