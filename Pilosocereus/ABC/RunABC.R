##load required libraries
library(abc)

##change working directory if necessary
#setwd("/manolo/GitHub/Pilosocereus/ABC")

##load models for each simulation
models<-scan("models.txt")

##load simulated summary statistics
sust<-read.table("SuSt.txt")

##inform total number of simulations
numsim=70000

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
#Standard deviation     2.8724 1.8074 1.4145 0.43063 0.31516
#Proportion of Variance 0.5893 0.2333 0.1429 0.01325 0.00709
#Cumulative Proportion  0.5893 0.8226 0.9656 0.97881 0.98590
#                           PC6     PC7     PC8     PC9
#Standard deviation     0.23888 0.20749 0.20650 0.16895
#Proportion of Variance 0.00408 0.00308 0.00305 0.00204
#Cumulative Proportion  0.98998 0.99305 0.99610 0.99814
#                          PC10    PC11    PC12    PC13
#Standard deviation     0.09859 0.09524 0.07373 0.04313
#Proportion of Variance 0.00069 0.00065 0.00039 0.00013
#Cumulative Proportion  0.99883 0.99948 0.99987 1.00000
#                            PC14
#Standard deviation     7.048e-08
#Proportion of Variance 0.000e+00
#Cumulative Proportion  1.000e+00

##use the first seven axes, as they contained at least 99% of the variance
pcasust<-pcasust$x[,1:7]
##use the same axes for the empirical data
pcaemp<-pcaemp[1:7]

##run cross-validation with varying number of simulations per model
pca.cv.modsel10K <- cv4postpr(models, pcasust, nval=100, tols =c(.05,.1,.2), method="neuralnet")
pca.cv.modsel0.5K <- cv4postpr(models[c(1:500,10001:10500,20001:20500,30001:30500,40001:40500,50001:50500,60001:60500)], pcasust[c(1:500,10001:10500,20001:20500,30001:30500,40001:40500,50001:50500,60001:60500),], nval=100, tol=.2, method="neuralnet")
pca.cv.modsel1K <- cv4postpr(models[c(1:1000,10001:11000,20001:21000,30001:31000,40001:41000,50001:51000,60001:61000)], pcasust[c(1:1000,10001:11000,20001:21000,30001:31000,40001:41000,50001:51000,60001:61000),], nval=100, tol=.2, method="neuralnet")
pca.cv.modsel2.5K <- cv4postpr(models[c(1:2500,10001:12500,20001:22500,30001:32500,40001:42500,50001:52500,60001:62500)], pcasust[c(1:2500,10001:12500,20001:22500,30001:32500,40001:42500,50001:52500,60001:62500),], nval=100, tol=.2, method="neuralnet")

##visualize the results of the cross validation with the different number of simulations per model used
summary(pca.cv.modsel10K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.05
#1  2  3  4  5  6  7
#1 48  0 13  0  0 39  0
#2  0 64  0  3 11  0 22
#3  1  0 97  1  0  0  1
#4  0  1  4 52 43  0  0
#5  0  3  0  3 91  0  3
#6 36  0 11  0  0 53  0
#7  1 48  0  5 12  0 34
#
#$tol0.1
#1  2  3  4  5  6  7
#1 42  0 11  0  0 47  0
#2  0 61  0  4  2  0 33
#3  0  0 97  1  0  1  1
#4  0  0  3 61 35  0  1
#5  0  0  0  2 98  0  0
#6 23  0 11  0  0 65  1
#7  0 43  0  4  6  1 46
#
#$tol0.2
#1   2   3   4   5   6   7
#1  48   0   6   0   0  46   0
#2   0  62   0   5   0   0  33
#3   0   0  96   2   0   1   1
#4   0   0   3  89   7   0   1
#5   0   0   0   0 100   0   0
#6  12   0   9   0   0  78   1
#7   0  45   0   4   1   1  49
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.05
#1      2      3      4      5      6      7
#1 0.4581 0.0071 0.0944 0.0041 0.0000 0.4304 0.0059
#2 0.0006 0.5009 0.0002 0.0625 0.0597 0.0004 0.3758
#3 0.0854 0.0032 0.7756 0.0396 0.0000 0.0920 0.0043
#4 0.0018 0.0862 0.0326 0.5405 0.2246 0.0028 0.1114
#5 0.0000 0.1291 0.0003 0.2056 0.5093 0.0001 0.1555
#6 0.4253 0.0087 0.0912 0.0043 0.0006 0.4581 0.0119
#7 0.0052 0.4283 0.0001 0.0714 0.0690 0.0046 0.4214
#
#$tol0.1
#1      2      3      4      5      6      7
#1 0.4755 0.0047 0.0754 0.0040 0.0000 0.4341 0.0063
#2 0.0000 0.5273 0.0000 0.0426 0.0152 0.0001 0.4148
#3 0.0990 0.0031 0.7479 0.0412 0.0000 0.1037 0.0052
#4 0.0012 0.0702 0.0292 0.6587 0.1495 0.0013 0.0900
#5 0.0000 0.0766 0.0004 0.1373 0.6922 0.0000 0.0935
#6 0.4158 0.0064 0.0714 0.0035 0.0000 0.4935 0.0095
#7 0.0050 0.4438 0.0001 0.0559 0.0238 0.0051 0.4663
#
#$tol0.2
#1      2      3      4      5      6      7
#1 0.5428 0.0068 0.0446 0.0073 0.0035 0.3821 0.0129
#2 0.0000 0.5252 0.0000 0.0380 0.0080 0.0000 0.4287
#3 0.0955 0.0031 0.7467 0.0525 0.0000 0.0972 0.0050
#4 0.0001 0.0323 0.0310 0.8377 0.0556 0.0001 0.0432
#5 0.0000 0.0243 0.0001 0.0547 0.8924 0.0000 0.0285
#6 0.3609 0.0062 0.0583 0.0084 0.0042 0.5489 0.0132
#7 0.0059 0.4525 0.0001 0.0433 0.0056 0.0056 0.4869

summary(pca.cv.modsel0.5K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.2
#1  2  3  4  5  6  7
#1 46  0 14  0  0 40  0
#2  0 52  0  3  9  0 36
#3  0  0 98  2  0  0  0
#4  0  4  4 57 35  0  0
#5  0  0  0  1 99  0  0
#6 27  0  7  1  0 64  1
#7  0 41  0  1 13  1 44
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.2
#1      2      3      4      5      6      7
#1 0.4878 0.0147 0.0727 0.0109 0.0023 0.3864 0.0253
#2 0.0005 0.4723 0.0030 0.0484 0.0407 0.0007 0.4344
#3 0.1211 0.0098 0.7042 0.0547 0.0001 0.1028 0.0074
#4 0.0001 0.0921 0.0397 0.6293 0.1456 0.0003 0.0929
#5 0.0000 0.1271 0.0027 0.1569 0.5691 0.0008 0.1435
#6 0.3934 0.0177 0.0472 0.0189 0.0048 0.4853 0.0326
#7 0.0020 0.4156 0.0042 0.0405 0.0558 0.0084 0.4737

summary(pca.cv.modsel1K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.2
#1  2  3  4  5  6  7
#1 44  0 12  0  0 44  0
#2  0 57  0  2  6  1 34
#3  0  0 96  4  0  0  0
#4  0  3  2 55 37  0  3
#5  0  1  0  1 97  0  1
#6 21  0 15  0  0 63  1
#7  0 35  0  3  5  0 57
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.2
#1      2      3      4      5      6      7
#1 0.4951 0.0097 0.0707 0.0100 0.0030 0.3936 0.0178
#2 0.0046 0.5173 0.0018 0.0407 0.0242 0.0054 0.4059
#3 0.1198 0.0017 0.6928 0.0714 0.0000 0.1128 0.0016
#4 0.0001 0.0929 0.0294 0.6171 0.1468 0.0003 0.1136
#5 0.0000 0.0910 0.0011 0.1315 0.6780 0.0002 0.0982
#6 0.3828 0.0117 0.0861 0.0121 0.0037 0.4821 0.0214
#7 0.0001 0.4275 0.0002 0.0398 0.0186 0.0003 0.5135

summary(pca.cv.modsel2.5K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.2
#1  2  3  4  5  6  7
#1 36  0 10  0  0 54  0
#2  0 61  0  5  2  1 31
#3  0  0 97  1  0  2  0
#4  0  1  0 79 19  0  1
#5  0  3  0  1 96  0  0
#6 20  0 14  0  0 66  0
#7  0 33  1  2  2  0 62
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.2
#1      2      3      4      5      6      7
#1 0.4919 0.0079 0.0585 0.0063 0.0018 0.4202 0.0135
#2 0.0022 0.4951 0.0003 0.0542 0.0122 0.0070 0.4290
#3 0.1021 0.0015 0.7415 0.0481 0.0000 0.1054 0.0014
#4 0.0000 0.0495 0.0058 0.7838 0.0973 0.0001 0.0635
#5 0.0000 0.0604 0.0005 0.0774 0.8001 0.0001 0.0614
#6 0.3878 0.0076 0.0817 0.0125 0.0047 0.4908 0.0149
#7 0.0010 0.4169 0.0072 0.0272 0.0105 0.0020 0.5351

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
#  postpr.out$values (14000 posterior samples)
#Models a priori:
#  1, 2, 3, 4, 5, 6, 7
#Models a posteriori:
#  1, 2, 3, 4, 5, 6, 7
#
#Proportion of accepted simulations (rejection):
#  1      2      3      4      5      6      7 
#0.1787 0.2784 0.0311 0.0645 0.1294 0.1247 0.1931 
#
#Bayes factors:
#  1      2      3      4      5      6      7
#1 1.0000 0.6420 5.7385 2.7708 1.3808 1.4330 0.9253
#2 1.5576 1.0000 8.9381 4.3156 2.1507 2.2320 1.4412
#3 0.1743 0.1119 1.0000 0.4828 0.2406 0.2497 0.1612
#4 0.3609 0.2317 2.0711 1.0000 0.4983 0.5172 0.3339
#5 0.7242 0.4650 4.1560 2.0066 1.0000 1.0378 0.6701
#6 0.6978 0.4480 4.0046 1.9336 0.9636 1.0000 0.6457
#7 1.0807 0.6939 6.2018 2.9945 1.4923 1.5487 1.0000
#
#
#Posterior model probabilities (neuralnet):
#  1      2      3      4      5      6      7 
#0.9939 0.0001 0.0000 0.0000 0.0000 0.0060 0.0000 
#
#Bayes factors:
#  1            2            3            4            5
#1 1.000000e+00 1.349494e+04 1.792202e+06 3.009022e+07 1.091754e+08
#2 1.000000e-04 1.000000e+00 1.328055e+02 2.229741e+03 8.090095e+03
#3 0.000000e+00 7.500000e-03 1.000000e+00 1.678950e+01 6.091690e+01
#4 0.000000e+00 4.000000e-04 5.960000e-02 1.000000e+00 3.628300e+00
#5 0.000000e+00 1.000000e-04 1.640000e-02 2.756000e-01 1.000000e+00
#6 6.100000e-03 8.202830e+01 1.089381e+04 1.829018e+05 6.636166e+05
#7 0.000000e+00 3.589000e-01 4.766970e+01 8.003507e+02 2.903886e+03
#6            7
#1 1.645157e+02 3.759629e+04
#2 1.220000e-02 2.786000e+00
#3 1.000000e-04 2.100000e-02
#4 0.000000e+00 1.200000e-03
#5 0.000000e+00 3.000000e-04
#6 1.000000e+00 2.285271e+02
#7 4.400000e-03 1.000000e+00

##visualize running time
PCANN.2_time
#Output:
#Time difference of 1.059996 mins