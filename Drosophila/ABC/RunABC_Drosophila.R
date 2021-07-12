##load required libraries
library(abc)

##change working directory if necessary
#setwd("/manolo/GitHub/Drosophila/ABC")

##load models for each simulation
models<-scan("models_melano_simulans.txt")

##load simulated summary statistics
sust<-read.table("SuSt_melano_simulans.txt")

##inform total number of simulations
numsim=200000

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
#                          PC1    PC2    PC3    PC4       PC5       PC6 9.586e-15
#Standard deviation     1.6799 1.3690 1.1168 1.0279 2.621e-07 9.332e-08 0.000e+00
#Proportion of Variance 0.4032 0.2677 0.1782 0.1509 0.000e+00 0.000e+00 1.000e+00
#Cumulative Proportion  0.4032 0.6709 0.8491 1.0000 1.000e+00 1.000e+00 es, as th

##use the first four axes, as they contained at least 99% of the variance
pcasust<-pcasust$x[,1:4]
##use the same axes for the empirical data
pcaemp<-pcaemp[1:4]


##run cross-validation 
start_time <- Sys.time()
pca.cv.modsel10K <- cv4postpr(models, pcasust, nval=100, tol=.2, method="neuralnet")
end_time <- Sys.time()
cvPCANN.2_time = end_time - start_time
cvPCANN.2_time

##visualize the results of the cross validation
summary(pca.cv.modsel10K)
#Output:
#Confusion matrix based on 100 samples for each model.
#
#$tol0.2
#1  2
#1 88 12
#2 25 75
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.2
#1      2
#1 0.7445 0.2555
#2 0.2786 0.7214

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
#  postpr.out$values (40000 posterior samples)
#Models a priori:
#  1, 2
#Models a posteriori:
#  1, 2
#
#Proportion of accepted simulations (rejection):
#  1      2 
#0.2116 0.7884 
#
#Bayes factors:
#  1      2
#1 1.0000 0.2684
#2 3.7253 1.0000
#
#
#Posterior model probabilities (neuralnet):
#  1      2 
#0.0637 0.9363 
#
#Bayes factors:
#  1       2
#1  1.0000  0.0680
#2 14.7032  1.0000

##visualize running time
PCANN.2_time
#Output:
#Time difference of 1.264912 mins

######################################################################
#Now we repeat the same steps for the D.melanogaster-D.sechellia pair#
######################################################################

##load models for each simulation
models<-scan("models_melano_sechellia.txt")

##load simulated summary statistics
sust<-read.table("SuSt_melano_sechellia.txt")

##inform total number of simulations
numsim=200000

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
#                          PC1    PC2    PC3    PC4     PC5      PC6       PC7
#Standard deviation     1.6446 1.3046 1.1898 1.0155 0.38281 1.82e-13 3.769e-15
#Proportion of Variance 0.3864 0.2431 0.2022 0.1473 0.02094 0.00e+00 0.000e+00
#Cumulative Proportion  0.3864 0.6295 0.8317 0.9791 1.00000 1.00e+00 1.000e+00

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
#1 78 22
#2 27 73
#
#
#Mean model posterior probabilities (neuralnet)
#
#$tol0.2
#1      2
#1 0.6882 0.3118
#2 0.3118 0.6882

##run the rejection step of ABC with the empirical data, estimating the running time
start_time <- Sys.time()
PCANN.2<-postpr(pcaemp, models, pcasust, tol = 0.2, method = "neuralnet")
end_time <- Sys.time()
PCANN.2_time = end_time - start_time

##visualize the results of the rejection step
summary(PCANN.2)
#Output:
#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
    method = "neuralnet")
Data:
 postpr.out$values (40000 posterior samples)
Models a priori:
 1, 2
Models a posteriori:
 1, 2

Proportion of accepted simulations (rejection):
     1      2 
0.2439 0.7561 

Bayes factors:
       1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
    1     2 
0.474 0.526 

Bayes factors:
       1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#
#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#
#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#
#
#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#
#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
       method = "neuralnet")
Data:
  postpr.out$values (40000 posterior samples)
Models a priori:
  1, 2
Models a posteriori:
  1, 2

Proportion of accepted simulations (rejection):
  1      2 
0.2439 0.7561 

Bayes factors:
  1      2
1 1.0000 0.3226
2 3.0996 1.0000


Posterior model probabilities (neuralnet):
  1     2 
0.474 0.526 

Bayes factors:
  1      2
1 1.0000 0.9012
2 1.1097 1.0000

#Call: 
#postpr(target = pcaemp, index = models, sumstat = pcasust, tol = 0.2, 
#       method = "neuralnet")
#Data:
#  postpr.out$values (40000 posterior samples)
#Models a priori:
#  1, 2
#Models a posteriori:
#  1, 2
#
#Proportion of accepted simulations (rejection):
#  1      2 
#0.2439 0.7561 
#
#Bayes factors:
#  1      2
#1 1.0000 0.3226
#2 3.0996 1.0000
#
#
#Posterior model probabilities (neuralnet):
#  1     2 
#0.474 0.526 
#
#Bayes factors:
#  1      2
#1 1.0000 0.9012
#2 1.1097 1.0000


##visualize running time
PCANN.2_time
#Output:
#Time difference of 1.981118 mins