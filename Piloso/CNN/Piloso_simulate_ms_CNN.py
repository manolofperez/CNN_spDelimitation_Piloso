#!/usr/bin/python

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html.

##import all required modules.
import random
import os
import math
import shlex, subprocess
import numpy as np
from sklearn.neighbors import NearestNeighbors

##define a function to read ms' simulations and transform then into a NumPy array.    
def ms2nparray(xfile):
	g = list(xfile)
	k = [idx for idx,i in enumerate(g) if len(i) > 0 and i.startswith('//')]
	f = []
	for i in k:
	    L = g[i+4:i+nDNANsam+4]
	    q = []
	    for i in L:
	    	i = [int(j) for j in list(i)]
	    	i = np.array(i, dtype=np.int8)
	        q.append(i)
	    q = np.array(q)
	    q = q.astype("int8")
	    f.append(np.array(q))   
	return f

### variable declarations

#define the number of simulations
Priorsize = 10000

###Define sample sizes as number of alleles (2 per specimen).
## nDNA sample size of JFE.
nDNAJFE = 8
## nDNA sample size of ITA.
nDNAITA = 6
## nDNA sample size of PMN.
nDNAPMN = 6
## nDNA sample size of EDB.
nDNAEDB = 6
## nDNA sample size of BOV.
nDNABOV = 8
## nDNA sample size of COC.
nDNACOC = 8
## nDNA sample size of INA.
nDNAINA = 8
## nDNA sample size of ODA.
nDNAODA = 6
## nDNA sample size of MEN.
nDNAMEN = 8
###Combine sample sizes for populations that are merged in one or more models.
## nDNA sample size of ITAPMN.
nDNAITAPMN = nDNAITA + nDNAPMN
## nDNA sample size of BOVEDB.
nDNABOVEDB = nDNAEDB + nDNABOV
## nDNA sample size of Central.
nDNACentral = nDNACOC + nDNAINA + nDNAODA + nDNAMEN
## nDNA sample sizes (number of alleles).
nDNANsam = nDNAITA + nDNAPMN + nDNAEDB + nDNABOV + nDNAJFE + nDNACOC + nDNAINA + nDNAODA + nDNAMEN

## number of years per generation
genlen = 15

#number of segregating sites for each marker
segsites = [26,3,12,9,21,8,13,22,18,7,7,24,14,11]

#initialize a list that will contain the simulations for each model.
simModel1 = []
simModel2 = []
simModel3 = []
simModel4 = []
simModel5 = []

## create a file to store parameters
parameters = file("parameters.txt","w")

### Clade ES Geneland Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior set to 0 in this model.
	coalRootDivTime = 0
	## nDNA ms's command
	com=subprocess.Popen("./ms %d 15 -s 10 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	#simModel1.append(sort_min_diff(np.array(ms2nparray(output)).swapaxes(0,1).reshape(62,-1)).T)
	simModel1.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(72,-1).T)

	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))
	models.write("1\n")

simModel1=np.array(simModel1)
#simModel1=simModel1.swapaxes(1,2)
np.savez_compressed('simModel1.npz', simModel1=simModel1)
del(simModel1)


### Clade ES Splitter Model
for i in range(Priorsize):

	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)
	## nDNA ms's command
	com=subprocess.Popen("./ms %d 15 -s 10 -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalT2, coalT1, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	#simModel2.append(sort_min_diff(np.array(ms2nparray(output)).swapaxes(0,1).reshape(62,-1)))
	simModel2.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(72,-1).T)

	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))
	models.write("2\n")

simModel2=np.array(simModel2)
#simModel2=simModel2.swapaxes(1,2)
np.savez_compressed('simModel2.npz', simModel2=simModel2)

### Clade ES Lumper Model
for i in range(Priorsize):

	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)

	## nDNA ms's command
	com=subprocess.Popen("./ms %d 15 -s 10 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	#simModel2.append(sort_min_diff(np.array(ms2nparray(output)).swapaxes(0,1).reshape(62,-1)))
	simModel3.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(72,-1).T)

	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))
	models.write("3\n")

simModel3=np.array(simModel3)
np.savez_compressed('simModel3.npz', simModel3=simModel3)

