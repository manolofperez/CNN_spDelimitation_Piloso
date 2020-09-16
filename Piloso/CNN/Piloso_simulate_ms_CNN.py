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

### Splitter Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	#subsequent divergence times, from until the previous divergence times, making sure they are younger than the root.
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)

	## migration priors set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDB=0
	mJFE_Central=0
	mITAPMN_BOVEDB=0
	mITAPMN_Central=0
	mBOVEDB_Central=0
	
	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalT2, coalT1, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)
	
	## cpDNA markers
	com=subprocess.Popen("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalT2/2, coalT1/2, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	cp_output = np.array(ms2nparray(com.read().splitlines()))
	## mtDNA markers
	com=subprocess.Popen("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalT2/2, coalT1/2, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	mt_output = np.array(ms2nparray(com.read().splitlines()))

	#concatenate segregating sites from all markers
	simModel1.append(np.concatenate((nc_output, cp_output, mt_output),axis=2).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))

#save NumPy arrays
simModel1=np.array(simModel1)
np.save('simModel1.npy', simModel1)
del(simModel1)


### Splitter Model + Migration
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	#subsequent divergence times, from until the previous divergence times, making sure they are younger than the root.
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)

	## migration priors following an uniform distribution from 0 to 5.
	mJFE_ITAPMN=random.uniform(0,5)
	mJFE_BOVEDB=random.uniform(0,5)
	mJFE_Central=random.uniform(0,5)
	mITAPMN_BOVEDB=random.uniform(0,5)
	mITAPMN_Central=random.uniform(0,5)
	mBOVEDB_Central=random.uniform(0,5)

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -m 1 2 %f -m 2 1 %f -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 2 3 %f -m 3 2 %f -m 2 4 %f -m 4 2 %f -m 3 4 %f -m 4 3 %f -ej %f 2 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_BOVEDB,  mJFE_Central, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_BOVEDB, mITAPMN_Central, mITAPMN_Central, mBOVEDB_Central, mBOVEDB_Central, coalT2, coalT1, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	## cpDNA markers
	com=subprocess.Popen("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -m 1 2 %f -m 2 1 %f -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 2 3 %f -m 3 2 %f -m 2 4 %f -m 4 2 %f -m 3 4 %f -m 4 3 %f -ej %f 2 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_BOVEDB,  mJFE_Central, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_BOVEDB, mITAPMN_Central, mITAPMN_Central, mBOVEDB_Central, mBOVEDB_Central, coalT2/2, coalT1/2, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	cp_output = np.array(ms2nparray(com.read().splitlines()))
	## mtDNA markers
	com=subprocess.Popen("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -m 1 2 %f -m 2 1 %f -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 2 3 %f -m 3 2 %f -m 2 4 %f -m 4 2 %f -m 3 4 %f -m 4 3 %f -ej %f 2 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_BOVEDB,  mJFE_Central, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_BOVEDB, mITAPMN_Central, mITAPMN_Central, mBOVEDB_Central, mBOVEDB_Central, coalT2/2, coalT1/2, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	mt_output = np.array(ms2nparray(com.read().splitlines()))

	#concatenate segregating sites from all markers
	simModel2.append(np.concatenate((nc_output, cp_output, mt_output),axis=2).swapaxes(0,1).reshape(nDNANsam,-1).T)

	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))

#save NumPy arrays
simModel2=np.array(simModel2)
np.save('simModel2.npy', simModel2)
del(simModel2)

### Lumper Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)

	## migration priors set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDB=0
	mJFE_Central=0
	mITAPMN_BOVEDB=0
	mITAPMN_Central=0
	mBOVEDB_Central=0

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	## cpDNA markers
	com=subprocess.Popen("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	cp_output = np.array(ms2nparray(com.read().splitlines()))
	## mtDNA markers
	com=subprocess.Popen("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	mt_output = np.array(ms2nparray(com.read().splitlines()))

	#concatenate segregating sites from all markers
	simModel3.append(np.concatenate((nc_output, cp_output, mt_output),axis=2).swapaxes(0,1).reshape(nDNANsam,-1).T)

	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))

#save NumPy arrays
simModel3=np.array(simModel3)
np.save('simModel3.npy', simModel3)
del(simModel3)

### Lumper Model + Migration
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)

	## migration prior between population JFE and the remaining following an uniform distribution from 0 to 5.
	mJFE_Central=random.uniform(0,5)

	## Other migration priors set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDB=0
	mITAPMN_BOVEDB=0
	mITAPMN_Central=0
	mBOVEDB_Central=0

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -m 1 4 %f -m 4 1 %f -ej %f 1 4" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, mJFE_Central, mJFE_Central, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	## cpDNA markers
	com=subprocess.Popen("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -m 1 4 %f -m 4 1 %f -ej %f 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, mJFE_Central, mJFE_Central, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	cp_output = np.array(ms2nparray(com.read().splitlines()))
	## mtDNA markers
	com=subprocess.Popen("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -m 1 4 %f -m 4 1 %f -ej %f 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, mJFE_Central, mJFE_Central, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	mt_output = np.array(ms2nparray(com.read().splitlines()))

	#concatenate segregating sites from all markers
	simModel4.append(np.concatenate((nc_output, cp_output, mt_output),axis=2).swapaxes(0,1).reshape(nDNANsam,-1).T)

	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))

#save NumPy arrays
simModel4=np.array(simModel4)
np.save('simModel4.npy', simModel4)
del(simModel4)

### Geneland Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior set to 0 in this model.
	coalRootDivTime = 0

	## migration prior set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDB=0
	mJFE_Central=0
	mITAPMN_BOVEDB=0
	mITAPMN_Central=0
	mBOVEDB_Central=0

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	## cpDNA markers
	com=subprocess.Popen("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral), shell=True, stdout=subprocess.PIPE).stdout
	cp_output = np.array(ms2nparray(com.read().splitlines()))
	## mtDNA markers
	com=subprocess.Popen("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral), shell=True, stdout=subprocess.PIPE).stdout
	mt_output = np.array(ms2nparray(com.read().splitlines()))

	#concatenate segregating sites from all markers
	simModel5.append(np.concatenate((nc_output, cp_output, mt_output),axis=2).swapaxes(0,1).reshape(nDNANsam,-1).T)

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))
	models.write("5\n")

#save NumPy arrays
simModel5=np.array(simModel5)
np.save('simModel5.npy', simModel5)