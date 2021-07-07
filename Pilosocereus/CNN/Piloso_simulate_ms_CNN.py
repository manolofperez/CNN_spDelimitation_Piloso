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
import time


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
## nDNA sample size of BOVEDBODA.
nDNABOVEDB = nDNAEDB + nDNABOV
## nDNA sample size of BOVEDBODA.
nDNABOVEDBODA = nDNAEDB + nDNABOV + nDNAODA
## nDNA sample size of INAMEN.
nDNAINAMEN = nDNAINA + nDNAMEN
## nDNA sample sizes (number of alleles).
nDNANsam = nDNAITA + nDNAPMN + nDNAEDB + nDNABOV + nDNAJFE + nDNACOC + nDNAINA + nDNAODA + nDNAMEN

## number of years per generation
genlen = 15

#number of segregating sites for each nuclear marker
segsites = [26,3,12,9,21,8,13,22,18,7,7,24,14,11]

#initialize a list that will contain the simulations for each model.
simModel1 = []
simModel2 = []
simModel3 = []
simModel4 = []
simModel5 = []
simModel6 = []
simModel7 = []

## create a file to store parameters
parameters = file("parameters.txt","w")

start = time.time()
## Splitter Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	#subsequent divergence times, from 0 until the previous divergence times, making sure they are younger than the root.
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)
	coalT3=random.uniform(0,coalT2)


	## migration prior set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDBODA=0
	mJFE_INAMEN=0
	mJFE_COC=0
	mITAPMN_BOVEDBODA=0
	mITAPMN_INAMEN=0
	mITAPMN_COC=0
	mBOVEDBODA_INAMEN=0
	mBOVEDBODA_COC=0
	mINAMEN_COC=0
	
	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 5 %d %d %d %d %d -ej %f 2 4 -ej %f 5 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDBODA, nDNAINAMEN, nDNACOC, coalT3, coalT2, coalT1, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)
	
	#Transpose the matrix
	simModel1.append(nc_output.swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))

#save NumPy arrays
simModel1=np.array(simModel1)
np.save('simModel1.npy', simModel1)
del(simModel1)
print ('Time: ')
print (time.time() - start)

start = time.time()
### Splitter Model + Migration
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	#subsequent divergence times, from 0 until the previous divergence times, making sure they are younger than the root.
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)
	coalT3=random.uniform(0,coalT2)

	## divergence time prior following an uniform distribution from 0 to 5.
	mJFE_ITAPMN=random.uniform(0,5)
	mJFE_BOVEDBODA=random.uniform(0,5)
	mJFE_INAMEN=random.uniform(0,5)
	mJFE_COC=random.uniform(0,5)
	mITAPMN_BOVEDBODA=random.uniform(0,5)
	mITAPMN_INAMEN=random.uniform(0,5)
	mITAPMN_COC=random.uniform(0,5)
	mBOVEDBODA_INAMEN=random.uniform(0,5)
	mBOVEDBODA_COC=random.uniform(0,5)
	mINAMEN_COC=random.uniform(0,5)

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 5 %d %d %d %d %d -m 1 2 %f -m 2 1 %f -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 1 5 %f -m 5 1 %f -m 2 3 %f -m 3 2 %f -m 2 4 %f -m 4 2 %f -m 2 5 %f -m 5 2 %f -m 3 4 %f -m 4 3 %f -m 3 5 %f -m 5 3 %f -m 4 5 %f -m 5 4 %f -ej %f 2 4 -ej %f 5 4 -ej %f 3 4 -ej %f 1 4" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDBODA, nDNAINAMEN, nDNACOC, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_BOVEDBODA,  mJFE_INAMEN, mJFE_INAMEN, mJFE_COC, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_INAMEN, mITAPMN_COC, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mBOVEDBODA_COC, mINAMEN_COC, mINAMEN_COC, coalT3, coalT2, coalT1, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)
	
	#Transpose the matrix
	simModel2.append(nc_output.swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))

#save NumPy arrays
simModel2=np.array(simModel2)
np.save('simModel2.npy', simModel2)
del(simModel2)
print ('Time: ')
print (time.time() - start)

start = time.time()
### Lumper Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	#subsequent divergence times, set to 0 in this model.
	coalT1=0
	coalT2=0
	coalT3=0

	## migration prior set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDBODA=0
	mJFE_INAMEN=0
	mJFE_COC=0
	mITAPMN_BOVEDBODA=0
	mITAPMN_INAMEN=0
	mITAPMN_COC=0
	mBOVEDBODA_INAMEN=0
	mBOVEDBODA_COC=0
	mINAMEN_COC=0

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 5 %d %d %d %d %d -ej 0 2 4 -ej 0 5 4 -ej 0 3 4 -ej %f 1 4" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDBODA, nDNAINAMEN, nDNACOC, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	#Transpose the matrix
	simModel3.append(nc_output.swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))

#save NumPy arrays
simModel3=np.array(simModel3)
np.save('simModel3.npy', simModel3)
del(simModel3)
print ('Time: ')
print (time.time() - start)

start = time.time()
### Lumper Model + Migration
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	#subsequent divergence times, set to 0 in this model.
	coalT1=0
	coalT2=0
	coalT3=0

	## migration prior, set to 0 in this model, except for the migration between JFE and the remaining populations.
	mJFE_ITAPMN=0
	mJFE_BOVEDBODA=random.uniform(0,5)
	mJFE_INAMEN=0
	mJFE_COC=0
	mITAPMN_BOVEDBODA=0
	mITAPMN_INAMEN=0
	mITAPMN_COC=0
	mBOVEDBODA_INAMEN=0
	mBOVEDBODA_COC=0
	mINAMEN_COC=0

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 5 %d %d %d %d %d -ej 0 2 4 -ej 0 5 4 -ej 0 3 4 -m 1 4 %f -m 4 1 %f -ej %f 1 4" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDBODA, nDNAINAMEN, nDNACOC, mJFE_BOVEDBODA, mJFE_BOVEDBODA, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	#Transpose the matrix
	simModel4.append(nc_output.swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))

#save NumPy arrays
simModel4=np.array(simModel4)
np.save('simModel4.npy', simModel4)
del(simModel4)
print ('Time: ')
print (time.time() - start)

start = time.time()
### Geneland Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior set to 0 in this model.
	coalRootDivTime = 0
	#subsequent divergence times, set to 0 in this model.
	coalT1=0
	coalT2=0
	coalT3=0

	## migration prior set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDBODA=0
	mJFE_INAMEN=0
	mJFE_COC=0
	mITAPMN_BOVEDBODA=0
	mITAPMN_INAMEN=0
	mITAPMN_COC=0
	mBOVEDBODA_INAMEN=0
	mBOVEDBODA_COC=0
	mINAMEN_COC=0

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 5 %d %d %d %d %d -ej 0 2 4 -ej 0 5 4 -ej 0 3 4 -ej 0 1 4" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDBODA, nDNAINAMEN, nDNACOC), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	#Transpose the matrix
	simModel5.append(nc_output.swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))

#save NumPy arrays
simModel5=np.array(simModel5)
np.save('simModel5.npy', simModel5)
del(simModel5)
print ('Time: ')
print (time.time() - start)

start = time.time()
### BPP_noGDI
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	#subsequent divergence times, from 0 until the previous divergence times, making sure they are younger than the root.
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)
	coalT3=0


	## migration prior set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDBODA=0
	mJFE_INAMEN=0
	mJFE_COC=0
	mITAPMN_BOVEDBODA=0
	mITAPMN_INAMEN=0
	mITAPMN_COC=0
	mBOVEDBODA_INAMEN=0
	mBOVEDBODA_COC=0
	mINAMEN_COC=0
	
	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 2 5 -ej %f 3 4 -ej %f 6 5 -ej %f 4 5 -ej %f 1 5" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC, coalT3, coalT2, coalT1, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	#Transpose the matrix
	simModel6.append(nc_output.swapaxes(0,1).reshape(nDNANsam,-1).T)

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))

#save NumPy arrays
simModel6=np.array(simModel6)
np.save('simModel6.npy', simModel6)
del(simModel6)
print ('Time: ')
print (time.time() - start)

start = time.time()
### BPP_noGDI + Migration
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	#subsequent divergence times, from 0 until the previous divergence times, making sure they are younger than the root.
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)
	coalT3=0

	## divergence time prior following an uniform distribution from 0 to 5.
	mJFE_ITAPMN=0
	mJFE_BOVEDBODA=random.uniform(0,5)
	mJFE_INAMEN=random.uniform(0,5)
	mJFE_COC=random.uniform(0,5)
	mITAPMN_BOVEDBODA=0
	mITAPMN_INAMEN=0
	mITAPMN_COC=0
	mBOVEDBODA_INAMEN=random.uniform(0,5)
	mBOVEDBODA_COC=random.uniform(0,5)
	mINAMEN_COC=random.uniform(0,5)

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 2 5 -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 1 5 %f -m 5 1 %f -m 1 6 %f -m 6 1 %f -m 3 4 %f -m 4 3 %f -m 3 5 %f -m 5 3 %f -m 3 6 %f -m 6 3 %f -m 4 5 %f -m 5 4 %f -m 4 6 %f -m 6 4 %f -m 5 6 %f -m 6 5 %f -ej %f 3 4 -ej %f 6 5 -ej %f 4 5 -ej %f 1 5" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_BOVEDBODA,  mJFE_INAMEN, mJFE_INAMEN, mJFE_COC, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_INAMEN, mITAPMN_COC, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mBOVEDBODA_COC, mINAMEN_COC, mINAMEN_COC, coalT3, coalT2, coalT1, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	#Transpose the matrix
	simModel7.append(nc_output.swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))

#save NumPy arrays
simModel7=np.array(simModel7)
np.save('simModel7.npy', simModel7)
print ('Time: ')
print (time.time() - start)
