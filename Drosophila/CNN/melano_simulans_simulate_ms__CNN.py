#!/usr/bin/python

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

##import all required modules.
import random
import os
import math
import shlex, subprocess
import numpy as np

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

## nDNA sample size of Dmelanogaster.
nDNADme = 3
## nDNA sample size of Dsimulans.
nDNADsi = 3

## combined sample sizes (number of alleles).
nDNANsam = nDNADme + nDNADsi

#number of segregating sites for each nuclear marker
segsites = [36,36]

#initialize a list that will contain the simulations for each model.
simModel1 = []
simModel2 = []

## create a file to store parameters.
parameters = file("parameters.txt","w")

### One Species Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(5,300)
	## divergence time prior set to 0.
	coalRootDivTime = 0
	
	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 2 %d %d -ej %f 1 2" % (nDNANsam, segsites[s],Theta, nDNADme, nDNADsi, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)
	
	## mtDNA marker
	com=subprocess.Popen("./ms %d 1 -s 553 -t %f -I 2 %d %d -ej %f 1 2" % (nDNANsam, Theta/4, nDNADme, nDNADsi, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	mt_output = np.array(ms2nparray(com.read().splitlines()))

	#concatenate segregating sites from all markers
	simModel1.append(np.concatenate((nc_output, mt_output),axis=2).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))

#save NumPy arrays
simModel1=np.array(simModel1)
np.savez_compressed('simModel1.npz', simModel1=simModel1)
del(simModel1)

### Two Species Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 5 to 300
	Theta = random.uniform(5,300)

	## divergence time prior following an uniform distribution from 0.01 to 0.5.
	coalRootDivTime = random.uniform(0.01,0.5)

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	nc_output = np.empty([1,nDNANsam,0])
	for s in range(len(segsites)):
		## nDNA markers
		com=subprocess.Popen("./ms %d 1 -s %d -t %f -I 2 %d %d -ej %f 1 2" % (nDNANsam, segsites[s],Theta, nDNADme, nDNADsi, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
		nc_output = np.append(nc_output, np.array(ms2nparray(com.read().splitlines())),axis=2)

	## mtDNA marker
	com=subprocess.Popen("./ms %d 1 -s 553 -t %f -I 2 %d %d -ej %f 1 2" % (nDNANsam, Theta/4, nDNADme, nDNADsi, coalRootDivTime/2), shell=True, stdout=subprocess.PIPE).stdout
	mt_output = np.array(ms2nparray(com.read().splitlines()))

	#concatenate segregating sites from all markers
	simModel2.append(np.concatenate((nc_output, mt_output),axis=2).swapaxes(0,1).reshape(nDNANsam,-1).T)

	## save parameter values
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))

#save NumPy arrays
simModel2=np.array(simModel2)
np.savez_compressed('simModel2.npz', simModel2=simModel2)