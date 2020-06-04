#!/usr/bin/python

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

import random
import os
import math
import shlex, subprocess
import numpy as np
from sklearn.neighbors import NearestNeighbors

def sort_min_diff(amat):
    '''this function takes in a SNP matrix with indv on rows and returns the same matrix with indvs sorted by genetic similarity.
    this problem is NP-hard, so here we use a nearest neighbors approx.  it's not perfect, but it's fast and generally performs ok.
    assumes your input matrix is a numpy array'''
    mb = NearestNeighbors(len(amat), metric='manhattan').fit(amat)
    v = mb.kneighbors(amat)
    smallest = np.argmin(v[0].sum(axis=1))
    return amat[v[1][smallest]]
    
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
Priorsize = 100000

## nDNA sample size of ART2.
nDNAART2 = 8
## nDNA sample size of ART1.
nDNAART1 = 8
## nDNA sample size of CRI.
nDNACRI = 8
## nDNA sample size of GIL.
nDNAGIL = 8
## nDNA sample size of MOC.
nDNAMOC = 8
## nDNA sample size of POS2.
nDNAPOS2 = 8
## nDNA sample size of UNA.
nDNAUNA = 8
## nDNA sample size of BUR.
nDNABUR = 8
## nDNA sample size of DBO.
nDNADBO = 8
## nDNA sample size of FMS.
nDNAFMS = 8
## nDNA sample size of APA1-4.
nDNAAPAs = 32
## nDNA sample size of MAG and TEG.
nDNAMAG_TEG = 16
## nDNA sample size of URU.
nDNAURU = 8
## nDNA sample size of PIR.
nDNAPIR = 8
## nDNA sample size of GOV.
nDNAGOV = 8
## nDNA sample size of FOR and DEL.
nDNAFOR_DEL = 16
## nDNA sample size of PET and ALC.
nDNAPET_ALC = 16
## nDNA sample size of MIN.
nDNAMIN = 8
## nDNA sample size of PGO.
nDNAPGO = 8
## nDNA sample size of AQU and RVE.
nDNAAQU_RVE = 16


## nDNA sample sizes (number of alleles).
nDNANsam = nDNAART2 + nDNAART1 + nDNACRI + nDNAGIL + nDNAMOC + nDNAPOS2 + nDNAUNA + nDNABUR + nDNADBO + nDNAFMS + nDNAAPAs + nDNAMAG_TEG + nDNAURU + nDNAPIR + nDNAGOV + nDNAFOR_DEL + nDNAPET_ALC + nDNAMIN + nDNAPGO + nDNAAQU_RVE
## number of years per generation
genlen = 15

simModel1 = []
simModel2 = []
simModel3 = []
## create a file to store parameters and one to store the models
parameters = file("parameters.txt","w")
models = file("models.txt","w")

#Define default values for priors absent in some models.
T5=0
T6=0

### Clade CO Lumper Model
for i in range(Priorsize):
	### Define parameters
	## Ne prior following a uniform distribution from 50 to 10000
	Ne = random.uniform(50000, 500000)
	## mutation rate according to Ossowski et al. (2010)
	mutrate =(7.0E-9)
	## use Ne and mutrate values to obtain the value of theta (required by ms)
	Theta = 4*Ne*mutrate*450
	## divergence time prior following an uniform distribution from 650k to 350k years ago.
	RootDivTime = random.uniform(350000, 650000)
	## subsequent divergence events with prior following an uniform distribution from the time of the previous event until the present.
	T4=random.uniform(0,RootDivTime)
	T3=random.uniform(0,T4)
	T2=random.uniform(0,T3)
	T1=random.uniform(0,T2)
	## use the DivTime in years to calculte divergence time in coalescent units (required by ms)
	coalRootDivTime = RootDivTime/(genlen*4*Ne)
	coalT4 = T4/(genlen*4*Ne)
	coalT3 = T3/(genlen*4*Ne)
	coalT2 = T2/(genlen*4*Ne)
	coalT1 = T1/(genlen*4*Ne)
	## nDNA ms's command
	com=subprocess.Popen("./ms %d 26 -s 5 -t %f -I 20 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d -ej 0 19 5 -ej 0 13 5 -ej 0 9 5 -ej 0 15 4 -ej 0 14 4 -ej 0 20 7 -ej 0 16 7 -ej 0 8 7 -ej 0 10 6 -ej 0 12 2 -ej 0 11 2 -ej %f 5 4 -ej %f 18 1 -ej %f 3 2 -ej %f 17 1 -ej %f 6 2 -ej %f 2 1 -ej %f 4 1 -ej %f 7 1" % (nDNANsam, Theta,  nDNAGIL, nDNAART1, nDNAART2, nDNAMAG_TEG, nDNAAPAs, nDNADBO, nDNAPET_ALC, nDNAMIN, nDNAFOR_DEL, nDNABUR, nDNAUNA, nDNAFMS, nDNAURU, nDNAGOV, nDNAPIR, nDNAAQU_RVE, nDNAMOC, nDNAPOS2, nDNACRI, nDNAPGO, coalT1, coalT1, coalT1, coalT2, coalT2, coalT3, coalT4, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel1.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Ne, RootDivTime, T1, T2, T3, T4, T5, T6))
	models.write("1\n")

simModel1=np.array(simModel1)
np.savez_compressed('simModel1.npz', simModel1=simModel1)
del(simModel1)

### Clade CO Splitter Model
for i in range(Priorsize):

	### Define parameters
	## Ne prior following a uniform distribution from 50 to 10000
	Ne = random.uniform(50000, 500000)
	## mutation rate according to Ossowski et al. (2010)
	mutrate =(7.0E-9)
	## use Ne and mutrate values to obtain the value of theta (required by ms)
	Theta = 4*Ne*mutrate*450
	## divergence time prior following an uniform distribution from 650k to 350k years ago.
	RootDivTime = random.uniform(350000, 650000)
	## subsequent divergence events with prior following an uniform distribution from the time of the previous event until the present.
	T5=random.uniform(0,RootDivTime)
	T4=random.uniform(0,T5)
	T3=random.uniform(0,T4)
	T2=random.uniform(0,T3)
	T1=random.uniform(0,T2)

	## use the DivTime in years to calculte divergence time in coalescent units (required by ms)
	coalRootDivTime = RootDivTime/(genlen*4*Ne)
	coalT5 = T5/(genlen*4*Ne)
	coalT4 = T4/(genlen*4*Ne)
	coalT3 = T3/(genlen*4*Ne)
	coalT2 = T2/(genlen*4*Ne)
	coalT1 = T1/(genlen*4*Ne)

	## nDNA ms's command
	com=subprocess.Popen("./ms %d 26 -s 5 -t %f -I 20 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d -ej 0 15 14 -ej 0 13 9 -ej 0 20 8 -ej 0 10 6 -ej 0 12 2 -ej %f 5 4 -ej %f 8 7 -ej %f 18 1 -ej %f 3 2 -ej %f 14 7 -ej %f 9 4 -ej %f 17 1 -ej %f 16 7 -ej %f 6 2 -ej %f 11 1 -ej %f 19 2 -ej %f 7 4 -ej %f 4 2 -ej %f 2 1" % (nDNANsam, Theta,  nDNAGIL, nDNAART1, nDNAART2, nDNAMAG_TEG, nDNAAPAs, nDNADBO, nDNAPET_ALC, nDNAMIN, nDNAFOR_DEL, nDNABUR, nDNAUNA, nDNAFMS, nDNAURU, nDNAGOV, nDNAPIR, nDNAAQU_RVE, nDNAMOC, nDNAPOS2, nDNACRI, nDNAPGO, coalT1, coalT1, coalT1, coalT1, coalT1, coalT2, coalT2, coalT2, coalT2, coalT3, coalT3, coalT4, coalT5, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel2.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Ne, RootDivTime, T1, T2, T3, T4, T5, T6))
	models.write("2\n")

simModel2=np.array(simModel2)
np.savez_compressed('simModel2.npz', simModel2=simModel2)
del(simModel2)


### Clade CO Geneland Model
for i in range(Priorsize):

	### Define parameters

	### Define parameters
	## Ne prior following a uniform distribution from 50 to 10000
	Ne = random.uniform(50000, 500000)
	## mutation rate according to Ossowski et al. (2010)
	mutrate =(7.0E-9)
	## use Ne and mutrate values to obtain the value of theta (required by ms)
	Theta = 4*Ne*mutrate*450
	## divergence time prior following an uniform distribution from 650k to 350k years ago.
	RootDivTime = random.uniform(350000, 650000)
	## subsequent divergence events with prior following an uniform distribution from the time of the previous event until the present.
	T6=random.uniform(0,RootDivTime)
	T5=random.uniform(0,T6)
	T4=random.uniform(0,T5)
	T3=random.uniform(0,T4)
	T2=random.uniform(0,T3)
	T1=random.uniform(0,T2)

	## use the DivTime in years to calculte divergence time in coalescent units (required by ms)
	coalRootDivTime = RootDivTime/(genlen*4*Ne)
	coalT6 = T6/(genlen*4*Ne)
	coalT5 = T5/(genlen*4*Ne)
	coalT4 = T4/(genlen*4*Ne)
	coalT3 = T3/(genlen*4*Ne)
	coalT2 = T2/(genlen*4*Ne)
	coalT1 = T1/(genlen*4*Ne)

	## nDNA ms's command
	com=subprocess.Popen("./ms %d 26 -s 5 -t %f -I 20 %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d -ej 0 3 2 -ej 0 12 6 -ej %f 19 1 -ej %f 18 1 -ej %f 10 6 -ej %f 5 4 -ej %f 15 13 -ej %f 20 8 -ej %f 17 11 -ej %f 13 4 -ej %f 16 8 -ej %f 11 1 -ej %f 14 4 -ej %f 8 7 -ej %f 6 1 -ej %f 9 7 -ej %f 4 1 -ej %f 7 1 -ej %f 2 1" % (nDNANsam, Theta,  nDNAGIL, nDNAART1, nDNAART2, nDNAMAG_TEG, nDNAAPAs, nDNADBO, nDNAPET_ALC, nDNAMIN, nDNAFOR_DEL, nDNABUR, nDNAUNA, nDNAFMS, nDNAURU, nDNAGOV, nDNAPIR, nDNAAQU_RVE, nDNAMOC, nDNAPOS2, nDNACRI, nDNAPGO, coalT1, coalT1, coalT1, coalT1, coalT1, coalT1, coalT2, coalT2, coalT2, coalT3, coalT3, coalT3, coalT4, coalT4, coalT5, coalT6, coalRootDivTime), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	simModel3.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	## save parameter values
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Ne, RootDivTime, T1, T2, T3, T4, T5, T6))
	models.write("3\n")

simModel3=np.array(simModel3)
np.savez_compressed('simModel3.npz', simModel3=simModel3)