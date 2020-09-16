#!/usr/bin/python

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

import random
import os
import math

### variable declarations

#define the number of simulations
Priorsize = 10000

## nDNA sample size of ITA.
nDNAITA = 8
## nDNA sample size of PMN.
nDNAPMN = 8
## nDNA sample size of EDB.
nDNAEDB = 8
## nDNA sample size of BOV.
nDNABOV = 8
## nDNA sample size of JFE.
nDNAJFE = 8
## nDNA sample size of COC.
nDNACOC = 8
## nDNA sample size of INA.
nDNAINA = 8
## nDNA sample size of ODA.
nDNAODA = 8
## nDNA sample size of MEN.
nDNAMEN = 8

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
## create a file to store parameters and one to store the models
parameters = file("parameters.txt","w")
models = file("models.txt","w")

#Define default values for priors absent in some models.
RootDivTime=0
T1=0
T2=0

### Clade ES Geneland Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior set to 0 in this model.
	coalRootDivTime = 0
	## nDNA ms's command

	## nDNA ms's command
	os.system("./ms %d 15 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4 | perl msSS.pl >> simModel1.txt" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral))

	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, RootDivTime))
	models.write("1\n")

### Clade ES Splitter Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)
	coalT1=random.uniform(0,coalRootDivTime)
	coalT2=random.uniform(0,coalT1)


	## nDNA ms's command
	os.system("./ms %d 15 -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4 | perl msSS.pl >> simModel2.txt" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral, coalT2, coalT1, coalRootDivTime))

	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, RootDivTime))
	models.write("2\n")

### Clade ES Lumper Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)


	## nDNA ms's command
	os.system("./ms %d 15 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4 | perl msSS.pl >> simModel3.txt" % (nDNANsam, Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral, coalRootDivTime))

	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, RootDivTime))
	models.write("3\n")
