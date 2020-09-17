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

## nDNA sample size of Dmelanogaster.
nDNADme = 3
## nDNA sample size of Dsechellia.
nDNADse = 2

## nDNA sample sizes (number of alleles).
nDNANsam = nDNADme + nDNADse

#number of segregating sites for each marker
segsites = [13,91,81]

## create a file to store parameters and one to store the models
parameters = file("parameters.txt","w")
models = file("models.txt","w")

### One Species Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(5,300)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = 0
	
	## ms commands
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 2 %d %d -ej %f 1 2 | sed '/prob/d' | perl msSS.pl >> simModel1.txt" % (nDNANsam, segsites[s],Theta, nDNADme, nDNADse, coalRootDivTime))

	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))
	models.write("1\n")

### Two Species Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(5,300)
	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.01,0.5)

	## ms commands
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 2 %d %d -ej %f 1 2 | sed '/prob/d' | perl msSS.pl >> simModel2.txt" % (nDNANsam, segsites[s],Theta, nDNADme, nDNADse, coalRootDivTime))
	
	## save parameter values and models
	parameters.write("%f\t%f\n" % (Theta, coalRootDivTime))
	models.write("2\n")

os.system("cat simModel* > SuSt_melano_sechellia.txt")