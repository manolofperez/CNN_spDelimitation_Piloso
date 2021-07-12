#!/usr/bin/python

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

##import all required modules.
import random
import os
import math
import time

### variable declarations

#define the number of simulations
Priorsize = 100000

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
#nDNABOVEDBODA = nDNAEDB + nDNABOV + nDNAODA
## nDNA sample size of BOVEDBODA.
nDNABOVEDB = nDNAEDB + nDNABOV
## nDNA sample size of INAMEN.
nDNAINAMEN = nDNAINA + nDNAMEN
## nDNA sample sizes (number of alleles).
nDNANsam = nDNAITA + nDNAPMN + nDNAEDB + nDNABOV + nDNAJFE + nDNACOC + nDNAINA + nDNAODA + nDNAMEN

## number of years per generation
genlen = 15

#number of segregating sites for each nuclear marker
segsites = [26,3,12,9,21,8,13,22,18,7,7,24,14,11]

## create a file to store parameters and one to store the models
parameters = file("parameters.txt","w")
models = file("models.txt","w")

start = time.time()
### Splitter Model
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
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 3 4 -ej %f 2 5 -ej %f 6 5 -ej %f 4 5 -ej %f 1 5 | sed '/prob/d' | perl msSS.pl >> simModel1.txt" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC, coalT3, coalT2, coalT1, coalRootDivTime))

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))
	models.write("1\n")
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
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 3 4 -m 1 2 %f -m 2 1 %f -m 1 4 %f -m 4 1 %f -m 1 5 %f -m 5 1 %f -m 1 6 %f -m 6 1 %f -m 2 4 %f -m 4 2 %f -m 2 5 %f -m 5 2 %f -m 2 6 %f -m 6 2 %f -m 4 5 %f -m 5 4 %f -m 4 6 %f -m 6 4 %f -m 5 6 %f -m 6 5 %f -ej %f 2 5 -ej %f 6 5 -ej %f 4 5 -ej %f 1 5 | sed '/prob/d' | perl msSS.pl >> simModel2.txt" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_BOVEDBODA,  mJFE_INAMEN, mJFE_INAMEN, mJFE_COC, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_INAMEN, mITAPMN_COC, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mBOVEDBODA_COC, mINAMEN_COC, mINAMEN_COC, coalT3, coalT2, coalT1, coalRootDivTime))
	
	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))
	models.write("2\n")

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
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 3 4 -ej 0 2 5 -ej 0 6 5 -ej 0 4 5 -ej %f 1 5 | sed '/prob/d' | perl msSS.pl >> simModel3.txt" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC, coalRootDivTime))
	
	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))
	models.write("3\n")

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
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 3 4 -ej 0 2 5 -ej 0 6 5 -ej 0 4 5 -m 1 5 %f -m 5 1 %f -ej %f 1 5 | sed '/prob/d' | perl msSS.pl >> simModel4.txt" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC, mJFE_BOVEDBODA, mJFE_BOVEDBODA, coalRootDivTime))
	
	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))
	models.write("4\n")

print ('Time: ')
print (time.time() - start)

start = time.time()	
### BPP+GDI Model
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

	## ms commands
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 3 4 -ej 0 2 5 -ej 0 6 5 -ej 0 4 5 -ej 0 1 5 | sed '/prob/d' | perl msSS.pl >> simModel5.txt" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC))
	
	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))
	models.write("5\n")

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
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 2 5 -ej %f 3 4 -ej %f 6 5 -ej %f 4 5 -ej %f 1 5 | sed '/prob/d' | perl msSS.pl >> simModel6n.txt" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC, coalT3, coalT2, coalT1, coalRootDivTime))

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))
	models.write("6\n")

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
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 6 %d %d %d %d %d %d -ej 0 2 5 -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 1 5 %f -m 5 1 %f -m 1 6 %f -m 6 1 %f -m 3 4 %f -m 4 3 %f -m 3 5 %f -m 5 3 %f -m 3 6 %f -m 6 3 %f -m 4 5 %f -m 5 4 %f -m 4 6 %f -m 6 4 %f -m 5 6 %f -m 6 5 %f -ej %f 3 4 -ej %f 6 5 -ej %f 4 5 -ej %f 1 5 | sed '/prob/d' | perl msSS.pl >> simModel7n.txt" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNAODA, nDNAINAMEN, nDNACOC, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_BOVEDBODA,  mJFE_INAMEN, mJFE_INAMEN, mJFE_COC, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_INAMEN, mITAPMN_COC, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mBOVEDBODA_COC, mINAMEN_COC, mINAMEN_COC, coalT3, coalT2, coalT1, coalRootDivTime))
	
	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, coalT3, mJFE_ITAPMN, mJFE_BOVEDBODA, mJFE_INAMEN, mJFE_COC, mITAPMN_BOVEDBODA, mITAPMN_INAMEN, mITAPMN_COC, mBOVEDBODA_INAMEN, mBOVEDBODA_COC, mINAMEN_COC))
	models.write("7\n")

print ('Time: ')
print (time.time() - start)
	
os.system("cat simModel* > SuSt.txt")