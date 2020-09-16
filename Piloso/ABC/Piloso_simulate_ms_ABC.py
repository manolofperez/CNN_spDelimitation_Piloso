#!/usr/bin/python

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

##import all required modules.
import random
import os
import math

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

#number of segregating sites for each nuclear marker
segsites = [26,3,12,9,21,8,13,22,18,7,7,24,14,11]

## create a file to store parameters and one to store the models
parameters = file("parameters.txt","w")
models = file("models.txt","w")

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

	## migration prior set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDB=0
	mJFE_Central=0
	mITAPMN_BOVEDB=0
	mITAPMN_Central=0
	mBOVEDB_Central=0
	
	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel1.txt" % (nDNANsam, segsites[s],Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalT2, coalT1, coalRootDivTime))

	## cpDNA markers
	com=os.system("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel1.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalT2/2, coalT1/2, coalRootDivTime/2))
	## mtDNA markers
	com=os.system("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -ej %f 2 4 -ej %f 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel1.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalT2/2, coalT1/2, coalRootDivTime/2))

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))
	models.write("1\n")

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

	## divergence time prior following an uniform distribution from 0 to 5.
	mJFE_ITAPMN=random.uniform(0,5)
	mJFE_BOVEDB=random.uniform(0,5)
	mJFE_Central=random.uniform(0,5)
	mITAPMN_BOVEDB=random.uniform(0,5)
	mITAPMN_Central=random.uniform(0,5)
	mBOVEDB_Central=random.uniform(0,5)

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -m 1 2 %f -m 2 1 %f -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 2 3 %f -m 3 2 %f -m 2 4 %f -m 4 2 %f -m 3 4 %f -m 4 3 %f -ej %f 2 4 -ej %f 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel2.txt" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_BOVEDB,  mJFE_Central, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_BOVEDB, mITAPMN_Central, mITAPMN_Central, mBOVEDB_Central, mBOVEDB_Central, coalT2, coalT1, coalRootDivTime))
	
	## cpDNA markers
	com=os.system("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -m 1 2 %f -m 2 1 %f -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 2 3 %f -m 3 2 %f -m 2 4 %f -m 4 2 %f -m 3 4 %f -m 4 3 %f -ej %f 2 4 -ej %f 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel2.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_BOVEDB,  mJFE_Central, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_BOVEDB, mITAPMN_Central, mITAPMN_Central, mBOVEDB_Central, mBOVEDB_Central, coalT2/2, coalT1/2, coalRootDivTime/2))
	## mtDNA markers
	com=os.system("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -m 1 2 %f -m 2 1 %f -m 1 3 %f -m 3 1 %f -m 1 4 %f -m 4 1 %f -m 2 3 %f -m 3 2 %f -m 2 4 %f -m 4 2 %f -m 3 4 %f -m 4 3 %f -ej %f 2 4 -ej %f 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel2.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB, nDNACentral, mJFE_ITAPMN, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_BOVEDB,  mJFE_Central, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_BOVEDB, mITAPMN_Central, mITAPMN_Central, mBOVEDB_Central, mBOVEDB_Central, coalT2/2, coalT1/2, coalRootDivTime/2))

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))
	models.write("2\n")

### Lumper Model
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)

	## migration prior set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDB=0
	mJFE_Central=0
	mITAPMN_BOVEDB=0
	mITAPMN_Central=0
	mBOVEDB_Central=0

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel3.txt" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalRootDivTime))
	
	## cpDNA markers
	com=os.system("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel3.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalRootDivTime/2))
	## mtDNA markers
	com=os.system("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel3.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, coalRootDivTime/2))

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))
	models.write("3\n")
	
### Lumper Model + Migration
for i in range(Priorsize):

	### Define parameters
	## Theta values from 1 to 15
	Theta = random.uniform(1,15)

	## divergence time prior following an uniform distribution from 0.5 to 5.
	coalRootDivTime = random.uniform(0.5,5)

	## migration prior set to 0 in this model.
	mJFE_ITAPMN=0
	mJFE_BOVEDB=0
	mJFE_Central=0
	mITAPMN_BOVEDB=0
	mITAPMN_Central=0
	mBOVEDB_Central=0

	### ms commands
	##nuclear markers, using the number of segregating sites for each loci
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -m 1 4 %f -m 4 1 %f -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel4.txt" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, mJFE_Central, mJFE_Central, coalRootDivTime))
	
	## cpDNA markers
	com=os.system("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -m 1 4 %f -m 4 1 %f -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel4.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, mJFE_Central, mJFE_Central, coalRootDivTime/2))
	## mtDNA markers
	com=os.system("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -m 1 4 %f -m 4 1 %f -ej %f 1 4 | sed '/prob/d' | perl msSS.pl >> simModel4.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral, mJFE_Central, mJFE_Central, coalRootDivTime/2))

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))
	models.write("4\n")
	
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

	## ms commands
	for s in range(len(segsites)):
		## nDNA markers
		com=os.system("./ms %d 1 -s %d -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4 | sed '/prob/d' | perl msSS.pl >> simModel5.txt" % (nDNANsam, segsites[s], Theta, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral))
	
	## cpDNA markers
	com=os.system("./ms %d 1 -s 15 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4 | sed '/prob/d' | perl msSS.pl >> simModel5.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral))
	## mtDNA markers
	com=os.system("./ms %d 1 -s 4 -t %f -I 4 %d %d %d %d -ej 0 2 4 -ej 0 3 4 -ej 0 1 4 | sed '/prob/d' | perl msSS.pl >> simModel5.txt" % (nDNANsam, Theta/4, nDNAJFE, nDNAITAPMN ,nDNABOVEDB , nDNACentral))

	## save parameter values and models
	parameters.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, coalRootDivTime, coalT1, coalT2, mJFE_ITAPMN, mJFE_BOVEDB, mJFE_Central, mITAPMN_BOVEDB, mITAPMN_Central, mBOVEDB_Central))
	models.write("5\n")
	
com=os.system("cat simModel* > SuSt.txt")