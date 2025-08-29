#Date: 27th of Nov, 2024
#Location: Blacksburg, VA (USA)
#Purpose: Genotype strains from Metastrainer CHain 8
#Input: Metastrainer output, source reference, linkage group file
#output: N strain fasta files

#Location: Blacksburg, VA and Raleigh, NC (USA)
#Date: Over several weeks
#Only print strains if pairwise identity is more than threshold, in descending order
#Input: Read MetaStrainer key file and default name Non-Reference 100% alleles
#Also produce a new genotype composition file (with new order of strain name regardless)

#Date: 6th of May, 2025
#Location: Raleigh, NC (USA)
#Purpose: Debug and Track Non-Reference allele genotyping
#Input: Nothing New
#output: Log identical reference matches to Error File for debugging



#Changes
#Date: Jan13

import os
import argparse
import sys
import copy
import traceback
import time

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

timestr = time.strftime("%m%d%Y_%H%M")

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", required=True, help="Louis-Marie's MCMC output")
parser.add_argument("-r","--reference", required=True,help="Reference fasta file")
parser.add_argument("-g","--genemap", required=True, help="map matrix family name to core gene id in cluster file")
parser.add_argument("-l","--linkgage", required=True, default="linkage_groups.txt",help="gene linkage file")
parser.add_argument("-o","--output",default="SampleX_Statistics_",help="Name of output statistics")
#parser.add_argument("-N","--Number", default=3, type=int,help="Number of strain peaks in MetaStrainer composition file")
parser.add_argument("-S","--SampleName",default="",help="SampleName to name output strains as per Louis-Marie request")
parser.add_argument("--singleton",default="../Preprocess/singles_singleton_noref.txt",help="Non-reference alleles file")
parser.add_argument("-s","--StrainThreshold",default=99.5,type=float,help="Threshold for defining a strain for merging")
parser.add_argument("-f","--FrequencyFile",required=True,help="File with initial sample composition")
parser.add_argument('--logfile',help="Logging Errors",default="Genotyping_errorLog.txt")
args = parser.parse_args()


#python /nas/longleaf/home/hazem/ReconAssScripts/GenotypeMetaStrainEight.py -i ../genotypes_$1.txt  -o Genotype -S $1 -f ../key_genotypes_$1.txt -g /nas/longleaf/home/hazem/core_genomes/Bifidobacterium_asteroides_core2/Basteroides_JACFPA01_family2gene_mapping_150.txt -r /nas/longleaf/home/hazem/core_genomes/Bifidobacterium_asteroides_core2/core_Bifidobacterium_asteroides_JACFPA01_0.fasta -l ../linkage_groups.txt
#	

if not os.path.exists(args.singleton):
	sys.exit("Error reading Singleton Non-Reference allele file %s"%(args.singleton))

geneMap = {}
with open(args.genemap) as genemap_file:
	lines = genemap_file.readline()
	while lines:
		lineContents = lines.strip("\n").split("\t")
		#Map to original core genes (no flank)
		geneMap[lineContents[1]] = lineContents[0]
		lines = genemap_file.readline()


sampleFrequency  = {}
#Not a smart way but I restricted MetaStrainer to 3
strainNum = 3 
with open(args.FrequencyFile) as input_FrequencyFile:
	lines = input_FrequencyFile.readline()
	while lines:
		if strainNum == 0:
			system.exit("Fix this script. MetaStrainer has been upgraded to produce a variable number of strains")
		lineContents = lines.strip("\n").split("\t")
		strainName = "strain" + str(strainNum)
		sampleFrequency[strainName] = float(lineContents[1])
		strainNum = strainNum - 1
		lines = input_FrequencyFile.readline()



#Get reference sequence (0-based)
referenceSeq = dict()
seqID = ""
sequence = ""
totalRefLength = 0#for updated stats genomeLength
for record in SeqIO.parse(args.reference, "fasta"):
	seqID = record.id
	sequence =record.seq
	totalRefLength = totalRefLength + len(sequence)
	referenceSeq[seqID] = list(sequence)


#Deconstruct positions back into original based on LM script
contigGroup={}
f=open(args.linkgage,"r")
for l in f:
	a=l.strip("\n").split("\t")
	contigGroup[a[0]]={}
	position=0
	for name in a[1:]:
		contigGroup[a[0]][name]=position
		sub=name.split("_")
		deb,fin=int(sub[-2]),int(sub[-1])
		position= position + (fin - deb)
		if (a[0]=="group651"):
			print("%s\t\tName %s Position %s"%(a[0],name,position))
f.close()
counter = 0
trackID = {}
for gene in contigGroup:
	counter = counter + len(contigGroup[gene].keys())
print("We have %s Genes in %s contigs "%(counter,len(contigGroup.keys())))
# for name in contigGroup:
# 	for gene in contigGroup[name]:
# 		print(name,"\t",gene,"\t",contigGroup[name][gene])

def convertGroupToPairGenes(contigGroupName,pos):
	genePointer = ""
	genePos = 0
	for gene in contigGroup[contigGroupName]:
		if (contigGroup[contigGroupName][gene] == 0):
			genePointer = gene
			genePos = contigGroup[contigGroupName][gene]
		#if (contigGroupName == "group2" or contigGroupName == "group703"):
			#print("Inside Converter")
			#print(genePointer,"\t",genePos)
		#should be <= to account for last position
		if (pos <= contigGroup[contigGroupName][gene]):
			#print("returning early")
			return(genePointer, pos - genePos)
		else:
			genePointer = gene
			genePos = contigGroup[contigGroupName][gene]
	
	#if loops ends and no return then position is in last gene/contig sequentially.
	genePointer = gene
	genePos = contigGroup[contigGroupName][gene]
	if (contigGroupName == "group6"  or contigGroupName == "group703"):
		pass
		#Debugging Conversion from Group/Contig back to Genes
		#print("Inside Converter End of Loop")
		#print(genePointer,"\t",genePos)

	return(genePointer, pos - genePos)

strain1 = copy.deepcopy(referenceSeq)
strain2 = copy.deepcopy(referenceSeq)
strain3 = copy.deepcopy(referenceSeq)


comitted_changes_count = N_comitted_changes_count = 0
strain1Changes = 0
strain2Changes = 0
strain3Changes = 0
ReferenceChange = {}
line_counter = 0
debugN_writer = open("Debug_N_genotype.txt","w")
stayN = 0
counterEncounter = 0
UnKnownErrors = 0
ErrorLoggerName = args.logfile + "_ErrorLog_" + timestr + ".txt"
ErrorLogger = open(ErrorLoggerName,"w",buffering = 20)
debugExceptionCount = 0
with open(args.input) as input_file:
	lines = input_file.readline()
	while lines:
		line_counter = line_counter + 1
		#group1	group1&83&G	0.5628140703517588	3	333333
		#group1	group1&126&G	0.9504310344827587	6	6666666666
		#group1	group1&83&T	0.4371859296482412	4	444444


		lineContents = lines.strip("\n").split("\t")
		
		#geneID = lineContents[0]
		alleleInfo = lineContents[1].split("&")

		#quick and dirty to accomodate single alleles (nopair)
		if (len(alleleInfo) == 3):
			groupID = alleleInfo[0]
			groupPosition = int(alleleInfo[1])
			allele = alleleInfo[2]
		elif(len(alleleInfo) == 4):
			groupID = alleleInfo[1]
			groupPosition = int(alleleInfo[2])
			allele = alleleInfo[3]

		#bug in output of MetaStrainer: temporary handling
		try:
			peak = int(lineContents[3])
		except:
			
			ErrorLogger.write("Peak Exception Error")
			ErrorLogger.write("MetaStrainer output: " + lines)
			ErrorLogger.write("lineContents[2]\t\t"+ lineContents[2] + "\n")
			ErrorLogger.write("lineContents[3](peak)\t"+ lineContents[2] + "\n")
			ErrorLogger.write("Probability (4)\t\t"+ lineContents[3] + "\n")
			
			debugExceptionCount = debugExceptionCount + 1
			lines = input_file.readline()
			continue
			

		probability = float(lineContents[4])
		if (len(alleleInfo) == 3):
			geneID,position = convertGroupToPairGenes(groupID,groupPosition)
		elif (len(alleleInfo) == 4):
			geneID = groupID
			position = groupPosition
		if geneID not in trackID:
			trackID[geneID] = "seen"

		#Debugging Conversion from Group/Contig back to Genes
		#if (groupID == "group6" or groupID =="group703"):
		#	print("We just printed Coverter with %s and got %s"%(groupPosition_1,position_1))
		#	print(lines)
		#	print("geneID1: ",geneID1,"\tposition1: ",position_1,"\tgeneID2: ",geneID2,"\tposition_2: ",position_2)

		##########################
		##########################
		###Please check MetaStrainer to above


		if geneID not in ReferenceChange:
			ReferenceChange[geneID] = {}
		if position not in ReferenceChange[geneID]:
			ReferenceChange[geneID][position] = {}
		
		if peak not in ReferenceChange[geneID][position]:
			if (probability >= 50):
				ReferenceChange[geneID][position][peak] = allele
			else:
				ReferenceChange[geneID][position][peak] = "N"
		else:
			ReferenceChange[geneID][position][peak] = "N"
		lines = input_file.readline()

comitted_changes_count = N_comitted_changes_count = commited_allstrains_count = comitted_allN_count = 0
strain1Changes =  N_strain1Changes = 0
strain2Changes =  N_strain2Changes = 0
strain3Changes =  N_strain3Changes = 0

#Apply non-reference alleles to all strains before writing other changes
singletonNonRefAlleleCounter = 0
#Count false Non-reference alleles if they are already identical to Reference!
falseNonRefAllele = 0
with open(args.singleton,"r") as single_noRef:
	lines = single_noRef.readline()
	while lines:
		lineContents = lines.strip("\n").split("\t")
		geneID = lineContents[0]
		position = int(lineContents[1])#This file is 0-based so do not adjust
		allele = lineContents[2]

		strain1[geneID][position] = allele
		strain2[geneID][position] = allele
		strain3[geneID][position] = allele

		#If Allele is identical to Reference then count and check alignment/BAM
		if (allele == referenceSeq[geneID][position]):
			falseNonRefllele = falseNonRefllele + 1
			ErrorLogger.write("NonRef Allele Singleton Error")
			ErrorLogger.log(lines)



		singletonNonRefAlleleCounter = singletonNonRefAlleleCounter + 1
		lines = single_noRef.readline()		



print("Loaded Singleton non-reference variant: %s" %(singletonNonRefAlleleCounter))
print("Wrong  Singleton non-reference variants: %s" %(falseNonRefAllele))
variantCounter = 0

geneCount = 0
for gene in strain2:
	geneCount+=1
	#print(gene,"\t",geneCount)
print("Keys strain1\t",len(strain1))
print("Keys strain2\t",len(strain2))
print("Keys strain3\t",len(strain3))
print("Keys Reference\t",len(ReferenceChange))
print("Debug count: Skipped lines with MS output error %s" %(debugExceptionCount))
print("Debug count: out of read %s" %(line_counter))

#print("Eh da Ref\t",ReferenceChange["MZNG01000001.1_40622_41242"])
strain12_sim = 0
strain13_sim = 0
strain23_sim = 0
peakCount = 0
#try:
for gene in ReferenceChange:
	for position in ReferenceChange[gene]:
		variantCounter = variantCounter + 1
		for peak in ReferenceChange[gene][position]:
			variantBase = ReferenceChange[gene][position][peak]
			peakCount = peakCount + 1
			try:
				if peak == 1:
					strain3[gene][position] = variantBase
				elif peak == 2:
					strain2[gene][position] = variantBase
				elif peak == 3:
					strain1[gene][position] = variantBase
				elif peak == 4:
					strain3[gene][position] = variantBase
					strain2[gene][position] = variantBase
				elif peak == 5:
					strain3[gene][position] = variantBase
					strain1[gene][position] = variantBase
				elif peak == 6:
					strain2[gene][position] = variantBase
					strain1[gene][position] = variantBase
				else:
					print("Peak = %s"%(peak))
					sys.exit("Input dictionary\\process is corrupt")
			except:
					print(gene," ",(position + 1)," length ",len(strain1[gene]))
					sys.exit("Input dictionary\\process is corrupt")


		if strain1[gene][position] != referenceSeq[gene][position]:
			if strain1[gene][position] == "N":
				N_strain1Changes = N_strain1Changes + 1
			else:
				strain1Changes = strain1Changes + 1

		if strain2[gene][position] != referenceSeq[gene][position]:
			if strain2[gene][position] == "N":
				N_strain2Changes = N_strain2Changes + 1
			else:
				strain2Changes = strain2Changes + 1

		if strain3[gene][position] != referenceSeq[gene][position]:
			if strain3[gene][position] == "N":
				N_strain3Changes = N_strain3Changes + 1
			else:
				strain3Changes = strain3Changes + 1

		if (strain1[gene][position]  == "N" or strain2[gene][position]  == "N" or strain3[gene][position]  == "N"):
			N_comitted_changes_count = N_comitted_changes_count + 1
			if (strain1[gene][position] == strain2[gene][position] == strain1[gene][position]) and strain1[gene][position] == "N":
				comitted_allN_count = comitted_allN_count + 1


		#if ((referenceSeq[gene][position - 1] != strain1[gene][position - 1] ) or (referenceSeq[gene][position - 1] != strain2[gene][position - 1]) or (referenceSeq[gene][position - 1] != strain3[gene][position - 1]) ) and (referenceSeq[gene][position - 1]  != "N"):
		if (strain1[gene][position]  != "N" and strain2[gene][position]  != "N" and strain3[gene][position]  != "N"):
			if ((referenceSeq[gene][position] != strain1[gene][position] ) or (referenceSeq[gene][position] != strain2[gene][position]) or (referenceSeq[gene][position] != strain3[gene][position]) ):
				comitted_changes_count = comitted_changes_count + 1

		if (strain1[gene][position] == strain2[gene][position]):
			strain12_sim = strain12_sim + 1

		if (strain1[gene][position] == strain3[gene][position]):
			strain13_sim = strain13_sim + 1

		if (strain2[gene][position] == strain3[gene][position]):
			strain23_sim = strain23_sim + 1

print("Strain1 Vs. Strain2 (by Variants) %s \t(by total ref) %s"%( (strain12_sim/variantCounter) * 100 ,(strain12_sim/totalRefLength) * 100 ) )
print("Strain1 Vs. Strain3 (by Variants) %s \t(by total ref) %s"%( (strain13_sim/variantCounter) * 100 ,(strain13_sim/totalRefLength) * 100 ) )
print("Strain2 Vs. Strain3 (by Variants) %s \t(by total ref) %s"%( (strain23_sim/variantCounter) * 100 ,(strain23_sim/totalRefLength) * 100 ) )
#sys.exit("Done Early to know strain similarlity")

#genomw-wide identical variants 
#strain_1v2 = (strain12_sim/totalRefLength) * 100
#strain_1v3 = (strain13_sim/totalRefLength) * 100
#strain_2v3 = (strain23_sim/totalRefLength) * 100

wholeRefID = totalRefLength - variantCounter
strain_1v2 = ((strain12_sim + wholeRefID)/totalRefLength) * 100
strain_1v3 = ((strain13_sim + wholeRefID)/totalRefLength) * 100
strain_2v3 = ((strain23_sim + wholeRefID)/totalRefLength) * 100

print("Debug totalRef %s" %(totalRefLength))
print("Debug whole Ref %s" %(wholeRefID))
print("Debug variantCounter %s" %(variantCounter))


print("Debug Pairwise final strain 1v2  (%s)" %(strain_1v2))
print("Debug Pairwise final strain 1v3  (%s)" %(strain_1v3))
print("Debug Pairwise final strain 2v3  (%s)" %(strain_2v3))


#identical variants 
strain_1v2_identicalVariants = (strain12_sim/variantCounter) * 100
strain_1v3_identicalVariants = (strain13_sim/variantCounter) * 100
strain_2v3_identicalVariants = (strain23_sim/variantCounter) * 100


print("Debug identical Strain1  (%s)" %((strain_1v2_identicalVariants + strain_1v3_identicalVariants)))
print("Debug identical Strain2  (%s)" %((strain_1v2_identicalVariants + strain_2v3_identicalVariants)))
print("Debug identical Strain3  (%s)" %((strain_1v3_identicalVariants + strain_2v3_identicalVariants)))

if (strain_1v2 >= args.StrainThreshold and strain_1v3 >= args.StrainThreshold and strain_2v3 >= args.StrainThreshold):
	print("Printing Only One Strain")
	sampleFrequency['strain1'] = sampleFrequency['strain1'] + sampleFrequency['strain2'] + sampleFrequency['strain3']
	del sampleFrequency['strain2']
	del sampleFrequency['strain3']
elif strain_1v2 >= args.StrainThreshold:
	print("Printing Strains 1 and 3")
	sampleFrequency['strain1'] = sampleFrequency['strain1'] + sampleFrequency['strain2']
	del sampleFrequency['strain2']
elif strain_1v3 >= args.StrainThreshold:
	print("Printing Strains 1 and 2")
	sampleFrequency['strain1'] = sampleFrequency['strain1'] + sampleFrequency['strain3']
	del sampleFrequency['strain3']
elif strain_2v3 >= args.StrainThreshold:
	print("Printing Strains 1 and 2 (not 3). Strain2 frequency is updated")
	sampleFrequency['strain2'] = sampleFrequency['strain3'] + sampleFrequency['strain2']
	del sampleFrequency['strain3']
else:
	print("Printing Strains Kollo Allesta")
#elif strain_1v2 >= args.StrainThreshold and strain_1v3 >= args.StrainThreshold
#	print("Only one strain detected at %s% threshold")
#	system.exit("Strain composition error. Exit")


output = open("strainRef.fasta","w")
i = 0
non_empty_sequence = 0
problem = 0
for gene in referenceSeq:
	geneID = gene
	sequence = "".join(referenceSeq[gene])
	if sequence == "" or len(sequence) < 200:
		problem = problem + 1
		continue
	else:
		non_empty_sequence += 1
	if i == 0:
		output.write(">"+geneID+"\n"+sequence)
	else:
		output.write("\n>"+geneID+"\n"+sequence)
	i = i + 1 #incrementing sequence count
output.close()

OutputSampleName = ""
if args.SampleName != "":
	OutputSampleName = args.SampleName + "_"

output = open(OutputSampleName + "strain1.fasta","w")
i = 0
non_empty_sequence = 0
problem = 0
for gene in strain1:
	geneID = "strain1&" + gene
	sequence = "".join(strain1[gene])
	if sequence == "" or len(sequence) < 200:
		problem = problem + 1
		continue
	else:
		non_empty_sequence += 1
	if i == 0:
		output.write(">"+geneID+"\n"+sequence)
	else:
		output.write("\n>"+geneID+"\n"+sequence)
	i = i + 1 #incrementing sequence count
output.close()
print("Written %s sequences and skipped %s empty sequences" %(non_empty_sequence,problem))
if (strain_1v2 >= args.StrainThreshold and strain_1v3 >= args.StrainThreshold and strain_2v3 >= args.StrainThreshold):
	print("Only one strain detected at %s%% threshold"%(args.StrainThreshold))
	#system.exit("Strain composition error. Printed only 1 strain. Quit...")


if strain_1v2 < args.StrainThreshold:
	output = open(OutputSampleName + "strain2.fasta","w")
	i = 0
	non_empty_sequence = 0
	problem = 0
	for gene in strain2:
		geneID = "strain2&" + gene
		sequence = "".join(strain2[gene])
		if sequence == "" or len(sequence) < 200:
			problem = problem + 1
			continue
		else:
			non_empty_sequence += 1
		if i == 0:
			output.write(">"+geneID+"\n"+sequence)
		else:
			output.write("\n>"+geneID+"\n"+sequence)
		i = i + 1 #incrementing sequence count
		
	output.close()
	print("Strain2: Written %s sequences and skipped %s empty sequences" %(non_empty_sequence,problem))


if (strain_1v3 < args.StrainThreshold and strain_2v3 < args.StrainThreshold):
	output = open(OutputSampleName + "strain3.fasta","w")
	i = 0
	non_empty_sequence = 0
	problem = 0
	for gene in strain3:
		geneID = "strain3&" + gene
		sequence = "".join(strain3[gene])
		if sequence == "" or len(sequence) < 200:
			problem = problem + 1
			continue
		else:
			non_empty_sequence += 1
		if i == 0:
			output.write(">"+geneID+"\n"+sequence)
		else:
			output.write("\n>"+geneID+"\n"+sequence)
		i = i + 1 #incrementing sequence count

	output.close()

	print("Strain3: Written %s sequences and skipped %s empty sequences" %(non_empty_sequence,problem))	

if os.path.exists("./gene_files/"):
	print("gene_file/ folder exists. Placing reconstructed genes/strains there...")
else:
	os.mkdir("./gene_files/")
	print("Created ./gene_files/ folder. Placing reconstructed genes/strains there...")


#Output a concatenated version for everything
concatRef = ""
concatStrain1 = ""
concatStrain2 = ""
concatStrain3 = ""

for gene in referenceSeq:

	i = 0#First line flag

	if gene in geneMap:
		outputGeneFile = "gene_files/" + geneMap[gene]
		output = open(outputGeneFile,"w")
	else:
		outputGeneFile = "gene_files/error.fa"
		output = open(outputGeneFile,"a")


	geneID = gene
	sequence = "".join(referenceSeq[gene])
	if sequence == "" or len(sequence) < 200:
		problem = problem + 1
		continue
	else:
		non_empty_sequence += 1
	if i == 0:
		output.write(">"+geneID+"\n"+sequence)
		i = i + 1
	else:
		output.write("\n>"+geneID+"\n"+sequence)
	concatRef = concatRef + sequence

	geneID = "strain1&" + gene
	sequence = "".join(strain1[gene])
	if sequence == "" or len(sequence) < 200:
		problem = problem + 1
		continue
	else:
		non_empty_sequence += 1
	if i == 0:
		output.write(">"+geneID+"\n"+sequence)
		i = i + 1
	else:
		output.write("\n>"+geneID+"\n"+sequence)
	concatStrain1 = concatStrain1 + sequence

	if strain_1v2 < args.StrainThreshold:
		geneID = "strain2&" + gene
		sequence = "".join(strain2[gene])
		if sequence == "" or len(sequence) < 200:
			problem = problem + 1
			continue
		else:
			non_empty_sequence += 1
		if i == 0:
			output.write(">"+geneID+"\n"+sequence)
			i = i + 1
		else:
			output.write("\n>"+geneID+"\n"+sequence)
		concatStrain2 = concatStrain2 + sequence

	if (strain_1v3 < args.StrainThreshold and strain_2v3 < args.StrainThreshold):
		geneID = "strain3&" + gene
		sequence = "".join(strain3[gene])
		if sequence == "" or len(sequence) < 200:
			problem = problem + 1
			continue
		else:
			non_empty_sequence += 1
		if i == 0:
			output.write(">"+geneID+"\n"+sequence)
			i = i + 1
		else:
			output.write("\n>"+geneID+"\n"+sequence)
		concatStrain3 = concatStrain3 + sequence

	output.close()


with open(OutputSampleName + "ConcatenatedStrains.fasta","w") as concatStrain_file:
	#concatStrain_file.write(">ConcatReference\n")
	#concatStrain_file.write(concatRef + "\n")
	concatStrain_file.write(">" + OutputSampleName + "ConcatStrain1\n")
	concatStrain_file.write(concatStrain1 + "\n")
	if strain_1v2 < args.StrainThreshold:
		concatStrain_file.write(">" + OutputSampleName + "ConcatStrain2\n")
		concatStrain_file.write(concatStrain2 + "\n")
	if (strain_1v3 < args.StrainThreshold and strain_2v3 < args.StrainThreshold):
		concatStrain_file.write(">" + OutputSampleName + "ConcatStrain3\n")
		concatStrain_file.write(concatStrain3 + "\n")





with open("newkey_" + OutputSampleName + ".txt","w") as newkeyfile:
	for strain in sampleFrequency:
		newkeyfile.write(strain +"\t"+str(sampleFrequency[strain]) + "\n")
	newkeyfile.close()


TotalReferenceChanges = comitted_changes_count + N_comitted_changes_count
TotalStrain1 = strain1Changes + N_strain1Changes
TotalStrain2 = strain2Changes + N_strain2Changes
TotalStrain3 = strain3Changes + N_strain3Changes
OutputLogName = args.SampleName + "_GenotypingStatistics.txt"
with open(OutputLogName,"w") as writeStat:
	writeStat.write("Pairwise strains by variants:\t%s\t%s\t%s\n"%(strain_1v2,strain_1v3,strain_1v3))
	writeStat.write("Genes with changes %s\n"%(len(ReferenceChange)))
	writeStat.write("Positions changes %s\n"%(variantCounter))
	writeStat.write("Processed peaks %s\n"%(peakCount))
	writeStat.write("Genotyped reference changes: %s\n"%(TotalReferenceChanges))
	writeStat.write("\tNoneN Changes: %s\n"%(comitted_changes_count))
	writeStat.write("\tN Changes: %s\n"%(N_comitted_changes_count))
	writeStat.write("\t\tAll Ns Changes: %s\n"%(comitted_allN_count))
	writeStat.write("Singleton None-Ref variants:  %s\n"%(singletonNonRefAlleleCounter))
	writeStat.write("Total genotype variants: %s\n"%(TotalReferenceChanges + singletonNonRefAlleleCounter))	
	writeStat.write("Strain1: %s changes. Ns: %s\tnon-N %s\n"%(TotalStrain1,N_strain1Changes,strain1Changes))
	writeStat.write("Strain2: %s changes. Ns: %s\tnon-N %s\n"%(TotalStrain2,N_strain2Changes,strain2Changes))
	writeStat.write("Strain3: %s changes. Ns: %s\tnon-N %s\n"%(TotalStrain3,N_strain3Changes,strain3Changes))
	writeStat.write("Final Genotyped strains: %s\n"%(len(sampleFrequency)))

with open("GenotypedStrainsCount.txt","w") as writeStrainCount:
	writeStrainCount.write("%s\n"%(len(sampleFrequency)))


print("Strain1 has %s changes, and N comitted_changes_count %s" % (strain1Changes, N_strain1Changes))
print("Strain2 has %s changes, and N comitted_changes_count %s" % (strain2Changes, N_strain2Changes))
print("Strain3 has %s changes, and N comitted_changes_count %s" % (strain3Changes, N_strain3Changes))
print("Genotyped strains: %s\n"%(len(sampleFrequency)))
print("Reference comitted_changes_count %s, and N comitted_changes_count %s" %(comitted_changes_count,N_comitted_changes_count))
print("looped through %s variants in total" %(variantCounter))
#print("From: %s encounters.   Resolved a total of %s Ns and kept %s"%(counterEncounter, globalNctr,stayN))
#print("Ns from genes %s from those Ns from positions %s from those Ns from peaks %s"%(geneN,posN,peakN))
#print("peak1 (strain3) Ns %s peak2 (strain2) Ns %s peak3 (strain1) Ns %s peak4 (strain2,3) Ns %s peak5 (strain1,3) Ns %s peak6 (strain1,2) Ns %s" %(peak1NCounter,peak2NCounter,peak3NCounter,peak4NCounter,peak5NCounter,peak6NCounter))
