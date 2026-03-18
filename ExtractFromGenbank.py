#Date: 6th  of February, 2024
#Location: Raleigh, NC (USA)
#Purpose: Extract CDS fasta for core genome in CoreCruncher format
#Input: Genbank file
#output: Fasta file

#Changes
#Date: 7th of Fberuary, 2025 Location: Raleigh, NC (USA)
#Output Negative strand in this script and purge the ParseStrand script. 
#Workflow have changed since late 2023 and time to minimise and consolidate smaller scripts
#added new option for strand output
#Change feature to Gene instead of CDS

#Date: 14th of August, 2025 Location: Blacksburg, VA (USA)
#Forked file for MetaStrainer release. Now extract GenBank Style Fasta in addition to non-flank full genome gene features fasta.
#Maintaining CoreCruncher style of appending gene coordinates to gene feature in Fasta IDs.
#Print coordinates file here (that was done in shell in previous workflow). 


#Date: 16th of March, 2026 Location: Blacksburg, VA (USA)
#Excluding some MGE features to reduce hypervariable genes adding noise and inflating strain count
#TODO: make it user customisable later

import os
import argparse
import sys


from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input", required=True, help="Input genbank file")
parser.add_argument("-o","--output", required=True, help="Output Fasta file")
parser.add_argument("-s","--strand", required=True, help="Strand information")
parser.add_argument("-a","--assembly", required=True, help="Output Genome chromosomes/contigs fasta file")
parser.add_argument("-c","--coordinates", default="coord.lst", help="Coordinates")

args = parser.parse_args()

if not (os.path.exists(args.input)):
	errorMessage = "GenBank file " + args.input+ " for -i/--input does not exist"
	sys.exit(errorMessage)

#Debuging GenBank files for discrepancies
geneFeaturesCount = 0
CDSFeaturesCount = 0
writtenGeneCount = 0


sequence = ""
sequence_name = ""

writeOut=open(args.output,"w")
writeStrand=open(args.strand,"w")
writeCoord=open(args.coordinates,"w")

records = SeqIO.parse(args.input, "genbank")
writeAssemblyFasta=open(args.assembly,"w")
contigCount = SeqIO.write(records, writeAssemblyFasta, "fasta")
print("Printed %s contigs/chromosomes to %s."%(contigCount,args.assembly))
writeAssemblyFasta.close()


#Mar 16 2026 block

#track skipped Pseudo genes
skipPseudoCount = 0

#track skipped potential
skipMGECount = 0

#Exclude list
#TODO: make it customisable by user
excludeMGESet = {
"transposase",
"integrase",
"phage",
"capsid",
"terminase"
} 

CDSProductDictionary = {}

#Mar 16 2026


for record in SeqIO.parse(args.input, "genbank"):
	#print("Contig %s "%(record.id))
	#>GCA_002088735.1_ASM208873v1_genomic.prot&NAHB01000058.1_44868_46430
	#if assembly_name == "":
	#	for link in record.dbxrefs:
	#		print("Links")
	#		print(link)
	#	assembly_name = "1"

	#Mar 16 2026
	#Need another loop to store CDS and products to skip MGEs instead of writing everything
	#Check that dictionary first
	for feature in record.features:
		contig_name = record.id
		start = 0
		end = 99999999999999

		if (feature.type == "CDS"):
			CDSFeaturesCount = CDSFeaturesCount + 1

			#The following block is copied as is from the gene feature loop to maintain consistency
			start = int(feature.location.start)
			end = int(feature.location.end)

			start_id = start + 1
			sequence_name = ">" + record.id + "_" + str(start_id) + "_" + str(end)
			productText = feature.qualifiers.get("product",[""])[0].lower()#Make it lower case for text comparisons

			CDSProductDictionary[sequence_name] = productText




	for feature in record.features:
		contig_name = record.id
		start = 0
		end = 99999999999999



		if (feature.type == "gene"):
			geneFeaturesCount = geneFeaturesCount + 1#if you put it above, it will count only number of contigs
			
			
			start = int(feature.location.start)
			end = int(feature.location.end)

			
			#start is 0-based so increment by 1 for output to help with downstream parsing
			start_id = start + 1
			sequence_name = ">" + record.id + "_" + str(start_id) + "_" + str(end)


			#Mar 16 2026
			############
			#Skip Pseudo Genes
			if "pseudo" in feature.qualifiers:
				skipPseudoCount = skipPseudoCount + 1
				continue

			#Skip MGE
			productName = CDSProductDictionary.get(sequence_name,"")

			found = False
			for MGE in excludeMGESet:
				if MGE in productName:
					found = True
					break
			if found:
				skipMGECount = skipMGECount + 1
				continue

			############
						
			if (feature.location.strand == -1):
				sequence = str(record.seq[start:end].reverse_complement())
				writeStrand.write(record.id + "_" + str(start_id) + "_" + str(end) + "\n")
			else:
				sequence = str(record.seq[start:end])

			writeOut.write(sequence_name + "\n" + sequence + "\n")
			writeCoord.write(record.id +"\t" +str(start_id)+"\t"+str(end)+"\n")

			#Mar 16 2026
			############
			#Counting written features
			writtenGeneCount = writtenGeneCount + 1
			############
			
writeOut.close()
writeStrand.close()
writeCoord.close()

#Mar 16 2026
############
#Loggin and debugging
print("GenBank gene features: %s"%(geneFeaturesCount))
print("GenBank CDS features: %s"%(CDSFeaturesCount))
print("Skipped Pseudo Genes: %s"%(skipPseudoCount))
print("Skipped MGE Genes: %s"%(skipMGECount))
print("Written Mapping Reference Genes: %s"%(writtenGeneCount))