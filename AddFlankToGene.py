#Date: 4th of April, 2022
#Location: Blacksburg, VA (USA)
#Purpose: Modify core by adding flanking regions on 5' and 3' 
#Rationale: 
#Input: fasta file, coordiates and range
#output: modified fasta
#TODO: add log with boundary modifications 
#Bugfix Nov 30, 2022, Adjust locii exceeding boundary limits (e.g. start less than x bp from ATG)
#Bugfix Jan 4, 2023, Fix 5' pos id to adjust for 1-based notation and match original core crunch coordinates.
#Bugfix Mar 17,2023 Output flankless corrdinate to flank_modified coordinate
#Bugfix Sep 8, 2023 Add negative strand information in advace
#Bugfix Sep 8, 2023 Enact modification from Mar 17 along with fix for start 1

#Feb 7, 2025 (Raleigh, NC) Fork file to remove dependance on Core Cruncher
#

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq #reverse complement sequence
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f","--fasta", required=True, help="fasta input file")
parser.add_argument("-i","--input", required=True, help="Input file should be a tabular file coordinates.")
parser.add_argument("-o","--output", required=True, help="gene features having specified flanking region in fasta output file")
parser.add_argument("-r","--range", required=True, type = int, help="range of flanking region to be included")
parser.add_argument("-s","--strand", required=True, help="Negative strand information to reverse complement to match core cruncher output and alignments downstream")
parser.add_argument("-m","--mapping", required=True, help="Output original locus id to locus with flanking region id mapping with maping to core cruncher fam id")
parser.add_argument("-t","--trimming", required=True, help="Output trimming lengths information for flanking regions 5' and 3' ")


args = parser.parse_args()

#example
#python ../Code/BeeGutProject/getFastaRecords_subrange.py -i core_snod_NAHB01_coordinates.txt -f NAHB01.fasta -r 100 -o core_snod_NAHB01_100.fasta 
#python ../Code/BeeGutProject/getFastaRecords_subrange.py -i core_gapis_NAHB01_coordinates.txt -f NAHB01.fasta -r 200 -o core_Gilliamella_apicola_MZNG01_200.fasta -m core_Gilliamella_apicola_MZNG01_200_mapping.txt
#python ../../Code/BeeGutProject/getFastaRecords_subrange.py -i ../core_gili_coordinates.txt -f ../MZNG01.fasta -r 150 -o core_Gilliamella_apicola_MZNG01_150.fasta -m core_Gilliamella_apicola_MZNG01_150_mapping.txt -s proper_negative_strand_locii_0 -a /Users/hazemsharaf/Projects/Bee_variants_project/core_genome/Gapicola_MZNG_family_gene_0.txt -t Gapicola_MZNG_flankTrim_150.txt
#python ../../Code/BeeGutProject/getFastaRecords_subrange.py -i ../core_gili_coordinates.txt -f ../MZNG01.fasta -r 0 -o core_Gilliamella_apicola_MZNG01_0.fasta -m core_Gilliamella_apicola_MZNG01_0_mapping.txt -s proper_negative_strand_locii_0 -a /Users/hazemsharaf/Projects/Bee_variants_project/core_genome/Gapicola_MZNG_family_gene_0.txt -t Gapicola_MZNG_flankTrim_0.txt



inputFastaDict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
#SeqIO.write needs a list iterator
#for seq in inputFastaDict:
	#print(seq + " has a length of " + str(len(inputFastaDict[seq].seq)))

outputFastaList = []

negative_strand_set = set()
i = 0
with open(args.strand) as negative_strand_file:
	lines = negative_strand_file.readline()
	while lines:
		locus = lines.strip("\n")
		i = i + 1
		negative_strand_set.add(locus)
		lines = negative_strand_file.readline()
print("Added %s locii to negative strand set" %(i),end="\n")


#/Users/hazemsharaf/Projects/Bee_variants_project/Gapicola_sim150_3st_1.5g_2/genotype_comparison

mapping_file = open(args.mapping,"w")
trimming_file = open(args.trimming,"w")
with open(args.input) as core_genes:
	lines = core_genes.readline()
	i = 0
	while lines:
		i = i + 1
		trim_start = 0
		trim_end = 0

		readID, start, end = lines.strip("\n").split("\t")
		old_recordID =readID + "_" + start + "_" + end
		
		#TODO: inform user if several regions on the same sequence have been added before
		if (int(start) - args.range - 1) < 0:
			print(readID + " with start " + str(start) + " is set to 1")
			trim_start = (int(start)  - 1) 
			start = 0
		else:
			start = int(start) - args.range - 1
			trim_start = args.range
		if (int(end) + args.range >= len(inputFastaDict[readID].seq)):
			print(readID + " with end " + str(end) + " is set to " + str(len(inputFastaDict[readID].seq)))
			trim_end = len(inputFastaDict[readID].seq) - int(end)
			end = len(inputFastaDict[readID].seq)
		else:
			end = int(end) + args.range
			trim_end = args.range
		new_recordID = readID + "_" + str(int(start) + 1) + "_" + str(int(end))
		new_sequence = inputFastaDict[readID].seq[start:end]

		if old_recordID in negative_strand_set:
			new_sequence = new_sequence.reverse_complement()

		outputFastaList.append(SeqRecord(new_sequence,new_recordID,'',''))

		#Keeping backward compatibility with fam to gene mapping for preprocessing
		mapping_file.write("gene" + str(i) + "\t" + old_recordID + "\t" + new_recordID + "\n")
		trimming_file.write(new_recordID + "\t" + str(trim_start) + "\t" + str(trim_end) + "\n")
		

		lines = core_genes.readline()

SeqIO.write(outputFastaList,args.output,'fasta')

print("Sequences written: %s" % i)