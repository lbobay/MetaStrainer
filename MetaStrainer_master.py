#Date: 29th  of Jul, 2025
#Location: Blacksburg, VA (USA)
#Purpose: MetaStrainer central file

import sys
import os
import argparse

parser = argparse.ArgumentParser(
	description='''
	MetaStrainer requires read pair fastq input file, a reference GeneBank file for alignment.
	User is required to provide a region of DNA to be added upstream and downsteam of genes of interest.
	It is recommended to be equal to the read length
	SampleName is recommended especially when batch processing samples.''',
	formatter_class=argparse.RawDescriptionHelpFormatter
	)

#Required
parser.add_argument("-1","--firstfile", required=True, help="Forward or first of pair")
parser.add_argument("-2","--secondfile", required=True, help="Reverse or second of pair")
parser.add_argument("-r","--reference", required=True, help="genbank reference")
parser.add_argument("-f","--flankregion", required=True, type = int, help="flank region to be added up stream and downsteam of gene of interest.")
parser.add_argument("-o","--output", required=True, help="output foldername")
#Optional
parser.add_argument("-S","--SampleName", help="Sample name. Can be helpful when batch running several samples. Default is output folder name")
parser.add_argument("-s","--StrainThreshold",default=99.5,type=float,help="Threshold for defining a strain for merging")
parser.add_argument("--threads", default=1, type = int, help="Number of threads to be used")
parser.add_argument("--seed", default=12345, type = int, help="MetaStrainer seed")
parser.add_argument("--force", default=False, action='store_true',help="Force overwrite")

#Preprocess Options
#no. reads singleton?
#std. dev reads coverage

#MetaStrainer Options
#MetaStrainer seed

#Genotype Options
#strain threshold
#genotype confidence?

#Optional postprocessing

args = parser.parse_args()

#Check for programs (bowtie2)
#Check for libraries (numpy emcee)
try:
    import emcee
    import numpy
    from Bio import SeqIO
except ImportError:
    sys.exit("Insall packages emcee,numpy and Bio before proceeding")


if args.SampleName is None:
	args.SampleName = os.path.basename(args.output)
	

#Check for folder existence
if (os.path.isdir(args.output)):
	if (args.force):
		print("Folder %s Exists. Overwriting."%(args.output))
		#remove folder
		#make folder
	else:
		sys.exit("%s Exists. Rerun with --force or choose another folder/sample name for output.")
else:
	try:
		os.mkdir(args.output)
	except:
		sys.exit("Error creating output folder")

ExecutedCommands = open(args.output +"/ExecutedCommands.log","w")

#Copied from CoreCruncher
loc = ""
for stuff in sys.argv:
	if "MetaStrainer_master.py" in stuff:
		loc = stuff.split("MetaStrainer_master.py")[0]
		print("Location= ",loc)


#Ceate Logs folder
try:
	os.mkdir(args.output+"/Logs/")
except:
	if (os.path.isdir(args.output+"/Logs")):
		print("Logs folder already exists")
	else:
		sys.exit("Error creating Logs folder")
LogsFolder = args.output+"/Logs/"



#Prepare Mapping reference
#Create folder
try:
	os.mkdir(args.output+"/Reference")
except:
	if (os.path.isdir(args.output+"/Reference")):
		print("Reference folder already exists")
	else:
		sys.exit("Error creating Reference folder")
ReferenceFolder = args.output+"/Reference/"

#Get Reference name
if not os.path.exists(args.reference):
    sys.exit("Reference Genbank file does not exist%s" %(args.reference))
RefFileName = os.path.basename(args.reference)
#ExtraCDSfromGenbank
print("Extracting Genome sequence and gene features from %s"%(RefFileName))
RefDBName = os.path.splitext(RefFileName)[0]
#RefDBNameFasta = RefDBName + ".fasta"
RefDBNameGenomeFasta = args.output+"/Reference/" +RefDBName + "_FullGenome.fasta"#RefSeq Syle
#Gene Features with GeneBank coordinate (used downstream and in genotyping)
RefDBNameFasta = args.output+"/Reference/" +RefDBName + ".fasta"
#Gene Features with flanking 5' and 3'regions (used from alignment only)
RefDBNameFlankFasta = args.output+"/Reference/" +RefDBName + "_" + str(args.flankregion) +".fasta"
#Trimming information which may be modified from flank range if gene feature is at contig boundary
RefDBNameRangeTrimming = args.output+"/Reference/" +RefDBName + "_flankTrim_" + str(args.flankregion) +".txt"
#GeneID Mapping between standard and flank modified fasta IDs
RefDBNameRangeMapping = args.output+"/Reference/" +RefDBName + "_family2gene_mapping_" + str(args.flankregion) +".txt"
#Negative strand file
RefDBNameNegstrand = args.output+"/Reference/" +RefDBName + "_negstrand.lst"

#os.system("python %sExtractCDSFromGenbank_CoreCruncher.py -i %s -o %s/Reference/%s -s %s/Reference/%s_negstrand.lst"%(loc,args.reference,args.output,RefDBNameFasta,args.output,RefDBName))
retcode = os.system("python %sExtractFromGenbank.py -i %s -o %s -s %s -a %s"%(loc,args.reference,RefDBNameFasta,RefDBNameNegstrand,RefDBNameGenomeFasta))
if (retcode != 0):
	sys.exit("Error Extracting features from Genbank file")

ExecutedCommands.write("python %sExtractFromGenbank.py -i %s -o %s -s %s -a %s\n\n"%(loc,args.reference,RefDBNameFasta,RefDBNameNegstrand,RefDBNameGenomeFasta))
ExecutedCommands.flush()
#os.system("echo -e \"baba blacksheep\\n\\n\\n\\n\\nHave you any wool %s\""%(RefDBName))

retcode = os.system("python %sAddFlankToGene.py -f %s -i coord.lst -o %s -r %s -s %s -m %s -t %s"%(loc,RefDBNameGenomeFasta,RefDBNameFlankFasta,args.flankregion,RefDBNameNegstrand,RefDBNameRangeMapping,RefDBNameRangeTrimming))
if (retcode != 0):
	sys.exit("Error generating modified reference file. Failed to add flanking regions.")
ExecutedCommands.write("python %sAddFlankToGene.py -f %s -i coord.lst -o %s -r %s -s %s -m %s -t %s\n\n"%(loc,RefDBNameGenomeFasta,RefDBNameFlankFasta,args.flankregion,RefDBNameNegstrand,RefDBNameRangeMapping,RefDBNameRangeTrimming))
ExecutedCommands.flush()

retcode = os.system("bowtie2-build %s %s%s 1>%sBuildingBowtie.log"%(RefDBNameFlankFasta,ReferenceFolder,RefDBName,LogsFolder))
if (retcode != 0):
	sys.exit("Error building bowtie2 index")
ExecutedCommands.write("bowtie2-build %s %s%s 1>%sBuildingBowtie.log\n\n"%(RefDBNameFlankFasta,ReferenceFolder,RefDBName,LogsFolder))
ExecutedCommands.flush()
os.system("date")


#Map
#Create folder
try:
	os.mkdir(args.output+"/Mapping")
except:
	if (os.path.isdir(args.output+"/Mapping")):
		print("Mapping folder already exists")
	else:
		sys.exit("Error creating Mapping folder")
MappingFolder=args.output+"/Mapping/"
MappingName=MappingFolder+args.SampleName

os.system("bowtie2 -p %s -x %s%s --very-sensitive --no-unal -1 %s -2 %s 2>%sMapping.log | samtools view -@ %s -Sb - | samtools sort -@ %s -o %s_sort_sen_nounal.bam "
	%(args.threads,ReferenceFolder,RefDBName,args.firstfile,args.secondfile,LogsFolder,args.threads,args.threads,MappingName))
if (retcode != 0):
	sys.exit("Error mapping fastq files to reference")
ExecutedCommands.write("bowtie2 -p %s -x %s%s --very-sensitive --no-unal -1 %s -2 %s 2>%sMapping.log | samtools view -@ %s -Sb - | samtools sort -@ %s -o %s_sort_sen_nounal.bam\n\n"
	%(args.threads,ReferenceFolder,RefDBName,args.firstfile,args.secondfile,LogsFolder,args.threads,args.threads,MappingName))
ExecutedCommands.flush()


os.system("samtools view -@ %s -h %s_sort_sen_nounal.bam > %s_sort_sen_nounal.sam"%(args.threads,MappingName,MappingName))
if (retcode != 0):
	sys.exit("Error converting BAM to SAM")
ExecutedCommands.write("samtools view -@ %s -h %s_sort_sen_nounal.bam > %s_sort_sen_nounal.sam\n\n"%(args.threads,MappingName,MappingName))
ExecutedCommands.flush()





os.system("date")


#Preprocess
try:
	os.mkdir(args.output+"/Preprocess")
except:
	if (os.path.isdir(args.output+"/Preprocess")):
		print("Preprocess folder already exists")
	else:
		sys.exit("Error creating Preprocess folder")
PreprocessFolder=args.output+"/Preprocess/"

print("Processing Read Pairs")
os.system("python %sPairingReads.py -s %s_sort_sen_nounal.sam -m %s -r %s -t %s -n %s > %sPairingReads.log"
	%(loc,MappingName,RefDBNameRangeMapping,RefDBNameFasta,RefDBNameRangeTrimming,RefDBNameNegstrand,LogsFolder))
if (retcode != 0):
	sys.exit("Error generating variants list. Check SAM file or any of the reference database files.")
ExecutedCommands.write("python %sPairingReads.py -s %s_sort_sen_nounal.sam -m %s -r %s -t %s -n %s > %sPairingReads.log\n\n"
	%(loc,MappingName,RefDBNameRangeMapping,RefDBNameFasta,RefDBNameRangeTrimming,RefDBNameNegstrand,LogsFolder))
ExecutedCommands.flush()
os.system("date")

print("Generating Linkage groups")
os.system("python %sLinkageGroups.py> %sLinkgageGroups.log"%(loc,LogsFolder))
if (retcode != 0):
	sys.exit("Error generating linkage groups.")
ExecutedCommands.write("python %sLinkageGroups.py > %sLinkgageGroups.log\n\n"%(loc,LogsFolder))
ExecutedCommands.flush()

print("Generating linkgage groups variant pairs")
os.system("python %sMakeContigPairs.py > %sContigPairs.log"%(loc,LogsFolder))
if (retcode != 0):
	sys.exit("Error generating contig pairs.")
ExecutedCommands.write("python %sMakeContigPairs.py > %sContigPairs.log\n\n"%(loc,LogsFolder))
ExecutedCommands.flush()


#MetaStrainer
print("Running core MetaStrainer chain")
os.system("python %sCoreChain.py -S %s -s contig_singles.txt -p contig_pairs.txt --seed %s > %sMetaStrainerChain.log"%(loc,args.SampleName,args.seed,LogsFolder))
if (retcode != 0):
	sys.exit("Error running MetaStrainer chain")
ExecutedCommands.write("python %sCoreChain.py -S %s -s contig_singles.txt -p contig_pairs.txt --seed %s> %sMetaStrainerChain.log\n\n"%(loc,args.SampleName,args.seed,LogsFolder))
ExecutedCommands.flush()


#Genotyping
print("Producing strain genotypes")
os.system("python %sGenotypingStrains.py -i genotypes_%s.txt -o Genotype -S %s  -f key_genotypes_%s.txt -g %s -r %s -l linkage_groups.txt > %sGenotype.Log"
	%(loc,args.SampleName,args.SampleName,args.SampleName,RefDBNameRangeMapping,RefDBNameFasta,LogsFolder))
if (retcode != 0):
	sys.exit("Error Genotyping Strains")
else:
	print("Genotyping successful.")
	with open("GenotypedStrainsCount.txt","r") as strainCountFile:
		strainCount = strainCountFile.readline()
		if (strainCount):
			print("Strains genotyped: %s"%(strainCount.strip()))
os.system("date")

ExecutedCommands.write("python %sGenotypingStrains.py -i genotypes_%s.txt -o Genotype -S %s  -f key_genotypes_%s.txt -g %s -r %s -l linkage_groups.txt > %sGenotype.Log\n\n"
	%(loc,args.SampleName,args.SampleName,args.SampleName,RefDBNameRangeMapping,RefDBNameFasta,LogsFolder))
ExecutedCommands.flush()

#Optional




