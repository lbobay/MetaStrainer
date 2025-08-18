import os
import sys
import numpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s","--samfile", required=True, help="sample SAM File to parse")
parser.add_argument("-m","--mapfile", required=True, help="Family Gene and flanked gene Mapping file")
parser.add_argument("-r","--reference", required=True, help="Reference Genome WITHOUT FLANKS. This is new in Louis-Marie implementation this time")
parser.add_argument("-d","--depth", type=float,default=1.5, help="Depth Standard Deviations for filtering variants (not alleles)")
parser.add_argument("-l","--logfile", default="processing.log", help="Log removed Gene/Pos")
parser.add_argument("--stat", default="PreprocessingStat.log", help="Summary of filtered alleles")
#Bring them back, needed to properly trim genes at contigs 3' or 5' ends
parser.add_argument("-t","--trimflankfile",required=True, help="trim information")
parser.add_argument("-n","--negstrand",required=True, help="complementary strand information")
args = parser.parse_args()

logfile = open(args.logfile,"w",buffering=100)
summaryfile = open(args.stat,"w",buffering=100)

locus_trim_hash = {}
#with open("/Users/hazemsharaf/Projects/Bee_variants_project/core_genome/Gapicola_MZNG_flankTrim_150.txt") as input_file:
with open(args.trimflankfile) as input_file:
	lines = input_file.readline()
	while lines:
		lineContents = lines.strip("\n").split("\t")
		locus = lineContents[0]
		trim_5 = lineContents[1]
		trim_3 = lineContents[2]
		locus_trim_hash[locus] = trim_5 + "&" + trim_3
		lines = input_file.readline()

neg_hash = []
#with open("/Users/hazemsharaf/Projects/Bee_variants_project/core_genome/Gapicola_MZNG_flankTrim_150.txt") as input_file:
with open(args.negstrand) as input_file:
	lines = input_file.readline()
	while lines:
		lineContents = lines.strip("\n")
		if lineContents not in neg_hash:
			neg_hash.append(lineContents)

		lines = input_file.readline()


flags={}
dico={}
counter=0
tot=0
#f=open("../Mapped/M1_S1_core_Smarce_MEDA01_150_sort_sen_nounal.sam")
#f=open("../../Mapped/Sample5_core_Gapicola_MZNG01_150_3st_2g_col2_sort_sen_nounal.sam")
#f=open("../Mapped/Sample4_core_Gapicola_MZNG01_150_3st_2g_col2_sort_sen_nounal.sam")
#f=open("Sample4_core_Gapicola_MZNG01_150_3st_2g_col2_sort_sen_nounal_02.sam","r")
#f=open("../Sample2_core_Gapicola_MZNG01_150_3st_2g_col2_sort_sen_nounal.sam","r")
f=open(args.samfile,"r")
for l in f:
	if l[0] != "@":
		counter+=1
		a=l.strip("\n").split("\t")
		read = a[0]
		#print(a[:10])
		cigar=a[5]
		phred=a[10]
		if cigar == "*":
			tot+=1
		if a[2] != "*" and "M" in a[5]:
			if a[1] in flags:
				flags[a[1]]+=1
			else:
				flags[a[1]]=1
			if read not in dico:
				dico[read] =[[a[1], a[2],a[9],int( a[3] ),cigar,phred]]
			else:
				dico[read].append([a[1], a[2],a[9] ,int( a[3] ),cigar,phred])	
		#if counter==100000:
			#break		

f.close()


print("Cigar with *",tot)

counter=0
tot=0
for read in dico:	
	if len(dico[read])==2:
		cigar1,cigar2=dico[read][0][-2],dico[read][1][-2]
		tot+=1
		if cigar1 != "*" and cigar2 != "*":
			counter+=1

print("tot=",tot)
print("counter=",counter)
#un,deux,more=0,0,0
#for read in dico:
#	if len(dico[read]) == 1:
#		un+=1
#	elif len(dico[read]) == 2:
#		deux+=1
#	else:
#		more+=1


#print(un,deux,more)




ref={}
#f=open("../../ref.fa","r")
f=open(args.reference,"r")
for l in f:
	if l[0]==">":
		name = l.strip("\n").strip(">").split(" ")[0]
		ref[name]=""
	else:
		ref[name]+=l.strip("\n")

f.close()


shift={}
gene={}
rev={}
#f=open("../../../core_genome/core_Gilliamella_apicola_MZNG01_150_mapping.txt","r")
f=open(args.mapfile,"r")
for l in f:
	a=l.strip("\n").split("\t")
	gene[a[2]]=a[1]
	rev[a[1]] = a[2]
	pos1,pos2=int(a[1].split("_")[-2]),int(a[2].split("_")[-2])
	shift[a[1]] = abs(pos1-pos2)
	shift[a[2]] = abs(pos1-pos2)

	
	#May 9th, 2025 (Raleigh, NC)
	#Reverse flank trimming for negative strand genes
	#Needed to fix singleton allele counts and variant positions in genes on contig ends
	(trimming_3,trimming_5) = map(int,locus_trim_hash[a[2]].split("&"))

	#May 13th, 2025 (Blacksburg, VA)
	#Keep it a[1] and use the file file generated with the reference through the 
	if (a[1] in neg_hash):
		shift[a[1]] = trimming_5
		shift[a[2]] = trimming_5
	else:
		shift[a[2]] = abs(pos1-pos2)


	

f.close()





# Finds single polymorphisms

variants={}
for gene_id in ref:
	variants[gene_id]={}


for read in dico:
	for part in dico[read]:
		#print(part)
		name = part[1]
		seq=part[2]
		pos=int(part[3]) - 1
		cigar=part[-2]
		phred=part[-1]
		if "D" not in cigar and "I" not in cigar and gene[name] in ref: 
			rescale = pos - shift[name]
			tag=0
			i = 0
			adj=0
			if rescale<0:
				#print(part)
				#print("L=",L,"res=",rescale,"adj=",adj)
				adj = shift[name] - pos
				rescale = 0
				tag=1
			L=min(   [len(seq[adj:]),len(ref[gene[name]][rescale:])]  )
			#if tag ==1 and L < 40:
				#print(part)
				#print("L=",L,"res=",rescale,"adj=",adj)
				#print(">read")
				#print(seq[adj:])
				#print(">seq")
				#print(ref[gene[name]])
			score = 0.0
			while i < L:
				N1=seq[i+adj]
				N2 = ref[gene[name]][rescale+i]
				quality = ord(phred[i+adj]) - 33
				#print(N1,N2,quality)
				#Hazem: Set quality to 30
				if quality >= 30:
					try:
						variants[gene[name]][rescale+i].append(N1)							
					except KeyError:
						variants[gene[name]][rescale+i]=[N2,N1]
				if N1==N2:
					score+=1
				#print(N,N2)
				i+=1


#for name in variants:
#	print(len(variants[name].keys() ) ) 




alpha={}
alpha["A"]="y"
alpha["C"]="y"
alpha["G"]="y"
alpha["T"]="y"
	
#Hazem:sorting variants positions inside genes for easier debugging
variants = {key : dict(sorted(val.items(), key = lambda ele: ele[0]))
       for key, val in variants.items()}

#Will write file after removing anomalous coverage Jul 29th (Giza, Egypt)
#Write singletons down
#After meeting with Louis-Marie: Keep Track of singleton similar to Ref and not similar to Ref
file_singleton_ref=open("singles_singleton_ref.txt","w")
file_singleton_noref=open("singles_singleton_noref.txt","w",buffering=100)
salvage = 0
poly={}
singleton_allele_noref_count = 0
singleton_allele_noref_10reads_count = 0
for name in variants:
	poly[name]={}
	for pos in variants[name]:
		poly[name][pos]= list(set(variants[name][pos]))
		kick=[]
		if (name == "MZNG01000002.1_89086_93303" and (pos == 3939 or pos == 4215)):
			print("---------\nStarting Loop for problems")
			print("Reference is %s" % (ref[name][pos]))
			print("Position %s  (%s) has the following" %(pos,(int(pos)+151)))
			print("Poly %s"%(poly[name][pos]))
			print("Length Poly %s"%(len(poly[name][pos])))
			print("Variants length is %s" % (len(variants[name][pos][1:])))
			print("Variants no. items = %s " % len(set(variants[name][pos])))
			print("NoIndex Variants length %s" % (len(variants[name][pos])))
			print("Variants items = %s " % (set(variants[name][pos])))
			print("---------")
			
		if len(poly[name][pos]) > 1:
			#print(poly[name][pos])
			#print(variants[name][pos])
			tot=0.0
			train=[]
			nucleotides=[]
			tmp = variants[name][pos][1:]
			for N in poly[name][pos]:
				nb=tmp.count(N)
				if nb >= 2:																				# Minimum number of reads to call a variant
					train.append(nb)
					nucleotides.append(N)
					tot+=nb
				else:
					kick.append(N)

			#Hazem, move kick loop up and clean
			for N in kick:
				while N in poly[name][pos]:
					poly[name][pos].remove(N)
			kick = []

			if (name == "MZNG01000002.1_89086_93303" and (pos == 3939 or pos == 4215)):
				print("Done kicking...showing variants remaining")
				print(poly[name][pos])
				print(len(poly[name][pos]))
				print("See? Problem here")


			#if condition and keep singleton alleles
			if len(poly[name][pos]) == 1:
				if (poly[name][pos][0] != ref[name][pos]):
					#Needed the following or else it will be supurious
					nb=tmp.count(poly[name][pos][0])
					if (nb/tot == 1):
						file_singleton_noref.write(name + "\t" + str(pos) + "\t" + poly[name][pos][0] + "\t" +  str(nb/tot) + "\t" + str(tot) + "\t" + str(nb) + "\n")
						singleton_allele_noref_count = singleton_allele_noref_count + 1
						if tot > 10:
							singleton_allele_noref_10reads_count = singleton_allele_noref_10reads_count + 1
						#if (name == "MZNG01000002.1_89086_93303" and (pos == 3939 or pos == 4215)):
						#	print("Supposedly writing " + name + "\t" + str(pos) + "\t" + N + "\t" +  str(nb/tot) + "\t" + str(tot) + "\t" + str(nb))
					#else:
						#if (name == "MZNG01000002.1_89086_93303" and (pos == 3939 or pos == 4215)):
						#	print("\n***Did not write pos %s with NB %s and tot %s fake N %s and proper %s in %s ***\n"%(pos,nb,tot,N,poly[name][pos][0],poly[name][pos]))
				else:
					file_singleton_ref.write(name + "\t" + str(pos) + "\t" + N + "\t" +  str(nb/tot) + "\t" + str(tot) + "\t" + str(nb) + "\n")
				#Do not remove now
				#del poly[name][pos]

			#Hazem: Remove alleles <1% or more than >99%
			original_tot = tot
			for N in poly[name][pos]:
				nb=tmp.count(N)
				if nb >=2 and (nb/original_tot < 0.01 or nb/original_tot > 0.99):
					train.remove(nb)
					nucleotides.remove(N)
					tot-=nb
					kick.append(N)	


			#Hazem, kick loop after removing 1 %
			for N in kick:
				while N in poly[name][pos]:
					poly[name][pos].remove(N)
			i=0
			while i < len(train):
				nb=train[i]
				N = nucleotides[i]
				#Hazem: make sure single alleles are removed from singles but tracked in singleton output for later merge
				if N in alpha and len(train) > 1:
				#if N in alpha:
					#h.write(name + "\t" + str(pos) + "\t" + N + "\t" +  str(nb/tot) + "\t" + str(tot) + "\t" + str(nb) + "\n")
					pass
				i+=1					
			#Hazem: if only 1 allele exist but it is different than reference, keep track of it for reconstruction later on
			if len(poly[name][pos]) == 1:
				del poly[name][pos]
				logfile.write("%s\t%s\n"%(name,pos))
		elif len(poly[name][pos]) <= 1:
			if (len(poly[name][pos]) == 1 and poly[name][pos][0] != ref[name][pos]):
				salvage = salvage + 1
				#print that if salvage counter is more than 0
			del poly[name][pos]
			logfile.write("%s\t%s\n"%(name,pos))
		#Hazem: sometimes list are empty but position is not deleted
		if (pos in poly[name] and len(poly[name][pos]) <=1):
			#print("Deleting position in Name: ",name,"\tPos: ",pos,"\t Length: ",str(len(poly[name][pos])),"\t Content: ",poly[name][pos])
			del poly[name][pos]
			logfile.write("%s\t%s\n"%(name,pos))
	#variants[name]=""
			
print("Salvaged %s more positions"%(salvage))#if that is ever more than 0 then I will need to print them to file
summaryfile.write("Singletons\t%s\n"%(singleton_allele_noref_count))
summaryfile.write("Singletons (min. 10 reads)\t%s\n"%(singleton_allele_noref_10reads_count))
summaryfile.close()

file_singleton_ref.close()
file_singleton_noref.close()
#sys.exit()#Added that temporary since I just need to reprocess the singleton to print them correctly. Comment when running reprocessing a sample from scratch
#Filter anomalous depths (should come from command line)
depthSetting = args.depth
geneSet = {}
globalAlleleDepth = []
passingVariants = 0
filtered = 0
filteredAbove = 0
filteredBelow = 0
deletPosResu = set()
deleteCounter = 0
h = open("tmp_singles.txt","w")
for name in poly:
	if name not in geneSet:
		geneSet[name] = []

	for pos in poly[name]:
		#poly[name][pos]= list(set(variants[name][pos]))
		tot=0.0
		tmp = variants[name][pos][1:]
		for N in poly[name][pos]:
			nb=tmp.count(N)
			tot+=nb
		globalAlleleDepth.append(tot)
		geneSet[name].append(tot)

	deletePositions = []
	for pos in poly[name]:
		tot=0.0
		tmp = variants[name][pos][1:]
		for N in poly[name][pos]:
			nb=tmp.count(N)
			tot+=nb

		geneMeanDepth = numpy.mean(geneSet[name])
		geneDepthSD = numpy.std(geneSet[name])

		if tot >= (geneMeanDepth - geneDepthSD*depthSetting) and tot <= (geneMeanDepth + geneDepthSD*depthSetting):
			
			passingVariants = passingVariants + 1
		else:
			filtered = filtered + 1
			if tot < (geneMeanDepth - geneDepthSD*depthSetting):
				filteredBelow = filteredBelow + 1
			if tot > (geneMeanDepth + geneDepthSD*depthSetting):
				filteredAbove = filteredAbove + 1
			#Delete position
			deletePositions.append(pos)

	for pos in deletePositions:
		del poly[name][pos]
		logfile.write("Deleted\t%s\t%s\tDepth\n"%(name,pos))
		deleteCounter = deleteCounter + 1


#Deleting positions with anomalous coverage
print("Removed a further ", deleteCounter," anomalous depth positions")

#Loop again to avoid RuntimeError: dictionary changed size during iteration
#Writing singles
for name in poly:
	for pos in poly[name]:
		tot=0.0
		train=[]
		nucleotides=[]
		tmp = variants[name][pos][1:]
		for N in poly[name][pos]:
			nb=tmp.count(N)
			train.append(nb)
			nucleotides.append(N)
			tot+=nb
		i=0
		while i < len(train):
			nb=train[i]
			N = nucleotides[i]
			#Hazem: make sure single alleles are removed from singles but tracked in singleton output for later merge
			if N in alpha and len(train) > 1:
				h.write(name + "\t" + str(pos) + "\t" + N + "\t" +  str(nb/tot) + "\t" + str(tot) + "\t" + str(nb) + "\n")
			i+=1
	variants[name]=""

print("Wrote %s variants and filtered %s variants\n of which %s are above filtering thresholds and %s are below filtering threshold"%(passingVariants,filtered, filteredAbove, filteredBelow))
h.close()



def gap(cigar):
	final=[]
	resu=""
	tag=""
	i=0
	while i < len(cigar):
		try:
			digit = int(cigar[i])
			resu += cigar[i]
		except ValueError:
			tag=cigar[i]
			final.append([int(resu),tag])
			resu=""
		i+=1
	return final
			




assess=[]

missed_genes=[]


pair={}
for name1 in poly:
	pair[name1]={}
	pair[name1][name1]={}
	for pos1 in poly[name1]:
		pair[name1][name1][pos1]={}
		for pos2 in poly[name1]:
			if pos1 != pos2:
				pair[name1][name1][pos1][pos2]={}




for read in dico:
	if len(dico[read]) == 2:
		part1,part2=dico[read][0],dico[read][1]
		gene1,gene2=part1[1],part2[1]
		if gene1 != gene2:
			name1,name2=gene[gene1],gene[gene2]
			if  name1 in poly and name2 in poly and name2 not in pair[name1]:
				pair[name1][name2]={}
				for pos1 in poly[name1]:
					pair[name1][name2][pos1]={}
					for pos2 in poly[name2]:
						pair[name1][name2][pos1][pos2]={}

print("GO")

here=0
missed=0

counter=0
both_good=0
same_good=0
double=0
single=0
same=0
for read in dico:
	bad,good=0,0
	tag=0
	if len(dico[read]) == 1:
		single+=1
	if len(dico[read]) == 2:
		double+=1
		#print("new read",dico[read][0])
		part1,part2=dico[read][0],dico[read][1]
		cigar1,cigar2=part1[-2],part2[-2]
		gene1,gene2=part1[1],part2[1]
		read1,read2=part1[2],part2[2]
		phred1,phred2 = part1[-1],part2[-1]
		messedup1,messedup2=0,0
		debut1,fin1,debut2,fin2=0,0,0,0
		if "I" in cigar1 or "D" in cigar1:
			code1 = gap(cigar1)
			new1=""
			new_phred1=""
			I=0
			for thing in code1:
				
				i=0
				while i < thing[0]:
					#print(I,i,thing[1])
					if thing[1]=="M":
						try:
							new1+= read1[I]
							new_phred1 += phred1[I]
						except IndexError:
							pass
					elif thing[1]=="D" and i==0:
						try:
							new1+= read1[I]
							new_phred1 += phred1[I]
						except IndexError:
							pass
						new1+="-"
						new_phred1 += "-"
					elif thing[1] == "D" and i >0:
						new1+="-"
						new_phred1 += "-"
					i+=1
					I+=1
			if gene[gene1] in ref:
				seq1=ref[gene[gene1]]
				pos1=int(part1[3])-1
				if pos1 > shift[gene1] and len(code1) <= 3:
					#print(cigar1,pos1,len(new1))
					#print(new1)
					#print(seq1[pos1-shift[gene1]:])
					#exit()
					read1=new1
					phred1=new_phred1
					tag=1
			messedup1=1
		else:
			debut1=int(cigar1.strip("M"))
			fin1=len(read1)
		if "I" in cigar2 or "D" in cigar2:
			code2 = gap(cigar2)
			#print("\n")
			#print(cigar2)
			new2=""
			new_phred2=""
			I=0
			for thing in code2:
				
				i=0
				while i < thing[0]:
					#print(I,i,thing[1])
					if thing[1]=="M":
						try:
							new2+= read2[I]
							new_phred2 += phred2[I]
						except IndexError:
							pass
					elif thing[1]=="D" and i==0:
						try:
							new2+= read2[I]
							new_phred2 += phred2[I]
						except IndexError:
							#print(I,len(read2),read2)
							pass
							#exit()
						new2+="-"
						new_phred2 += "-"
					elif thing[1] == "D" and i >0:
						new2+="-"
						new_phred2 += "-"
					i+=1
					I+=1
			if gene[gene2] in ref:
				seq2=ref[gene[gene2]]
				pos2=int(part2[3])-1
				if pos2 > shift[gene2] and len(code2) <= 3:
					#print(cigar2,pos2,len(new2))
					#print(new2)
					#print(seq2[pos2-shift[gene2]:])
					#exit()
					read2=new2
					phred2=new_phred2
					tag=1
			messedup2=1
		else:
			debut2=int(cigar2.strip("M"))
			fin2=len(read2)
			
		
		
		
		pos1=int(part1[3])-1
		start1=pos1-shift[gene1]
		if start1 < 0:
			start1=0
			adj1=shift[gene1]-pos1
		else:
			adj1=0
		start2=pos2-shift[gene2]
		if start2 < 0:
			start2=0
			adj2=shift[gene2]-pos2
		else:
			adj2=0
		diff=abs(pos1-pos2+1)
		if gene[gene1] in ref:
			seq1=ref[gene[gene1]]
			
			sub_read1=read1[adj1:]
			sub_seq1=seq1[pos1-shift[gene1] + adj1:]                  				        # Rescaling position
			
			L=min([len(sub_read1),len(sub_seq1)])
			if L >= 10:	
				#print(cigar1,pos1,debut1,sub_read1[:debut1],sub_seq1)
				if 1==1:
				#	print("numbers:",pos1,shift[gene1], adj1)
					#length1,len(read1),sub_read1[len(read1)-fin1:],sub_seq1[len(read1)-fin1+1:])
					#print(cigar1,"pos1=",pos1,"adj1=",adj1)
					#print(sub_read1[:10])
					#print(sub_seq1[:10])
					#print(">read")
					#print(sub_read1)
					#print(">seq")
					#print(sub_seq1)
					score1=0.0
					i=0
					if messedup1 >=1:
						five_first1,five_first2=sub_read1[:5],sub_seq1[:5]
						#print(five_first1,five_first2)
						five_last1,five_last2=sub_read1[L-5:L],sub_seq1[L-5:L]
						#print(five_last1,five_last2)
						if five_first1 == five_first2 and five_last1== five_last2:
							while i < L:
								N1,N2=sub_read1[i],sub_seq1[i]
								if N1 in alpha and N2 in alpha:
									if N1 == N2:
										score1+=1
								i+=1
							if score1/L >= 0.90:
								good+=1
							else:
								bad+=1
							assess.append(score1/L)
						else:
							bad+=1
					else:
						if 1==1:
							while i < L:
								N1,N2=sub_read1[i],sub_seq1[i]
								if N1 in alpha and N2 in alpha:
									if N1 == N2:
										score1+=1
								i+=1
							if score1/L >= 0.80:
								#bad+=1
								good+=1
							elif score1/L >= 0.5:
								bad+=1
							else:
								bad+=1
								#print("shift",shift[gene1],"seq=",len(seq1))
								#print(">" + gene1)
								#print(sub_read1)
								#print(len(read1),score1/L)
							assess.append(score1/L)
						else:
							bad+=1
		else:
			if gene1 not in missed_genes:
				missed_genes.append(gene1)


		pos2=int(part2[3])-1
		start2=pos2-shift[gene2]
		if start2 < 0:
			start2=0
			adj2=shift[gene2]-pos2
		else:
			adj2=0
		start2=pos2-shift[gene2]
		if start2 < 0:
			start2=0
			adj2=shift[gene2]-pos2
		else:
			adj2=0
		diff=abs(pos2-pos2+1)
		if gene[gene2] in ref:
			seq2=ref[gene[gene2]]
			
			sub_read2=read2[adj2:]
			sub_seq2=seq2[pos2-shift[gene2] + adj2:]
			
			L=min([len(sub_read2),len(sub_seq2)])
			if L >= 10:	
				#print(cigar2,pos2,debut2,sub_read2[:debut2],sub_seq2)
				if 2==2:
				#	print("numbers:",pos2,shift[gene2], adj2)
					#length2,len(read2),sub_read2[len(read2)-fin2:],sub_seq2[len(read2)-fin2+2:])
					#print(cigar2,"pos2=",pos2,"adj2=",adj2)
					#print(sub_read2[:20])
					#print(sub_seq2[:20])
					#print(">read")
					#print(sub_read2)
					#print(">seq")
					#print(sub_seq2)
					score2=0.0
					i=0
					#Changes here copied from LM in slack
					if messedup2 >=1:
						five_first1,five_first2=sub_read2[:5],sub_seq2[:5]
						#print(five_first2,five_first2)
						five_last1,five_last2=sub_read2[L-5:L],sub_seq2[L-5:L]
						#print(five_last2,five_last2)
						if five_first1 == five_first2 and five_last1== five_last2:
							while i < L:
								N1,N2=sub_read2[i],sub_seq2[i]
								if N1 in alpha and N2 in alpha:
									if N1 == N2:
										score2+=1
								i+=1
							if score2/L >= 0.90:
								good+=1
							else:
								bad+=1
							assess.append(score2/L)
						else:
							bad+=1
					else:
						if 2==2:
							while i < L:
								N1,N2=sub_read2[i],sub_seq2[i]
								if N1 in alpha and N2 in alpha:
									if N1 == N2:
										score2+=1
								i+=1
							if score2/L >= 0.80:
								#bad+=2
								good+=1
							elif score2/L >= 0.5:
								bad+=1
							else:
								bad+=1
								#print("shift",shift[gene2],"seq=",len(seq2))
								#print(">" + gene2)
								#print(sub_read2)
								#print(len(read2),score2/L)
							assess.append(score2/L)
						else:
							bad+=1
		else:
			if gene2 not in missed_genes:
				missed_genes.append(gene2)
				
		
		
		if gene[gene1] in poly:																			# Compare Within READ1
			for pos in poly[gene[gene1]]:
				#print(rescale1, pos1)
				if pos >= pos1 -shift[gene1] and  pos < pos1 -shift[gene1] + len(read1):
					#print(pos1)
					move=shift[gene[gene1]] + pos - pos1 - 1
					#print("pos=",pos,"pos1=",pos1,"shift=",shift[gene[gene1]],"adj1=",adj1,"move=",move)
					#print(">read", len(read1[adj1:]) )
					#print(read1[adj1:])
					#print(read1)
					#print(">ref")
					#print(ref[gene[gene1]][pos1-shift[gene1] + adj1:])
					#print(ref[gene[gene1]])
					#print(read1[adj1:][pos-shift[gene1]-1])
					#print(ref[gene[gene1]][pos])
					#print(poly[gene[gene1]][pos])
					N1 = read1[shift[gene[gene1]] + pos - pos1 ]
					Q1 = ord(phred1[shift[gene[gene1]] + pos - pos1 ]) - 33
					if Q1 >= 20:
						N = ref[gene[gene1]][pos]
						if N1 not in poly[gene[gene1]][pos]:
							missed+=1
						else:
							here+=1
							#Shift the four loop below to be inside the else so that binome is constructed only if allele is in position (from singles)
							for POS in  poly[gene[gene1]]:
								if POS >= pos1 -shift[gene1] and  POS < pos1 -shift[gene1] + len(read1) and pos !=POS:
									N2 = read1[shift[gene[gene1]] + POS - pos1 ]
									#print("11",cigar1,len(read1), len(phred1) )
									Q2 = ord(phred1[shift[gene[gene1]] + POS - pos1 ]) - 33
									if Q2 >= 20:
										N2 = ref[gene[gene1]][POS]
										if N2 not in poly[gene[gene1]][POS]:
											missed+=1
										else:
											here+=1
											#make the shift in the second else of the N2 of binome if not in POS
											name1=gene[gene1]
											binome = N1 + N2
											#print(name1,name1,pos,POS,binome)
											#print(pair[name1][name1].keys())
											if binome not in pair[name1][name1][pos][POS]:
												pair[name1][name1][pos][POS][binome] =1
											else:
												pair[name1][name1][pos][POS][binome] +=1

											#if "MZNG01000095.1_7031_7990" == name1 and pos == 446 and POS == 545 and "MZNG01000095.1_7031_7990" in pair["MZNG01000095.1_7031_7990"] and 446 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"] and 545 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][446]:
											#	print("FIRST",pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][446][545])
											#if "MZNG01000095.1_7031_7990" == name1 and pos == 545 and POS == 446 and "MZNG01000095.1_7031_7990" in pair["MZNG01000095.1_7031_7990"] and 545 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"] and 446 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][545]:
											#	print("FIRST_rev",pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][545][446])




		if gene[gene2] in poly:																			# Compare WIthin READ2
			for pos in poly[gene[gene2]]:
				#print(rescale1, pos1)
				if pos >= pos2 -shift[gene2] and  pos < pos2 -shift[gene2] + len(read2):
					#print(pos1)
					move=shift[gene[gene2]] + pos - pos2 - 1
					#print("pos=",pos,"pos1=",pos1,"shift=",shift[gene[gene1]],"adj1=",adj1,"move=",move)
					#print(">read", len(read1[adj1:]) )
					#print(read1[adj1:])
					#print(read1)
					#print(">ref")
					#print(ref[gene[gene1]][pos1-shift[gene1] + adj1:])
					#print(ref[gene[gene1]])
					#print(read1[adj1:][pos-shift[gene1]-1])
					#print(ref[gene[gene1]][pos])
					#print(poly[gene[gene1]][pos])
					N2 = read2[shift[gene[gene2]] + pos - pos2 ]
					Q2 = ord(phred2[shift[gene[gene2]] + pos - pos2 ]) - 33
					if Q2 >= 20:
						N = ref[gene[gene2]][pos]
						if N2 not in poly[gene[gene2]][pos]:
							missed+=1
						else:
							here+=1
							#Shift the four loop below to be inside the else so that binome is constructed only if allele is in position (from singles)
							for POS in  poly[gene[gene2]]:
								if POS >= pos2 -shift[gene2] and  POS < pos2 -shift[gene2] + len(read2) and pos != POS:
									N1 = read2[shift[gene[gene2]] + POS - pos2 ]
									Q1 = ord(phred2[shift[gene[gene2]] + POS - pos2 ]) - 33
									if Q1 >= 20:
										N = ref[gene[gene2]][POS]
										if N1 not in poly[gene[gene2]][POS]:
											missed+=1
										else:
											here+=1
											#make the shift in the second else of the N2/N1 of binome if not in POS	
											name2=gene[gene2]
											binome = N2 + N1
											#print(name1,name1,pos,POS,binome)
											#print(pair[name1][name1].keys())
											#if abs(pos - POS) > 160:
											#	print(pos,POS)
											#print(pos,POS)
											if binome not in pair[name2][name2][pos][POS]:
												pair[name2][name2][pos][POS][binome] =1
											else:
												pair[name2][name2][pos][POS][binome] +=1

											#if "MZNG01000095.1_7031_7990" == name2 and pos == 446 and POS == 545 and "MZNG01000095.1_7031_7990" in pair["MZNG01000095.1_7031_7990"] and 446 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"] and 545 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][446]:
											#	print("SECOND",pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][446][545])
											#if "MZNG01000095.1_7031_7990" == name2 and pos == 545 and POS == 446 and "MZNG01000095.1_7031_7990" in pair["MZNG01000095.1_7031_7990"] and 545 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"] and 446 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][545]:
											#	print("SECOND_rev",pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][545][446])


		if gene[gene1] in poly and gene[gene2] in poly:	
			name1=gene[gene1]																		# Compare BETWEEN READS
			for pos in poly[gene[gene1]]:
				#print(rescale1, pos1)
				if pos >= pos1 -shift[gene1] and  pos < pos1 -shift[gene1] + len(read1):
					#print(pos1)
					move=shift[gene[gene1]] + pos - pos1 - 1
					#print("pos=",pos,"pos1=",pos1,"shift=",shift[gene[gene1]],"adj1=",adj1,"move=",move)
					#print(">read", len(read1[adj1:]) )
					#print(read1[adj1:])
					#print(read1)
					#print(">ref")
					#print(ref[gene[gene1]][pos1-shift[gene1] + adj1:])
					#print(ref[gene[gene1]])
					#print(read1[adj1:][pos-shift[gene1]-1])
					#print(ref[gene[gene1]][pos])
					#print(poly[gene[gene1]][pos])
					N1 = read1[shift[gene[gene1]] + pos - pos1 ]
					Q1 = ord(phred1[shift[gene[gene1]] + pos - pos1 ]) - 33
					if Q1 >= 20:
						N = ref[gene[gene1]][pos]
						if N1 not in poly[gene[gene1]][pos]:
							missed+=1
						else:
							here+=1
							#Shift the four loop below to be inside the else so that binome is constructed only if allele is in position (from singles)
							for POS in  poly[gene[gene2]]:
								if POS >= pos2 -shift[gene2] and  POS < pos2 -shift[gene2] + len(read2) and pos != POS:
									N2 = read2[shift[gene[gene2]] + POS - pos2 ]
									Q2 = ord(phred2[shift[gene[gene2]] + POS - pos2 ]) - 33
									if Q2 >= 20:
										N = ref[gene[gene2]][POS]
										if N2 not in poly[gene[gene2]][POS]:
											missed+=1
										else:
											here+=1	
											#make the shift in the second else of the N2 of binome if not in POS		
											name2=gene[gene2]
											binome = N1 + N2
											#print(name1,name1,pos,POS,binome)
											#print(pair[name1][name1].keys())
											#if abs(pos - POS) > 160:
											#print(pos,POS)
											#print(pair[name1][name2].keys())
											#print(pair[name2][name1].keys())
											if binome not in pair[name1][name2][pos][POS]:
												pair[name1][name2][pos][POS][binome] =1												
											else:
												pair[name1][name2][pos][POS][binome] +=1


for pos1 in pair[name1][name1]:
	for pos2 in pair[name1][name1][pos1]:
		#print(pos1,pos2,pair[name1][name1][pos1][pos2])
		pass
						


seen={}
rev={}
"writing output"
h=open("tmp_pairs.txt","w")
for name1 in pair:
	for name2 in pair[name1]:
		for pos1 in pair[name1][name2]:
			for pos2 in pair[name1][name2][pos1]:
				#Skip if there is only one binome (singleton pair)
				tot=0.0
				nogap_binomes = 0
				#avoid unnecessary comparison into loops and saves up ~2 minutes.
				if len(pair[name1][name2][pos1][pos2]) == 0:
					continue
				for binome in pair[name1][name2][pos1][pos2]:
					#Changes copied from LM in slack to exclude writing gaps
					if pair[name1][name2][pos1][pos2][binome] > 1 and "-" not in binome:
						tot+=pair[name1][name2][pos1][pos2][binome]
						nogap_binomes = nogap_binomes + 1
				#do not waste time and skip. Erroneous position pairs. May be check for emptiness first.
				if (tot == 0.0):
					continue
				#skip position pairs with just one binome
				if (nogap_binomes == 1):
					continue
				#skip pairs of positions covered by less than 20 reads
				if (tot < 20.0):
					continue
				for binome in pair[name1][name2][pos1][pos2]:
					nb = pair[name1][name2][pos1][pos2][binome]
					#Changes copied from LM in slack to exclude writing gaps
					if nb >1 and "-" not in binome:
						resu1 = [name1 , name2 , str(pos1) , str(pos2) , binome[0] + binome[1]]
						resu2 = [name2 , name1 , str(pos2) , str(pos1) , binome[1] + binome[0]]
						resu1 = "\t".join(resu1)
						resu2 = "\t".join(resu2)
						short1,short2=name1 + name2 + str(pos1) + str(pos2), name2 + name1 + str(pos2) + str(pos1) 
						if short2 not in rev:
							rev[short1]="y"
							h.write(resu1 + "\t" + str(nb/tot) + "\t" + str(tot) + "\t" + str(nb) + "\n")
							seen[resu1]="y"
							seen[resu2]="y"
					

h.close()



#if "MZNG01000095.1_7031_7990" in pair and "MZNG01000095.1_7031_7990" in pair["MZNG01000095.1_7031_7990"] and 446 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"] and 545 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][446]:
#	print(pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][446][545])
#if "MZNG01000095.1_7031_7990" in pair and "MZNG01000095.1_7031_7990" in pair["MZNG01000095.1_7031_7990"] and 545 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"] and 446 in pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][545]:
#	print(pair["MZNG01000095.1_7031_7990"]["MZNG01000095.1_7031_7990"][545][446])


print(missed,"/",here)

print("single=",single)
print("double=",double)
print("good ones=",good)
print("bad ones=",bad)

print("missing genes",len(missed_genes))
print("ASSESS:",sum(assess)/len(assess),len(assess))
logfile.close()




