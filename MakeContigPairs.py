
# To run after linkage_groups.py
# transform the files of single and paired alleles from a gene coordinates to contig coordinates (i.e., group of genes linked)

parent={}
group={}
rescale={}
f=open("linkage_groups.txt","r")
for l in f:
	a=l.strip("\n").split("\t")
	group[a[0]]=a[1:]
	position=0
#Testing with group651 from Serratia Reference1 M1_S1
#	if a[0] == "group651":
#		print("Group651",end="\t")
	for name in a[1:]:
		parent[name]=a[0]
		rescale[name]=position
		sub=name.split("_")
		deb,fin=int(sub[-2]),int(sub[-1])
		#Dec23rd, 2024: position needs to be incremented with every iteration
		position=position + (fin -deb)
		#if a[0] == "group651":
		#	print(name,end="\t")
		#	print(position,end="\t")
	#if a[0] == "group651":
	#	print("")


f.close()



h=open("contig_pairs.txt","w")
f=open("tmp_pairs.txt","r")
for l in f:
	a=l.strip("\n").split("\t")
	id1,id2=a[0],a[1]
	contig = parent[id1]
	if contig!=parent[id1]:
		print("ABORTING, these genes should be part of the same group")
		exit()
	pos1,pos2=int(a[2]),int(a[3])
	new1,new2=pos1+rescale[id1],pos2+rescale[id2]
	binome=a[4]
	N1,N2=binome[0],binome[1]
	if new1 > new2:
		resu = [contig, str(new2) ,str(new1) , N2 + N1 ,a[5],a[6]]
	else:
		resu = [contig, str(new1) ,str(new2) , N1 + N2 ,a[5],a[6]]
# Testing Serratia Reference1 M1_S1
#	if id1 == "NZ_PQGI01000007.1_97689_99365" or id2 == "NZ_PQGI01000007.1_97689_99365":
#		print(resu)
	h.write("\t".join(resu) + "\n")

f.close()
h.close()



seen={}
I=0
h=open("contig_singles.txt","w")
f=open("tmp_singles.txt","r")
for l in f:
	a=l.strip("\n").split("\t")
	name=a[0]
	pos=int(a[1])#Sep 26 2024: move pos to avoid scope error
	
	#Issue when processing Shiva's lines
	#Sep 26 2024: Change try and catch to if else. NoPair are variant positions that are not paired (do not have a pair in tmp_pairs and hence no contig -> gene mapping in linkage)
	if name in parent:
		contig = parent[name]
		new = pos + rescale[name]
		h.write(contig + "\t" + str(new) + "\t" + "\t".join(a[2:]) + "\n" )
	else:
		if name not in seen:
			seen[name]="y"
			I+=1
		#Dec 3 2024: append gene to nopair to allow mapping to original gene when genotyping
		h.write("nopair" + str(I) + "&" + name + "\t" + str(pos) + "\t" + "\t".join(a[2:]) + "\n" )#Use pos with nopair to match tmp_singles

f.close()
h.close()




