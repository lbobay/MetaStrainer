
# To run after the single and paired alleles files are created 
# Create the linkage groups from the polymorphism files


link={}
seen={}
dico={}
nb=0
f=open("tmp_pairs.txt","r")
for l in f:
	a=l.strip("\n").split("\t")
	id1,id2=a[0],a[1]
	resu=id1 + "&" + id2
	if resu not in seen:
		nb+=1
		dico[nb]={}
		dico[nb][id1]="y"
		dico[nb][id2]="y"
		seen[resu]="y"
		print("id1: %s\tid2: %s"%(id1,id2))
		print("dico nb",nb)




found=1
while found ==1:
	found=0
	delete={}
	for nb1 in dico:
		for nb2 in dico:
			if nb1!=nb2:
				for id1 in dico[nb1]:
					if id1 in dico[nb2]:
						found=1
						merge=[nb1,nb2]
						break
		if found==1:
			break
	if found==1:
		print("merging",merge)
		nb1=merge[0]
		nb2=merge[1]
		for id2 in dico[nb2]:
			dico[nb1][id2]="y"
			if "MZNG01000001.1_40622_41242" == id2:
				print(nb1)
		del dico[nb2]
	
	
h=open("linkage_groups.txt","w")
I=0
nb=0
while nb < max(dico.keys()):
	nb+=1
	if nb in dico:
		I+=1
		group = "group" + str(I)
		print(nb,group)
		tmp = list(dico[nb])
		tmp.sort()
		h.write(group + "\t" + "\t".join(tmp) + "\n")
	
exit()



