#Version3


import numpy as np
import emcee
import random
import sys
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-s","--singles", help="Input with singles")
parser.add_argument("-p","--pairs", required=True, help="import with pairs")
parser.add_argument("-S","--samplename", required=True, help="Samplename")
parser.add_argument("--seed", default=43,type=int, help="Samplename")
args = parser.parse_args()

if not os.path.exists(args.singles):
    sys.exit("Input variant singletons does not exist. %s" %(args.singles))
#Check if sam file exists
if not os.path.exists(args.pairs):
    sys.exit("Input variants pairs file does not exist.")






def closest(unique1,unique2,freq1,freq2,landscape):
	memo1=[2,""]
	memo2=[2,""]
	second1,second2=[2,""],[2,""]
	i=0
	while i < 6:
		nb=landscape[i]
		dist1 = abs(freq1 - nb)
		dist2 = abs(freq2 - nb)
		if dist1 < memo1[0]:
			memo1 = [dist1,i+1]
		if dist2 < memo2[0]:
			memo2 = [dist2,i+1]
		i+=1
	i=0
	while i < 6:
		nb=landscape[i]
		dist1 = abs(freq1 - nb)
		dist2 = abs(freq2 - nb)
		if dist1 < second1[0] and dist1 >= memo1[0]:
			second1 = [dist1,i+1]
		if dist2 < second2[0] and dist1 >= memo2[0]:
			second2 = [dist2,i+1]
		i+=1
	return [memo1,memo2,second1,second2]


def closest_single(unique1,freq1,landscape):
	memo1=[2,""]
	second1=[2,""]
	third1=[2,""]
	i=0
	while i < 6:
		nb=landscape[i]
		dist1 = abs(freq1 - nb)
		if dist1 < memo1[0]:
			memo1 = [dist1,i+1]
		i+=1
	i=0
	while i < 6:
		nb=landscape[i]
		dist1 = abs(freq1 - nb)
		if dist1 < second1[0] and dist1 > memo1[0]:
			second1 = [dist1,i+1]
		i+=1
	i=0
	while i < 6:
		nb=landscape[i]
		dist1 = abs(freq1 - nb)
		if dist1 < third1[0] and dist1 > memo1[0] and dist1 > second1[0]:
			third1 = [dist1,i+1]
		i+=1
	return [memo1,second1,third1]



def does_complement(liste):
	liste=[str(liste[0]),str(liste[1])]
	liste.sort()
	peak1,peak2 = str(liste[0]),str(liste[1])
	if peak1 == peak2:
		answer=0
	elif peak1 in ["1","2","3"] and peak2 in ["1","2","3"]:
		answer=1
	elif liste == ["1","6"]:
		answer=1
	elif liste == ["2","5"]:
		answer=1
	elif liste == ["3","4"]:
		answer=1
	else:
		answer=0
	return answer


complement={}
seen={}
allele={}
#f=open("Sample" + str(sample_number) + "_150bp.txt","r")
#f=open("Sample" + str(sample_number) + "_core_Gapicola_MZNG01_150_DP10-No1AD_noFlank_frequency_singles_MClust.txt","r")
f=open(args.singles,"r")
#l=f.readline() no longer needed, I skipped the header since I no longer use R
for l in f:
	a=l.strip("\n").split("\t")
	name = a[0] + "&" + a[1] + "&" + a[2]
	resu = a[0] + "&" + a[1] 
	if resu not in seen:
		seen[resu]= name
	else:
		complement[name]=seen[resu]
		complement[seen[resu]] = name
		
	freq = float(a[3])
	if freq >= 0.01 or 1 - freq >= 0.01:
		allele[name]=freq

f.close()

print("allele:",name,allele[name])



locus={}
binomes={}
new_binomes={}
if 1==1:
	#f=open("Sample" + str(sample_number) + "_core_Gapicola_MZNG01_150_DP10-No1AD_noFlank_frequency_pairs_MClust_singles.txt","r")
	#Sample9_core_Gapicola_MZNG01_150_2g_3st_col2_DP10-No1AD_noFlank_frequency_singles.txt
	f=open(args.pairs,"r")
	for l in f:
		a=l.strip("\n").split("\t")
		name=a[0]
		pos1,pos2=a[1],a[2]
		pair= a[3]
		resu = name + "&" + pos1 + "&" + pos2 + "&" + pair 
		freq=float(a[4])
		if freq > 0.02 and 1 - freq > 0.02:
			binomes[resu]=freq
			resu2 = name + "&" + pos1 + "&" + pos2
			if resu2 not in new_binomes:
				new_binomes[resu2 ] = [pair]
			else:
				new_binomes[resu2 ].append(pair)
			resu3=name + "&" + pos1 + "&" + pos2 
			if resu3 not in locus:
				locus[resu3]=[]
			locus[resu3].append([pair,freq])
	f.close()


#print("binomes:",resu,binomes[resu]) #Check later why resu gives a Key Exception 
print("new_binomes:",resu2, new_binomes[resu2])

impossible={}
array={}
for resu in binomes:
	a=resu.split("&")
	group,pos1,pos2,pair=a[0],a[1],a[2],a[3]
	N1,N2=pair[0],pair[1]
	if group not in array:
		array[group]={}
		impossible[group]={}
	id1,id2=group + "&" + pos1 + "&" + N1,group + "&" + pos2 + "&" + N2
	if group + "&" + pos1 == group + "&" + pos2:
		if id1 not in impossible[group]:
			impossible[group][id1]=[id2]
		else:
			impossible[group][id1].append(id2)
		if id2 not in impossible[group]:
			impossible[group][id2]=[id1]
		else:
			impossible[group][id2].append(id1)
	if id1 not in array[group]:
		array[group][id1] = [id2]
	else:
		array[group][id1].append(id2)
	if id2 not in array[group]:
		array[group][id2] = [id1]
	else:
		array[group][id2].append(id1)


for group in array:
	if 1==2:
		genotypeA,genotypeB,genotypeC=[],[],[]
		all_ids = list(array[group].keys())
		all_ids.sort()
		i=0
		while i < len(all_ids)-1:
			id1=all_ids[i]
			id2=all_ids[i+1]
			if id1 != id2:
				list1 = list(array[group][id1])
				list2 = list(array[group][id2])
				yes=0
				no=0
				for variant1 in list1:
					#print(variant1)
					if variant1 in list2:
						yes=1
					if complement[variant1] in list2:
						no=1
				print(yes,no,variant1,complement[variant1] ,list1, list2 )
				
			i+=1





def log_prob(x, mu, cov):
	diff = x - mu
	return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))



#def log_prob(x, ivar):
    #return -0.5 * np.sum(ivar * x ** 2)


ndim = 1 # Number of dimensions

np.random.seed(args.seed)
#seed is also needed here as well
random.seed(args.seed)
means = np.random.rand(ndim)

cov = 0.5 - np.random.rand(ndim**2).reshape((ndim, ndim))
cov = np.triu(cov)
cov += cov.T - np.diag(cov.diagonal())
cov = np.dot(cov, cov)


nwalkers = 32
p0 = np.random.rand(nwalkers, ndim)

print("p0=",p0.shape)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=[means, cov])

print("log prob of (po[0], means, cov) in line below")
print(log_prob(p0[0], means, cov))

print("p0 = np.random.rand(nwalkers, ndim) = %s" % (p0))

state = sampler.run_mcmc(p0, 500)
samples= sampler.get_chain()
small,big=min(samples[:,0,0]),max(samples[:,0,0])
L = big - small
sampler.reset()

p0 = np.random.rand(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=[means, cov])
	

sampler.run_mcmc(state, 1)


samples = sampler.get_chain()    # MAKE FLAT: samples = sampler.get_chain(flat=True)

print(samples.shape)


jumps=[5,10,100,"new"]
micro=[0.33,0.33,0.34]

jump=1


toto={}
for I in range(30):
	toto[I]=[]

ALL=[]
i=0
while i <=3000:
	sampler.run_mcmc(state, 1)
	samples = sampler.get_chain()
	state = samples[i,:]
	move1 = samples[-1,0,0]
	for I in range(30):
		move1 = samples[-1,I,0]
		toto[I].append(move1)
		ALL.append(move1)
	i+=1



together={}
together["1"]=["1","2","3","4","5"]
together["2"]=["1","2","3","4","6"]
together["3"]=["1","2","3","5","6"]
together["4"]=["1","2"]
together["5"]=["1","3"]
together["6"]=["2","3"]


paired_peaks={}
paired_peaks["1"]=["1","4","5"]
paired_peaks["2"]=["2","4","6"]
paired_peaks["3"]=["3","5","6"]
paired_peaks["4"]=["1","2","4"]
paired_peaks["5"]=["1","3","5"]
paired_peaks["6"]=["2","3","6"]

fix = np.mean(ALL)

sampler.reset()

individual_seen={}





sampler.run_mcmc(state, 1)
action="move"
series=0
SCORE=100000000000000000000000000
i=1
I=0
converge_steps = 0
old_tmp_micro = []
#temporary adding buffering=1 since I let job timeout on LongLeaf
h=open("moves_" + args.samplename + ".txt","w",buffering=1)
while i <=50000:
	h.write(str(i) + "\t" + str(micro[0]) + "\t" + str(micro[1]) + "\t" + str(micro[2]) + "\n")
	mono_attribution={}
	haplotype={}
	bad_peak={}
	miss_pairing=0
	if i%10==0:
		big_jump="yes"
	else:
		big_jump="no"
	score=0
	TAG=0
	sampler.run_mcmc(state, 1)
	samples = sampler.get_chain()
	#print(samples)
	#print("state",state)
	#print("1",samples[-1,0,0],samples[-1,1,0])
	if action == "move":
		jump=1
		
	if jump== "new":
		tag=1
		while tag==1:
			a = random.uniform(0.0, 1.0)
			b = random.uniform(0.0, 1.0)
			if a + b < 1:
				tag=0
		c = 1-a-b
		tmp = [a,b,c]
		random.shuffle(tmp)
			#print("NEW")
		a,b,c=tmp[0],tmp[1],tmp[2]
	else:
		a,b,c=micro[0],micro[1],micro[2]
		sign1=random.uniform(0.0, 1.0)
		move1 = 0.05 * (samples[-1,0,0] - fix) * jump
		if sign1 < 0.5:
			move1 = float(-move1)
		sign2=random.uniform(0.0, 1.0)
		move2 = 0.05 * (samples[-1,1,0] - fix) * jump
		if sign2 < 0.5:
			move2 = float(-move2)
	if abs(move1) < 1 and abs(move2) < 1:
		if series >= 20 or big_jump== "yes":
			jump=random.choice(jumps)
			
		else:
			jump=1
		if TAG > 10:
			jump="new"
		#print(move1,move2)
		if jump== "new":
			tag=1
			while tag==1:
				a = random.uniform(0.0, 1.0)
				b = random.uniform(0.0, 1.0)
				if a + b < 1:
					tag=0
			c = 1-a-b
			tmp = [a,b,c]
			random.shuffle(tmp)
			#print("NEW")
			a,b,c=tmp[0],tmp[1],tmp[2]
		a,b = a+move1 , b+move2
		if a < -0.001 or a > 1.001:
			a= a - move1
		if b < -0.001 or b > 1.001:
			b= b - move2
		if a + b <= 1.001 and a + b >= -0.001:
			c = 1-a-b
			memo_micro=list(micro)																		# CHANGED
			micro = [a , b , c]
			i+=1
		else:
			TAG+=1
	else:
		TAG+=1
	if TAG==0:
		landscape=list(micro)
		landscape.sort()
		tmp_micro= list(micro)
		tmp_micro.sort()
		one,two,three=tmp_micro[0],tmp_micro[1],tmp_micro[2]
		landscape.append(tmp_micro[0]+tmp_micro[1])
		landscape.append(tmp_micro[0]+tmp_micro[2])
		landscape.append(tmp_micro[1]+tmp_micro[2])
		assigned_peak={}
		assign={}
		for resu3 in locus:
			variantsA,variantsB=[],[]
			distancesA,distancesB=[],[]
			codeA,codeB=[],[]
			sub=resu3.split("&")
			group=sub[0]
			pos1,pos2=sub[1],sub[2]
			tmp=[]
			for doublon in locus[resu3]:
				N1,N2=doublon[0][0],doublon[0][1]
				freq = doublon[1]
				
				unique1,unique2=group + "&"  + pos1 + "&" + N1 ,  group + "&"  + pos2 + "&" + N2
				resu3=group + "&"  + pos1 + "&" + pos2
				out = closest(unique1,unique2,allele[unique1],allele[unique2],landscape)
				individual_seen[unique1],individual_seen[unique2]="y","y"
				if group not in haplotype:
					haplotype[group]={}
					bad_peak[group]={}
				peak1,peak2=str(out[0][1]),str(out[1][1])
				#if group + "&"  + pos1 == "group1&126" and N1 == "G":
					#print(unique1,unique2,freq,out[0][1],out[1][1],[allele[unique1],allele[unique2]]) 
				if unique1 not in haplotype[group]:
					haplotype[group][unique1]={}
					bad_peak[group][unique1]={}
				if unique2 not in haplotype[group]:
					haplotype[group][unique2]={}
					bad_peak[group][unique2]={}
				if peak1 not in haplotype[group][unique1]:
					haplotype[group][unique1][peak1]=0.0
				if peak2 not in haplotype[group][unique2]:
					haplotype[group][unique2][peak2]=0.0
				if peak1 not in bad_peak[group][unique1]:
					bad_peak[group][unique1][peak1]=0.0
				if peak2 not in bad_peak[group][unique2]:
					bad_peak[group][unique2][peak2]=0.0
				if peak1 in paired_peaks[peak2]:
					haplotype[group][unique1][peak1]+=1
					haplotype[group][unique2][peak2]+=1
				else:
					bad_peak[group][unique1][peak1]+=1
					bad_peak[group][unique2][peak2]+=1
				mono_attribution[unique1]=peak1
				mono_attribution[unique2]=peak2
				if peak1 != peak2:
					if peak1 not in paired_peaks[peak2]:
						miss_pairing += 1
						score += 1
				tmp.append(out[0][0])
				tmp.append(out[1][0])
				if N1 not in variantsA:
					variantsA.append(N1)
					codeA.append(out[0][1])
					distancesA.append(out[0][0])
				if N2 not in variantsB:
					variantsB.append(N2)
					codeB.append(out[1][1])
					distancesB.append(out[1][0])
			if len(variantsA) == 2 and len(variantsB)==2:
				#allA1,allA2,allB1,allB2=group + "&" + pos1 + "&" + variantsA[0], group + "&" + pos1 + "&" + variantsA[1] , group + "&" + pos2 + "&" + variantsB[0] , group + "&" + pos2 + "&" + variantsB[1]
				answer1,answer2=does_complement(codeA),does_complement(codeB)
				#print(variantsA,variantsB,codeA,codeB,answer1,answer2,landscape)
				if answer1==0 or answer2==0:
					#print("Not together")
					score += len(tmp)*2
				else:
					score += sum(tmp)
	if action == "move":
		out=open("genotypes_" + args.samplename + ".txt","w")
		for group in haplotype:
			for variant in haplotype[group]:
				#tmp=list(set(haplotype[group][variant]))
				if 1==1:
					memo,MAX,tot="",0,0.0
					for peak in haplotype[group][variant]:
						tot+=haplotype[group][variant][peak]
						if haplotype[group][variant][peak] > MAX:
							MAX =  haplotype[group][variant][peak]
							memo=peak
					wrong=0
					for peak in bad_peak[group][variant]:
						tot+=1
					if tot == 0:
						percent="NA"
					else:
						percent=100*MAX/tot
					out.write(group + "\t" + variant + "\t"  + str(allele[variant]) + "\t" + memo + "\t" + str(percent)  + "\t" + str(int(tot)) + "\n")
		for single in allele:
			sub = single.split("&")
			group = sub[0]
			if group not in haplotype:
				what = closest(single,single,allele[single],allele[single],landscape)
				peak = what[0][1]
				out.write(group + "\t" + single + "\t" + str(allele[single]) + "\t" + str(peak) + "\t" + "100.00\t1\n")
			#elif single not in haplotype[group]:
				#print(single,"missing")
		out.close()
	if i==0:
		SCORE = 100000000000000000000000
	if score < SCORE and miss_pairing > 0:
		action="move"
		
		print(i,"score=",SCORE,action,micro,jump,"         moves=",move1,move2,miss_pairing)
		key=open("key_genotypes_" + args.samplename + ".txt","w")
		key.write("1\t" + str(tmp_micro[0]) + "\n2\t" + str(tmp_micro[1]) + "\n3\t" + str(tmp_micro[2]) + "\n" )
		key.close()
		#Dec 9 2024: Hazem Longleaf change
		#Terminate loop when solutions converged and movements are stable. Customise later
		#Doing temporary sort in order not to mess up storage
		move_diff = sum([abs(x - y) for x,y in zip( sorted(tmp_micro) ,sorted(old_tmp_micro) )])
		
		#if move_diff < 0.001:
		if move_diff == 0:
			converge_steps = converge_steps + 1
			print("Move diff %s" %(move_diff))
		else:
			converge_steps = 0
			print("Resetting movements count")
			print(old_tmp_micro)
			print(tmp_micro)
			print("Move diff %s" %(move_diff))
		#tracking change in optimum
		old_tmp_micro = tmp_micro
	else:
		action="don't move"
		print("don't move",jump)
		#Terminate loop when solutions converged and movements are stable. Customise later
		converge_steps = converge_steps + 1
		print("Steps %s miss_pairing %s" %(converge_steps, miss_pairing))
	if score < SCORE and score > 0 and miss_pairing > 0:
		I+=1
		state = samples[I,:]
		SCORE = float(score)
		series=0
	else:
		series +=1
		micro = list(memo_micro)																			#### CHANGED

	#if score < SCORE and score > 0 and miss_pairing > 0 and converge_steps > 50:
	if score > 0 and miss_pairing > 0 and converge_steps > 100:
		print("number of steps %s "%(i))
		sys.exit("Movements stable for 100 steps. Exit")
	











print("END")

exit()


