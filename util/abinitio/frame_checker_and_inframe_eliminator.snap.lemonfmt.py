#coding:utf-8

# Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology
# 
# This file is part of GINGER.
# 
# GINGER is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GINGER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with GINGER; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import sys
import random
import os
vr=sys.version_info.major
if vr == 2:
	import string
elif vr == 3:
	string=str

genome=open(sys.argv[1],"r")
faa=open(sys.argv[2],"r")
gff=open(sys.argv[3],"r")

g1_tag=str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))
g1=open("genome."+g1_tag+".g1","w")
n=0
for g in genome:
	g = g.strip()
	if len(g) > 0:
		if g[0] == ">":
			if n == 0:
				n+=1
				g1.write(g.strip(">")+"\t")
			else:
				n+=1
				g1.write("\n"+g.strip(">")+"\t")
		else:
			g1.write(g)
g1.write("\n")
genome.close()
g1.close()
##########################################
dict_genome={}
g2=open("genome."+g1_tag+".g1","r")
for e in g2:
	e = e.strip()
	ee = e.split("\t")
	dict_genome[ee[0]]=ee[1]
g2.close()
os.system("rm genome."+g1_tag+".g1")
#########################################
f1_tag=str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))
f1=open("gene."+f1_tag+".f1","w")
m=0
#box=[]
for f in faa:
	f = f.strip()
	if len(f) > 0:
		if f[0] == ">":
			if m == 0:
				m+=1
				f1.write(f.strip(">")+"\t")
			else:
				m+=1
				f1.write("\n"+f.strip(">")+"\t")
		else:
			f1.write(f)
f1.write("\n")
faa.close()
f1.close()
##########################################
f_f=open("gene."+f1_tag+".f1","r")
e_box=[]
for z in f_f:
	z = z.strip()
	zz = z.split("\t")
	ids=zz[0]
	seqe=zz[1]
	if "*" in seqe[:-1]:
		e_box.append(ids)
f_f.close()
#print(e_box)
##########################################
dict_gene={}
f2=open("gene."+f1_tag+".f1","r")
for a in f2:
	a = a.strip()
	aa = a.split("\t")
	#box.append(aa[0])
	dict_gene[aa[0]]=aa[1]
f2.close()
os.system("rm gene."+f1_tag+".f1")
#########################################
#for t in box:
#	print(">"+t)
#	print(dict_gene[t])
########################################
#########################################
tmp_tag=str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))
#print(tmp_tag)
#tmp_tag=sys.argv[1].split("/")[-1]+"."+sys.argv[2].split("/")[-1]+"."+sys.argv[3].split("/")[-1]
tmp=open("tmp."+tmp_tag+".tmp","w")
dict_gs={}
box=[]
for s in gff:
	s = s.strip()
	if len(s) > 0:
		ss = s.split("\t")
		sca = ss[0]
		typ = ss[2]
		st = int(ss[3]) - 1
		en = int(ss[4])
		di = ss[6]
		if typ == "mRNA":
			tmp.write(s+"\n")
		elif typ == "CDS":
			tmp.write(s+"\n")
			iD = ".".join(ss[8].split(";")[0].split("=cds")[1].split(".")[1:])
			#iD = ss[8].split(";")[0].split("=cds.model.")[1]
			if di=="+":
				if iD in box:
					dict_gs[iD]+=dict_genome[sca][st:en]
				elif iD not in box:
					box.append(iD)
					dict_gs[iD]=dict_genome[sca][st:en]
			elif di=="-":
				if iD in box:
					pre=dict_genome[sca][st:en]
					rev=pre[::-1]
					tran= rev.translate(string.maketrans("ATGCatgc","TACGtacg"))
					dict_gs[iD]+=tran
				elif iD not in box:
					box.append(iD)
					pre=dict_genome[sca][st:en]
					rev=pre[::-1]
					tran= rev.translate(string.maketrans("ATGCatgc","TACGtacg"))
					dict_gs[iD]=tran
tmp.close()
gff.close()
#########################################
#for t in box:
#       print(">"+t)
#       print(dict_gs[t])
########################################
DNA2Protein = {
        'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
        'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
        'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
        'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

        'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R',
        'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R',
        'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R',
        'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R',

        'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S',
        'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S',
        'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R',
        'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R',

        'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
        'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
        'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
        'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G',

        'CTN' : 'L', 'GTN' : 'V', 'TCN' : 'S', 'CCN' : 'P',
        'ACN' : 'T', 'GCN' : 'A', 'CGN' : 'R', 'GGN' : 'G',
}
#######################################
dict_tran={}
for l in box:
	dict_tran[l]={}
	seq0=dict_gs[l][0:]
	seq1=dict_gs[l][1:]
	seq2=dict_gs[l][2:]
	p1=""
	for i in range(0, int(len(seq0) / 3)):
		codon = seq0[3 * i : 3 * i + 3]
		if codon.upper() in DNA2Protein:
			p1 += DNA2Protein[codon.upper()]
		else:
			p1 += 'X'
	p2=""
	for x in range(0, int(len(seq1) / 3)):
		codon1 = seq1[3 * x : 3 * x + 3]
		if codon1.upper() in DNA2Protein:
			p2 += DNA2Protein[codon1.upper()]
		else:
			p2 += 'X'
	p3=""
	for y in range(0, int(len(seq2) / 3)):
		codon2 = seq2[3 * y : 3 * y + 3]
		if codon2.upper() in DNA2Protein:
			p3 += DNA2Protein[codon2.upper()]
		else:
			p3 += 'X'
	if p1 == dict_gene[l]:
		dict_tran[l]= 0
	elif p2 == dict_gene[l]:
		dict_tran[l]= 1
	elif p3 == dict_gene[l]:
		dict_tran[l]= 2
	else:
		dict_tran[l]= 3
		print("#>"+l+"_snap\n#"+dict_gene[l])
		print("#>"+l+"_trans_frame_1\n#"+p1)
		print("#>"+l+"_trans_frame_2\n#"+p2)
		print("#>"+l+"_trans_frame_3\n#"+p3)
		print("#>"+l+"_nucl\n#"+seq0)
	#print(l+"\t"+str(dict_tran[l]))
#######################################
gff_2nd=open("tmp."+tmp_tag+".tmp","r") 
#print(tmp_tag)
check_box=[]
trans_table={2:1,1:2,0:0}
amari=0
s_w = 0
for k in gff_2nd:
	k = k.strip()
	kk = k.split("\t")
	Tp=kk[2]
	St=int(kk[3])
	En=int(kk[4])
	if Tp == "mRNA":
		ID=kk[8].split(";")[0].split("=")[1]
		#ID=kk[8].split(";")[0].split("=model.")[1]
		if ID not in e_box:
			if s_w == 0:
				s_w += 1
				print(k)
			elif s_w > 0:
				print("\n"+k)
	elif Tp == "CDS":
		ID=".".join(kk[8].split(";")[0].split("=cds")[1].split(".")[1:])
		#ID=kk[8].split(";")[0].split("=cds.model.")[1]
		if ID not in e_box:
			if ID not in check_box:
				check_box.append(ID)
				kk[7]=str(dict_tran[ID])
				TSt=St+int(dict_tran[ID])
				size=En-TSt+1
				amari=size%3
				print("\t".join(kk))
			elif ID in check_box:
				slid=trans_table[amari]
				kk[7]=str(slid)
				TSt=St+slid
				size=En-TSt+1
				amari=size%3
				print("\t".join(kk))
gff_2nd.close()
os.system("rm tmp."+tmp_tag+".tmp")
####################################
