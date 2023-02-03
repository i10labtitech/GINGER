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
import os

#塩基配列
fna=open(sys.argv[1],"r")
faa=open(sys.argv[2],"r")
gff=open(sys.argv[3],"r")

faa_dict={}
data1=""
for a in faa:
	a = a.strip()
	if len(a) > 0:
		if a[0] == ">":
			data1=a.strip(">")
			faa_dict[data1]=["",0]
		else:
			faa_dict[data1][0]=faa_dict[data1][0]+a
			faa_dict[data1][1]=faa_dict[data1][1]+len(a)
faa.close()

fna_dict={}
data2=""
for n in fna:
	n=n.strip()
	if len(n) > 0:
		if n[0] == ">":
			data2=n.strip(">")
			fna_dict[data2]=["",0]
		else:
			fna_dict[data2][0]=fna_dict[data2][0]+n
			fna_dict[data2][1]=fna_dict[data2][1]+len(n)
fna.close()

length_fna=len(fna_dict)
length_faa=len(faa_dict)
rm_dict={}
if length_fna != length_faa:
	print("The number of nucleotide sequences and amino acid sequences is different !!")
	sys.exit()
for i in faa_dict:
	nucl_len=fna_dict[i][1]
	amin_len=faa_dict[i][1]*3
	if nucl_len != amin_len:
		dif=nucl_len-amin_len
		if "*" in faa_dict[i][0]:
			rm_dict[i]=dif
		#print(i+"\t"+str(dif))
		#print("#nucl:"+fna_dict[i][0])
		#print("#amin:"+faa_dict[i][0])

n=0
tmp=open(sys.argv[3]+".tmp","w")
for g in gff:
	g=g.strip()
	g_type=g.split("\t")[2]
	if g_type == "mRNA" and n==0:
		n=1
		tmp.write(g)
	elif g_type == "mRNA" and n!=0:
		tmp.write("\n"+g)
	else:
		tmp.write("@@_@@"+g)
tmp.write("\n")
gff.close()
tmp.close()
tmp2=open(sys.argv[3]+".tmp","r")
for t in tmp2:
	t=t.strip()
	tt=t.split("@@_@@")
	word=tt[0].split("\t")[8].split("ID=")[1]
	if word in rm_dict:
		strand=tt[0].split("\t")[6]
		if strand=="+":
			mrna=tt[0].split("\t")
			num=int(mrna[4]) - rm_dict[word]
			mrna[4]=str(num)
			tt[0]="\t".join(mrna)
			cds=tt[-1].split("\t")
			num=int(cds[4]) - rm_dict[word]
			cds[4]=str(num)
			tt[-1]="\t".join(cds)
		elif strand=="-":
			mrna=tt[0].split("\t")
			num=int(mrna[3]) + rm_dict[word]
			mrna[3]=str(num)
			tt[0]="\t".join(mrna)
			cds=tt[-1].split("\t")
			num=int(cds[3]) + rm_dict[word]
			cds[3]=str(num)
			tt[-1]="\t".join(cds)
		sys.stdout.write("\n".join(tt)+"\n")
	else:
		sys.stdout.write("\n".join(tt)+"\n")
tmp2.close()
os.system("rm "+sys.argv[3]+".tmp")
