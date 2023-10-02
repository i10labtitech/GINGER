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

snap=open(sys.argv[1],"r")
tmp=open(sys.argv[1]+".tMp","w")
sw=0
mem=""

for i in snap:
	i = i.strip()
	ii = i.split("\t")
	ID=ii[-1]
	if mem!=ID and sw==0:
		tmp.write(i)
		sw+=1
		mem=ID
	elif mem!=ID and sw>0:
		tmp.write("\n"+i)
		sw+=1
		mem=ID
	elif mem==ID:
		tmp.write("@@@"+i)
snap.close()
tmp.close()

TMP=open(sys.argv[1]+".tMp","r")
pos_box=[]
info_box=[]
seq_id=""
gene_id=""
strand=""
number=0

for l in TMP:
	pos_box=[]
	info_box=[]
	seq_id=""
	gene_id=""
	strand=""
	number=0
	l=l.strip()
	ll=l.split("@@@")
	for n in ll:
		number+=1
		nn=n.split("\t")
		ids=nn[2]
		st=int(nn[3])
		en=int(nn[4])
		gene_id=nn[-1]
		seq_id=nn[0]
		strand=nn[6]
		if ids=="Einit":
			pos_box.append(st)
			pos_box.append(en)
			nn[2]="CDS"
			nn[-1]="ID=cds"+str(number)+"."+nn[-1]+";Parent="+nn[-1]
			info_box.append("\t".join(nn))
		elif ids=="Eterm":
			pos_box.append(st)
			pos_box.append(en)
			nn[2]="CDS"
			nn[-1]="ID=cds"+str(number)+"."+nn[-1]+";Parent="+nn[-1]
			info_box.append("\t".join(nn))
		elif ids=="Exon":
			pos_box.append(st)
			pos_box.append(en)
			nn[2]="CDS"
			nn[-1]="ID=cds"+str(number)+"."+nn[-1]+";Parent="+nn[-1]
			info_box.append("\t".join(nn))
		elif ids=="Esngl":
			pos_box.append(st)
			pos_box.append(en)
			nn[2]="CDS"
			nn[-1]="ID=cds"+str(number)+"."+nn[-1]+";Parent="+nn[-1]
			info_box.append("\t".join(nn))
	pb=sorted(pos_box)
	#print(pb)
	mrna_st=str(pb[0])
	mrna_en=str(pb[-1])
	print(seq_id+"\tSNAP\tmRNA\t"+mrna_st+"\t"+mrna_en+"\t.\t"+strand+"\t.\tID="+gene_id)
	for b in info_box:
		print(b)
TMP.close()
os.system("rm "+sys.argv[1]+".tMp")
