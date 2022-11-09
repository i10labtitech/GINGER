#coding:utf-8

import sys
import os

gff = open(sys.argv[1],"r")
name = sys.argv[2]
#n_f = open(sys.argv[2]+".tmp.gff","w")

for i in gff:
	i = i.strip()
	if len(i) > 0 and i[0] != "#":
		ii = i.split("\t")
		seq_name=ii[0]
		if seq_name == name:
			#print(seq_name+"\t"+name)
			print(i)

gff.close()


##説明
#ann作成第一段階
#対象の配列に予測された遺伝子情報のみをgffまたはgtfから抜き出す。
