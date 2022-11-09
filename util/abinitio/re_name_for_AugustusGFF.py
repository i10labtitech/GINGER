#coding:utf-8

import sys

File=open(sys.argv[1],"r")
n=0
name1=""
name2=""


for i in File:
	i = i.strip()
	ii = i.split("\t")
	tag = ii[2]
	data = ii[8]
	if tag == "gene":
		n+=1
		name1 = "g" + str(n)
		iii = "\t".join(ii[:8])
		print(iii+"\t"+name1)
	elif tag == "transcript":
		name2 = "g"+str(n)+".t1"#+str(n)
		iii = "\t".join(ii[:8])
		print(iii+"\t"+name2)
	else:
		iii = "\t".join(ii[:8])
		print(iii+'\ttranscript_id "'+name2+'"; gene_id "'+name1+'";')
File.close()
