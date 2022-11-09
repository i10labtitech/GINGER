#coding:utf-8

import sys
import os

gff = open(sys.argv[1],"r")

s = 0

for i in gff:
	i = i.strip()
	if i[0] != "#":
		ii = i.split("\t")
		type = ii[2]
		if type == "mRNA" and s == 0:
			s += 1
			sys.stdout.write(i)
		elif type == "mRNA" and s > 0:
			s += 1
			sys.stdout.write("\n"+i)
		else:
			sys.stdout.write("KIREME"+i)
sys.stdout.write("\n")
