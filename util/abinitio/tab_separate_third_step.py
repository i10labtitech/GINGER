#coding:utf-8

import sys

File=open(sys.argv[1],"r")

for i in File:
	i = i.strip()
	ii = i.split("\t")
	print(ii[0])
	print(ii[1])
File.close()

