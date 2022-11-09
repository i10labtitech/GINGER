#coding:utf-8

import sys

line_file=open(sys.argv[1],"r")
box=[]

for i in line_file:
	i = i.strip()
	if i not in box:
		print(i)
		box.append(i)
line_file.close()

