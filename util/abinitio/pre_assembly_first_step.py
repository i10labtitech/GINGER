#coding:utf-8

import os
import sys

File=open(sys.argv[1],"r")
Num=int(sys.argv[2])

#fastaを一行にする

n=0
rmk=open("rmk.tMp","w")

for i in File:
	i = i.strip()
	if len(i) > 0:
		if i[0] == ">" and n == 0:
			n = 1
			rmk.write(i+"\t")
		elif i[0] == ">" and n > 0:
			n += 1
			rmk.write("\n"+i+"\t")
		elif i[0] != ">":
			rmk.write(i)
rmk.write("\n")
rmk.close()

del i
del n

#行頭に長さをつける。

rmk2=open("rmk.tMp","r")
List=open("list.tMp","w")

for i in rmk2:
	i = i.strip()
	ii = i.split("\t")[1]
	ln = str(len(ii))
	List.write(ln+"@"+i+"\n")
rmk2.close()

#出来上がったlist.tMpをソート
#sort -nrk 1 -t "@" list.tMp > list_sort.tMp

