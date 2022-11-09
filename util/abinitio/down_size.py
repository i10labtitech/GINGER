#coding:utf-8

import sys
import os

gb = open(sys.argv[1],"r")
gff = open(sys.argv[2],"r")

switch = 0
box = []

for b in gb:
	b = b.strip()
	if len(b)>5:
		if b[:3] == "CDS":
			switch = 1
		elif b[1:5] == "gene" and switch == 1:
			switch = 0
			bb = b.split('"')[1]
			if bb not in box:
				box.append(bb)
gb.close()

num = str(len(box))
res = open("down_size"+num+".gff","w")
err = open("err.log","w")
type_list = ["mRNA","CDS","exon"]
e = 0
ebox = []

for f in gff:
	f = f.strip()
	if f[0] != "#":
		ff = f.split("\t")
		kata = ff[2]
		idd = ff[-1]
		if kata in type_list:
			if kata == "mRNA":
				c1 = idd.split(";")[0].split("=")[1]
				if c1 in box:
					res.write(f+"\n")
			elif kata == "exon":
				c2 = idd.split(";")[1].split("=")[1]
				if c2 in box:
					res.write(f+"\n")
			elif kata == "CDS":
				c3 = idd.split(";")[1].split("=")[1]
				if c3 in box:
					res.write(f+"\n")
		elif kata not in type_list:
			if kata not in ebox: 
				e += 1
				err.write("["+str(e)+"]注意："+kata+"を検出しました。このプログラムではmRNA,CDS,exon以外扱うことができません。\n")
				ebox.append(kata)
			elif kata in ebox:
				e +=1
res.close()

if e == 0:
	err.write("エラーはありませんでした\n")
elif e > 0:
	err.write("\nエラーは合計"+num+"ありました\n")
err.close()

#注意：GFFの9カラム目がID、Parentの順に並んでいる必要あり。
