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
	if len(f) < 1:
                continue
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
				err.write("["+str(e)+"]warning!："+kata+" is detected\n")
				ebox.append(kata)
			elif kata in ebox:
				e +=1
res.close()

if e == 0:
	err.write("No error\n")
elif e > 0:
	err.write("\n"+num+" errors\n")
err.close()
