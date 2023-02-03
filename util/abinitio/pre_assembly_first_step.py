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

import os
import sys

File=open(sys.argv[1],"r")
Num=int(sys.argv[2])

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

rmk2=open("rmk.tMp","r")
List=open("list.tMp","w")

for i in rmk2:
	i = i.strip()
	ii = i.split("\t")[1]
	ln = str(len(ii))
	List.write(ln+"@"+i+"\n")
rmk2.close()
