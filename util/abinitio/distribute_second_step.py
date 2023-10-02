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

Num=int(sys.argv[2])

rang=range(Num)
for i in rang:
	nUm=i+1
	TF=open("sPlitfile"+str(nUm)+".tMp","w")
	TF.close()
	del TF

del i

m=1
x=Num
lIst=open(sys.argv[1],"r")
for i in lIst:
	i = i.strip()
	ii = i.split("@")[1:]
	iii="".join(ii)	
	if m <= Num:
		TF=open("sPlitfile"+str(m)+".tMp","a")
		TF.write(iii+"\n")
		m+=1
		TF.close()
	elif m > Num:
		if x > 0:
			TF=open("sPlitfile"+str(x)+".tMp","a")
			TF.write(iii+"\n")
			x=x-1
			TF.close()
		elif x == 0:
			m = 1
			x = Num
			TF=open("sPlitfile"+str(m)+".tMp","a")
			TF.write(iii+"\n")
			TF.close()
			m+=1
		#m=1
		#TF=open("sPlitfile"+str(m)+".tMp","a")
		#TF.write(iii+"\n")
	del TF
lIst.close()

