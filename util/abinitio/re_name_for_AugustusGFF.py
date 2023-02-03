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
