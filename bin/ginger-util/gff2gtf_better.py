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
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.i

import sys

infile = open(sys.argv[1],"r")
errlog = open("gff2gtf.err.log","w")
#outfile_name = sys.argv[1].split(".")[0] + ".gtf"
#outfile = open(outfile_name,"w")
s = 0
t0 = ""
t1 = ""
t2 = ""
t3 = ""
t4 = ""

for i in infile:
	i = i.strip()
	if len(i) > 0:
		if i[0] == "#":
			#outfile.write(i+"\n")
			print(i)
		elif i[0] != "#":
			ii = i.split("\t")
			if len(ii) > 8:
				scaffold_name = ii[0]
				tool_name = ii[1]
				data_type = ii[2]
				start_pos = ii[3]
				end_pos = ii[4]
				u1 = ii[5]
				strand = ii[6]
				u2 = ii[7]
				tag_name = ii[8]
				if data_type == "mRNA":
					t0 = tag_name.split(";")[0].split("=")[1]
					#t0 = tag_name.split(";")[0].split(":")[1]
					#t1 = tag_name.split(";")[1].split(":")[1]
					t2 = 'gene_id "'+t0+'"; transcript_id  "'+t0+'"'
					#outfile.write(ii[0]+"\t"+ii[1]+"\t"+ii[2]+"\t"+ii[3]+"\t"+ii[4]+"\t"+ii[5]+"\t"+ii[6]+"\t"+ii[7]+"\t"+t2+"\n")
					print(ii[0]+"\t"+ii[1]+"\t"+ii[2]+"\t"+ii[3]+"\t"+ii[4]+"\t"+ii[5]+"\t"+ii[6]+"\t"+ii[7]+"\t"+t2)
				if data_type != "mRNA":
					t3 = tag_name.split(";")[1].split("=")[1]
					t4 = 'gene_id "'+t3+'"; transcript_id  "'+t3+'"'
					#outfile.write(ii[0]+"\t"+ii[1]+"\t"+ii[2]+"\t"+ii[3]+"\t"+ii[4]+"\t"+ii[5]+"\t"+ii[6]+"\t"+ii[7]+"\t"+t4+"\n")
					print(ii[0]+"\t"+ii[1]+"\t"+ii[2]+"\t"+ii[3]+"\t"+ii[4]+"\t"+ii[5]+"\t"+ii[6]+"\t"+ii[7]+"\t"+t4)
			if len(ii) <= 8:
				s += 1
				errlog.write("Error"+str(s)+": seems to be a row with a small number of columns\n")
if s == 0:
	errlog.write("Finished\n")
elif s > 0:
	errlog.write("Total "+str(s)+" erros found\n")


