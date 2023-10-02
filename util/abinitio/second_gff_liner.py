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
