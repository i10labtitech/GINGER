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

infile = open(sys.argv[1],'r')
name = sys.argv[2]
s = 0

rmk = open('rmk.fasta','r')
for o in rmk:
	o = o.strip()
	oo = o.strip(">")
	if o[0] == '>' and name == oo:
		print(o)
		s = 1
	elif '>' not in o and s == 1:
		print(o)
		s = 0
rmk.close()
