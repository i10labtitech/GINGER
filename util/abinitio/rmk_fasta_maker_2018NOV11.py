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
rmkfile = open('rmk.fasta','w')
switch = 0

for i in infile:
	i = i.strip()
	if len(i) > 0:
		if i[0] == '>' and switch == 0:
			rmkfile.write(i)
			rmkfile.write('\n')
			switch = 1
		elif i[0] == '>' and switch == 1:
			rmkfile.write('\n' + i)
			rmkfile.write('\n')
		else:
			rmkfile.write(i)
rmkfile.close()
