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

args = sys.argv

if len(args) != 3:
    print("Usage:")
    print("python " + args[0] + " [gff] [seqkit fx2tab -nl output]")
    sys.exit()

total_dict = {}

with open(args[2], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        total_dict[line_table[0]] = int(line_table[1])

flag = False

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[2] == "gene":
            beg_pos = int(line_table[3])
            end_pos = int(line_table[4])
            if beg_pos > 0 and end_pos <= total_dict[line_table[0]]:
                flag = True
                print(line_trim)
            else:
                flag = False
        else:
            if flag:
                print(line_trim)

