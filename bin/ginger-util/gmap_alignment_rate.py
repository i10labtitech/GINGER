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

if len(args) != 2:
    print("Usage:")
    print("python gmap_alignment_rate.py [gmap_result.gmap]")
    sys.exit()

total_dict = {}

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        tmp_id = ""
        if line_trim:
            if line_trim[0] == ">":
                tmp_ID = line_trim[1:]
                total_dict[tmp_ID] = []
            else:
                line_trim_2 = line_trim.strip(" ")
                if line_trim_2[0] == "+" or line_trim_2[0] == "-":
                    line_table = line_trim_2.split()
                    total_dict[tmp_ID].append(line_table[2][:-1])

for key, val in total_dict.items():
    for i in range(len(val)):
        print(key + "\t" + str(i+1) + "\t" + val[i])
