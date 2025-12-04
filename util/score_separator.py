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

score_set = []

with open (args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8 and line_table[0][0] != "#":
            if line_table[2] == "CDS":
                score_set.append(float(line_table[5]))

score_set.sort(reverse=True)

#sep_list = score_set[::len(score_set)/10]
#print(sep_list)

sep_list = score_set[::int(len(score_set)/10)][1:-1]
map_result = map(str, sep_list)
print("#" + " ".join(map_result))

with open (args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8 and line_table[0][0] != '#':
            if line_table[2] == "mRNA":
                print(line_trim)
            elif line_table[2] == "CDS":
                tmp_score = float(line_table[5])
                counter = 9
                for item in sep_list:
                    if tmp_score >= item:
                        line_table[5] = str(counter)
                        print("\t".join(line_table))
                        break
                    counter -= 1
                if counter == 0:
                    line_table[5] = "0"
                    print("\t".join(line_table))


