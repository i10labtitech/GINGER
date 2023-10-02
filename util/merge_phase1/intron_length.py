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
from operator import itemgetter

args = sys.argv

def calc_len(data):
    data.sort(key=itemgetter(0))
    for idx, item in enumerate(data):
        if idx != 0:
            print(data[idx][0] - data[idx-1][1])

if len(args) != 4:
    print("Usage:")
    print("python " + args[0] + " [gff] [target] [delimiter]")
    sys.exit()

total_dict = {}
pos_set = []

tgt = args[2]
delim = args[3]

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8:
            if line_table[2] == tgt:
                if pos_set != []:
                    calc_len(pos_set)
                pos_set = []
            elif line_table[2] == delim:
                s_pos = int(line_table[3])
                e_pos = int(line_table[4])
                pos_set.append([s_pos, e_pos])
if pos_set != []:
    calc_len(pos_set)

