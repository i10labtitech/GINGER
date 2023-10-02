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

total_dict = {}
threshold = 30

def id_trim(table, name):
    tags = table[8].split(";")
    for item in tags:
        key_value = item.split("=")
        if key_value[0] == name:
            return key_value[1]

if len(args) != 2:
    print("Usage:")
    print("python " + args[0] + " [all_spalnresult.gff]")
    sys.exit()

#filtering
with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[2] == "mRNA":
            tmp_id = id_trim(line_table, "ID")[4:]
            total_dict[tmp_id] = [0, 0]
        if line_table[2] == "cds":
            tmp_id = id_trim(line_table, "Parent")[4:]
            tmp_len = float(line_table[4]) - float(line_table[3])
            tmp_score = float(line_table[5])
            total_dict[tmp_id][0] += tmp_len
            total_dict[tmp_id][1] += tmp_len * tmp_score

#scoring
with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[2] == "gene":
            tmp_id = id_trim(line_table, "ID")[4:]
            if total_dict[tmp_id][1] / total_dict[tmp_id][0] > threshold:
                print(line_trim)
        elif line_table[2] == "cds":
            tmp_id = id_trim(line_table, "Parent")[4:]
            if total_dict[tmp_id][1] / total_dict[tmp_id][0] > threshold:
                line_table[5] = "{0:.3f}".format(float(line_table[5]) / 100)
                print("\t".join(line_table))

