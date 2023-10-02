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
    print("python " + args[0] + " Group.gff")
    sys.exit()

def id_trim(table, tag):
    attribute = table[8].split(";")
    for item in attribute:
        key_value = item.split("=")
        if tag == key_value[0]:
            return key_value[1]

total_dict = {}

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 11:
            gff_data = line_table[3:]
            if gff_data[2] == "mRNA":
                tmp_id = id_trim(gff_data, "ID")
                total_dict.setdefault(tmp_id, [])
            elif gff_data[2] == "CDS":
                tmp_id = id_trim(gff_data, "Parent")
                total_dict[tmp_id].append(gff_data)

for key, value in total_dict.items():
    if len(value) == 1:
        print(value[0][0] + "\t" + value[0][1] + "\tmRNA\t" + value[0][3] + "\t"+ value[0][4] + "\t" + value[0][5] + "\t" + value[0][6] + "\t.\tID=" + key)
        print("\t".join(value[0]))


