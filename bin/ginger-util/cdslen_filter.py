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
    print("python " + args[0] + " [phase2_output.gff] [threshold(bp)]")
    sys.exit()

threshold = int(args[2])
tmp_mrna = ""
tmp_cds = []
cds_len = 0

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[2] == "mRNA":
            if tmp_cds != []:
                if threshold < cds_len:
                    print(tmp_mrna)
                    for item in tmp_cds:
                        print(item)
                tmp_cds = []
                cds_len = 0
            tmp_mrna = line_trim
        elif line_table[2] == "CDS":
            tmp_cds.append(line_trim)
            cds_len += (int(line_table[4]) - int(line_table[3]) + 1)

