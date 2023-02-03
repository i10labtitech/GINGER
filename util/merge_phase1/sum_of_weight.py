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

methods = {
    "MAPPING_WEIGHT":0,
    "DENOVO_WEIGHT":0,
    "HOMOLOGY_WEIGHT":0,
    "AUGUSTUS_WEIGHT":0,
    "SNAP_WEIGHT":0
}

max_weight = 0

with open (args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("=")
        if line_table[0] in methods:
            methods[line_table[0]] += float(line_table[1])
            if max_weight < float(line_table[1]):
                max_weight = float(line_table[1])
sum_of_score = 0

for val in methods.values():
    sum_of_score += val / max_weight

print("{:.2f}".format(sum_of_score))
