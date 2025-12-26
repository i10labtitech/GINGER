import sys

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

args = sys.argv
total_data = []

threshold = float(args[2])

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        if line_trim and line_trim[0] == ">":
            tmp_data = [line_trim]
            total_data.append(tmp_data)
        else:
            total_data[-1].append(line_trim)
for entry in total_data:
    id = 0.0
    for item in entry:
        item_table = item.strip(" ").split(" ")
        if item_table[0] == "Percent" and item_table[1] == "identity:":
            id = float(item_table[2])
    if id > threshold:
        for item in entry:
            print(item)

#print(total_data)

