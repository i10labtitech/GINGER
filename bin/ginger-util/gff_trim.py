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

def id_trim(table, name):
    values = line_table[8].split(";")
    for i in values:
        feat = i.split("=")
        if feat[0] == name:
            return feat[1]

if len(args) != 3:
    print("Usage:")
    print("python " + args[0] + " [gff] [list]")
    sys.exit()

total_dict = {}

with open(args[2], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        total_dict[line_trim] = 0

flag = False
with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[2] == "gene":
            tmp_ID = id_trim(line_table, "ID")
            if tmp_ID in total_dict:
                flag = False
            else:
                print(line_trim)
                flag = True
        else:
            if flag:
                print(line_trim)
