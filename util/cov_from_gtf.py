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

with open (args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8 and line_table[0][0] != "#":
            if line_table[2] == "exon":
                tmp_tag = line_table[8].split(";")
                tmp_id = ""
                tmp_cov = ""
                tmp_exon = ""
                for item in tmp_tag:
                    if "exon_number" in item:
                        tmp_exon = item.split('\"')[1]
                    elif "cov" in item:
                        tmp_cov = item.split('\"')[1]
                    elif "transcript_id" in item:
                        tmp_id = item.split('\"')[1]
                print(tmp_id + "\t" + tmp_exon + "\t" + tmp_cov)
            
            
