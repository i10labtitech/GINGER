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

total_data = []
tmp_data = []
tmp_id = ""

def scoring(item):
    longlen = max([(item[2] - item[1] + 1), (item[5] - item[4] + 1)*3])
    shortlen = min([(item[2] - item[1] + 1), (item[5] - item[4] + 1)*3])
    mismatch = item[6]*3
    return float(shortlen - mismatch) * 100 / float(longlen)

if len(args) != 3:
    print("Usage:")
    print("python o4_to_gff.py [spalnresult_o4] [ID prefix]")
    sys.exit()

prefix = args[2]

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        if not (line_trim[0] == "#" or line_trim[0] == "@"):
            line_table = line_trim.split("\t")
            beg_pos = int(line_table[8])
            end_pos = int(line_table[9])
            q_beg = int(line_table[6])
            q_end = int(line_table[7])
            MisMch = int(line_table[4])
            Indels = int(line_table[5])
            if beg_pos < end_pos :
                #scaffold, begin_pos, end_pos, strand, query_begin_pos, query_end_pos, mismatch, indels, query_ID
                tmp_data.append([line_table[1], beg_pos, end_pos, "+", q_beg, q_end, MisMch, Indels, line_table[0]])
            else :
                tmp_data.append([line_table[1], end_pos, beg_pos, "-", q_beg, q_end, MisMch, Indels, line_table[0]])
            #total_data.setdefault(line_table[0],[]).append([float(line_table[2]), int(line_table[3])])
        else:
            if tmp_data:
                total_data.append(tmp_data)
            tmp_data = []

counter = 1
c_counter = 1
frame = 0
for item in total_data:
    scaffold = item[0][0]
    strand = item[0][3]
    source = item[0][8]
    mrna_beg, mrna_end = 0, 0
    if strand == "+":
        item.sort(key=itemgetter(1))
        mrna_beg = item[0][1]
        mrna_end = item[-1][2] + 3
        item[-1][2] += 3 #stop codon
        item[0][5] += 3
    else:
        item.sort(key=itemgetter(1), reverse=True)
        mrna_beg = item[-1][1] - 3
        mrna_end = item[0][2]
        item[-1][1] -= 3 #stop codon
        item[0][5] += 3
    if mrna_beg == 0 or mrna_end == 0:
        continue
    gene = [scaffold, "ALN", "gene", str(mrna_beg), str(mrna_end), ".", strand, ".", "ID=gene"+str(counter)+"."+prefix+";Name=nil;Source="+source]
    mrna = [scaffold, "ALN", "mRNA", str(mrna_beg), str(mrna_end), ".", strand, ".", "ID=mRNA"+str(counter)+"."+prefix+";Parent=gene"+str(counter)+"."+prefix+";Name=nil"]
    print("\t".join(gene))
    print("\t".join(mrna))
    for i in item:
        tmp_score = scoring(i)
        cds = [scaffold, "ALN", "cds", str(i[1]), str(i[2]), "{0:.3f}".format(tmp_score), strand, str(frame), "ID=cds"+str(c_counter)+";Parent=mRNA"+str(counter)+"."+prefix+";Name=nil"]
        print("\t".join(cds))
        frame = (frame - (i[2] - i[1] + 1) % 3) % 3
        c_counter += 1
    counter += 1
    frame = 0
    
