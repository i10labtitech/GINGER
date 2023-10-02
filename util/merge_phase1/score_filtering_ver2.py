#coding:utf-8

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
import os
import re

def GroupLengthGet(group):
    count = 0
    stposlist = []
    edposlist = []
    out = {}
    for l in group:
        count = count + 1
        ll = l.split("\t")
        if count == 1:
            gro_num = ll[0]
        if gro_num != ll[0]:
            stposlist.sort()
            edposlist.sort(reverse=True)
            length = edposlist[0] - stposlist[0] + 1
            out[gro_num] = length
            stposlist = []
            edposlist = []
        if ll[5] == "CDS":
            stposlist.append(int(ll[6]))
            edposlist.append(int(ll[7]))
        gro_num = ll[0]
    stposlist.sort()
    edposlist.sort(reverse=True)
    length = edposlist[0] - stposlist[0] + 1
    out[gro_num] = length

    return out

def GroupIdScoreGet(gff):
    out = {}
    cds = {}
    for l in gff:
        ll = l.split("\t")
        if ll[2] == "mRNA":
            score = float(ll[5])
            info = ll[8]
            m = re.search(r"ID=group_num_([0-9]+)_",info)
            gro_num = m.group(1)
            out[gro_num] = score
        elif ll[2] == "CDS":
            info = ll[8]
            m = re.search(r"ID=group_num_([0-9]+_[0-9]+)\.mrna[0-9]+\.cds([0-9]+)",info)
            #ID=group_num_1_1.mrna1.cds4;Parent=group_num_1_1.mrna1
            mRNA_num = m.group(1)
            cds_count = int(m.group(2))
            if (mRNA_num in cds) and (cds[mRNA_num] < cds_count):
                cds[mRNA_num] = cds_count
            else:
                cds[mRNA_num] = 1
            
    return out, cds


if len(sys.argv) < 4:
    print ('score_filtering.py phase1.gff group.gff thr1[FLOAT](for single exon) the2[FLOAT](multi exon) out.gff')
    exit()

gff=open(sys.argv[1],"r")
group=open(sys.argv[2],"r")
wfile=open (sys.argv[5],"w")
thr1=float(sys.argv[3])
thr2=float(sys.argv[4])

print("Threshold:",thr1,thr2) 

gro2length = GroupLengthGet(group)
gro2score, cdsCount = GroupIdScoreGet(gff)
gff.seek(0)

for g in gff:
    if g[0] != "#":
        gg = g.split("\t")
        info = gg[8]
        m = re.search(r"ID=group_num_([0-9]+)_([0-9]+)",info)
        gro_num = m.group(1)
        sgro_num = m.group(2)
        mRNA_num = gro_num + "_" + sgro_num
        if cdsCount[mRNA_num] == 1:
            if gro2length[gro_num] * thr1 < gro2score[gro_num]:
                #print("1", gro2length[gro_num], thr1, gro2length[gro_num] * thr1, gro2score[gro_num])
                #print(g)
                wfile.write(g)
        else:
            if gro2length[gro_num] * thr2 < gro2score[gro_num]:
                #print("2", gro2length[gro_num], thr2, gro2score[gro_num])
                #print(g)
                wfile.write(g)

gff.close()
group.close()
wfile.close()
