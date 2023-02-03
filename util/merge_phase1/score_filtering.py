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
    for l in gff:
        ll = l.split("\t")
        if ll[2] == "mRNA":
            score = float(ll[5])
            info = ll[8]
            m = re.search(r"ID=group_num_([0-9]+)_",info)
            gro_num = m.group(1)
            out[gro_num] = score
            
    return out


if len(sys.argv) < 4:
    print ('score_filtering.py phase1.gff group.gff thr[FLOAT] out.gff')
    exit()

gff=open(sys.argv[1],"r")
group=open(sys.argv[2],"r")
wfile=open (sys.argv[4],"w")
thr=float(sys.argv[3])

gro2length = GroupLengthGet(group)
gro2score = GroupIdScoreGet(gff)
gff.seek(0)

for g in gff:
    if g[0] != "#":
        gg = g.split("\t")
        info = gg[8]
        m = re.search(r"ID=group_num_([0-9]+)_",info)
        gro_num = m.group(1)
        if gro2length[gro_num] * thr < gro2score[gro_num]:
            wfile.write(g)

gff.close()
group.close()
wfile.close()
