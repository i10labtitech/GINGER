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

def id_trim(values, name):
    records = values[8].split(";")
    for item in records:
        key_val = item.split("=")
        if key_val[0] == name:
            return key_val[1]

def translate(seq):
    answer = ""
    for i in range(len(seq)//3):
        codon = seq[i*3:i*3+3]
        if codon[0] == "A":
            if codon[1] == "A":
                if codon[2] == "A" or codon[2] == "G":
                    answer += "K"
                elif codon[2] == "C" or codon[2] == "T":
                    answer += "N"
                else:
                    answer += "X"
            elif codon[1] == "C":
                answer += "T"
            elif codon[1] == "G":
                if codon[2] == "A" or codon[2] == "G":
                    answer += "R"
                elif codon[2] == "C" or codon[2] == "T":
                    answer += "S"
                else:
                    answer += "X"
            elif codon[1] == "T":
                if codon[2] == "A" or codon[2] == "C" or codon[2] == "T":
                    answer += "I"
                elif codon[2] == "G":
                    answer += "M"
                else:
                    answer += "X"
            else:
                answer += "X"
        elif codon[0] == "C":
            if codon[1] == "A":
                if codon[2] == "A" or codon[2] == "G":
                    answer += "Q"
                elif codon[2] == "C" or codon[2] == "T":
                    answer += "H"
                else:
                    answer += "X"
            elif codon[1] == "C":
                answer += "P"
            elif codon[1] == "G":
                answer += "R"
            elif codon[1] == "T":
                answer += "L"
            else:
                answer += "X"
        elif codon[0] == "G":
            if codon[1] == "A":
                if codon[2] == "A" or codon[2] == "G":
                    answer += "E"
                elif codon[2] == "C" or codon[2] == "T":
                    answer += "D"
                else:
                    answer += "X"
            elif codon[1] == "C":
                answer += "A"
            elif codon[1] == "G":
                answer += "G"
            elif codon[1] == "T":
                answer += "V"
            else:
                answer += "X"
        elif codon[0] == "T":
            if codon[1] == "A":
                if codon[2] == "A" or codon[2] == "G":
                    answer += "*"
                elif codon[2] == "C" or codon[2] == "T":
                    answer += "Y"
                else:
                    answer += "*"
            elif codon[1] == "C":
                answer += "S"
            elif codon[1] == "G":
                if codon[2] == "A":
                    answer += "*"
                elif codon[2] == "C" or codon[2] == "T":
                    answer += "C"
                elif codon[2] == "G":
                    answer += "W"
                else:
                    answer += "X"
            elif codon[1] == "T":
                if codon[2] == "A" or codon[2] == "G":
                    answer += "L"
                elif codon[2] == "C" or codon[2] == "T":
                    answer += "F"
                else:
                    answer += "X"
            else:
                answer += "X"
        else:
            answer += "X"
    return answer

if len(args) != 4:
    print("Usage:")
    print("python " + args[0] + " [annotation.gff] [genome.fa] [output prefix]")
    sys.exit()

genome_dict = {} # fasta_id : seq
tmp_id = ""
tmp_seq = ""

with open(args[2], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        if line_trim[0] == ">":
            if tmp_seq != "":
                genome_dict[tmp_id] = tmp_seq
            tmp_id = line_trim[1:]
            tmp_seq = ""
        else:
            tmp_seq += line_trim
    if tmp_seq != "":
        genome_dict[tmp_id] = tmp_seq
    tmp_seq = ""

mrna_id = ""
cds_data = []

id_w_data = []
cds_w_data = []
pep_w_data = []

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[2] == "mRNA":
            if cds_data != []:
                if cds_data[0][3] == "+":
                    cds_data.sort(key=itemgetter(1),reverse=False)
                    for item in cds_data:
                        tmp_seq += genome_dict[item[0]][item[1]-1:item[2]]
                else:
                    cds_data.sort(key=itemgetter(2), reverse=True)
                    #print(cds_data)
                    for item in cds_data:
                        tmp_seq += genome_dict[item[0]][item[1]-1:item[2]][::-1]
                tmp_seq = tmp_seq[cds_data[0][4]:]
                id_w_data.append(">" + mrna_id)
                cds_w_data.append(tmp_seq)
                pep_w_data.append(translate(tmp_seq))
            mrna_id = id_trim(line_table, "ID")
            cds_data = []
            tmp_seq = ""
        if line_table[2] == "CDS":
            cds_data.append([line_table[0], int(line_table[3]), int(line_table[4]), line_table[6], int(line_table[7])])
    if cds_data != []:
            if cds_data[0][3] == "+":
                cds_data.sort(key=itemgetter(1),reverse=False)
                for item in cds_data:
                    tmp_seq += genome_dict[item[0]][item[1]-1:item[2]]
            else:
                cds_data.sort(key=itemgetter(2), reverse=True)
                for item in cds_data:
                    tmp_seq += genome_dict[item[0]][item[1]-1:item[2]][::-1]
            tmp_seq = tmp_seq[cds_data[0][4]:]
            id_w_data.append(">" + mrna_id)
            cds_w_data.append(tmp_seq)
            pep_w_data.append(translate(tmp_seq))

l_width = 60
with open(args[3]+".cds", "w") as data:
    for n, m in zip(id_w_data, cds_w_data):
        data.write(n+"\n")
        for i in range((len(m)//l_width)+1):
            data.write(m[i*l_width:i*l_width+l_width]+"\n")

with open(args[3]+".pep", "w") as data:
    for n, m in zip(id_w_data, pep_w_data):
        data.write(n+"\n")
        for i in range((len(m)//l_width)+1):
            data.write(m[i*l_width:i*l_width+l_width]+"\n")





