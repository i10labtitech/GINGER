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

if len(args) < 3:
    print("Usage:")
    print("python " + args[0] + " Group.gff [threshold(%)] [one_class_min_score]")
    sys.exit()

threshold = float(args[2])
one_class_min = float(args[3]) if len(args) > 3 else 0.0   # >0: also rescue 1-class groups at/above this weighted score
group_dict = {}
counter = 1


def score_compare(group):
    global counter
    scaffold = group[0][0]
    score_dict = {}
    new_score_dict = {}
    for item in group:
        key = (item[1], item[3], item[4], item[6], item[7])   # tuple key: source tags may contain "_"
        score = float(item[5])
        score_dict.setdefault(key, 0.0)
        if score_dict[key] < score:
            score_dict[key] = score
    for key, value in score_dict.items():
        new_key = key[1:]
        new_score_dict.setdefault(new_key, 0.0)
        new_score_dict[new_key] += value
    max_score = 0.0
    max_key = None
    for key, value in new_score_dict.items():
        if max_score < value:
            max_score = value
            max_key = key
    max_items = list(max_key)
    print(scaffold + "\tphase2\tmRNA\t" + max_items[0] + "\t" + max_items[1] + "\t" + "{:.3f}".format(max_score) + "\t" + max_items[2] + "\t.\tID=single_exon_gene_" + str(counter))
    print(scaffold + "\tphase2\tCDS\t" + max_items[0] + "\t" + max_items[1] + "\t" + "{:.3f}".format(max_score) + "\t" + max_items[2] + "\t" + max_items[3] + "\tID=single_exon_gene_" + str(counter) + ".1;Parent=single_exon_gene_" + str(counter))
    counter += 1

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 11:
            gff_data = line_table[3:]
            if gff_data[2] == "CDS":
                group_dict.setdefault(line_table[0], []).append(gff_data)

for value in group_dict.values():
    filtered_group = []
    max_pos = 0
    min_pos = 0
    for item in value:
        if max_pos == 0 or max_pos < int(item[4]):
            max_pos = int(item[4])
        if min_pos == 0 or min_pos > int(item[3]):
            min_pos = int(item[3])
    group_len = max_pos - min_pos + 1
    for item in value:
        item_len = int(item[4]) - int(item[3]) + 1
        if (item_len / float(group_len)) * 100 > threshold:
            filtered_group.append(item)
    rna_flag = 0
    homology_flag = 0
    abinitio_flag = 0
    best_class_score = {"rna": 0.0, "hom": 0.0, "abi": 0.0}
    for item in filtered_group:
        method = item[1]
        score = float(item[5])
        if method == "mappingbase" or method == "denovobase" or method.startswith("RNASEQ_OTHER"):
            rna_flag = 1
            best_class_score["rna"] = max(best_class_score["rna"], score)
        elif method == "homology" or method.startswith("HOMOLOGY_OTHER"):
            homology_flag = 1
            best_class_score["hom"] = max(best_class_score["hom"], score)
        elif method == "AUGUSTUS" or method == "SNAP" or method.startswith("ABINITIO_OTHER"):
            abinitio_flag = 1
            best_class_score["abi"] = max(best_class_score["abi"], score)
    nclass = rna_flag + homology_flag + abinitio_flag
    if nclass >= 2:
        score_compare(filtered_group)
    elif nclass == 1 and one_class_min > 0 and max(best_class_score.values()) >= one_class_min:
        score_compare(filtered_group)   # strong single-class rescue (pseudogene-prone, hence the score bar)


