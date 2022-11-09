import sys

args = sys.argv

if len(args) != 3:
    print("Usage:")
    print("python weighten.py [all.gff] [weight.txt]")
    sys.exit()

pre_weight_dict = {}
weight_dict = {}
max_weight = 0

with open(args[2], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        tmp_weight = float(line_table[1])
        pre_weight_dict[line_table[0]] = tmp_weight
        if max_weight < tmp_weight:
            max_weight = tmp_weight

for key, val in pre_weight_dict.items():
    weight_dict[key] = val / max_weight

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8 and line_table[0][0] != "#":
            if line_table[2] == "CDS":
                line_table[5] = str(float(line_table[5]) * weight_dict[line_table[1]])
            print("\t".join(line_table))
        else:
            print(line_trim)
                
