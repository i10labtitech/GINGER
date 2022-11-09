import sys

args = sys.argv

methods = {
    "MAPPING_WEIGHT":0,
    "DENOVO_WEIGHT":0,
    "HOMOLOGY_WEIGHT":0,
    "AUGUSTUS_WEIGHT":0,
    "SNAP_WEIGHT":0,
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
