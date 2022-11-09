import sys

args = sys.argv

def id_trim(values, name):
    records = values[8].split(";")
    for item in records:
        key_val = item.split("=")
        if key_val[0] == name:
            return key_val[1]

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8:
            if line_table[2] == "mRNA":
                print(id_trim(line_table, "normalized_score"))

