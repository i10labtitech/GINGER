import sys

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
