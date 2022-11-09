import sys

args = sys.argv

if len(args) != 2:
    print("Usage:")
    print("python " + args[0] + " Group.gff")
    sys.exit()

def id_trim(table, tag):
    attribute = table[8].split(";")
    for item in attribute:
        key_value = item.split("=")
        if tag == key_value[0]:
            return key_value[1]

total_dict = {}

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 11:
            gff_data = line_table[3:]
            if gff_data[2] == "mRNA":
                tmp_id = id_trim(gff_data, "ID")
                total_dict.setdefault(tmp_id, [])
            elif gff_data[2] == "CDS":
                tmp_id = id_trim(gff_data, "Parent")
                total_dict[tmp_id].append(gff_data)

for key, value in total_dict.items():
    if len(value) == 1:
        print(value[0][0] + "\t" + value[0][1] + "\tmRNA\t" + value[0][3] + "\t"+ value[0][4] + "\t" + value[0][5] + "\t" + value[0][6] + "\t.\tID=" + key)
        print("\t".join(value[0]))


