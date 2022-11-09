import sys

args = sys.argv

if len(args) != 2:
    print("Usage:")
    print("python " + args[0] + " [gff]")
    sys.exit()

def id_trim(values, name):
    records = values[11].split(";")
    for item in records:
        key_val = item.split("=")
        if key_val[0] == name:
            return key_val[1]

mrna_data = []
cds_dict = {}

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[5] == "mRNA":
            tmp_id = id_trim(line_table, "ID")
            mrna_data.append([tmp_id, line_table])
        elif line_table[5] == "CDS":
            tmp_id = id_trim(line_table, "Parent")
            cds_dict.setdefault(tmp_id, []).append(line_table)

for item in mrna_data:
    s_pos = [int(x[6]) for x in cds_dict[item[0]]]
    e_pos = [int(x[7]) for x in cds_dict[item[0]]]
    ret = item[1]
    ret[6] = str(min(s_pos))
    ret[7] = str(max(e_pos))
    print("\t".join(ret))
    for i in cds_dict[item[0]]:
        print("\t".join(i))
    

