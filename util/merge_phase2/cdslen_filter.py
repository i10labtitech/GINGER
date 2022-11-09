import sys

args = sys.argv

if len(args) != 3:
    print("Usage:")
    print("python " + args[0] + " [phase2_output.gff] [threshold(bp)]")
    sys.exit()

threshold = int(args[2])
tmp_mrna = ""
tmp_cds = []
cds_len = 0

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[2] == "mRNA":
            if tmp_cds != []:
                if threshold < cds_len:
                    print(tmp_mrna)
                    for item in tmp_cds:
                        print(item)
                tmp_cds = []
                cds_len = 0
            tmp_mrna = line_trim
        elif line_table[2] == "CDS":
            tmp_cds.append(line_trim)
            cds_len += (int(line_table[4]) - int(line_table[3]) + 1)

