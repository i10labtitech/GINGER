import sys

args = sys.argv

if len(args) != 3:
    print("Usage:")
    print("python " + args[0] + " [gff] [seqkit fx2tab -nl output]")
    sys.exit()

total_dict = {}

with open(args[2], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        total_dict[line_table[0]] = int(line_table[1])

flag = False

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if line_table[2] == "gene":
            beg_pos = int(line_table[3])
            end_pos = int(line_table[4])
            if beg_pos > 0 and end_pos <= total_dict[line_table[0]]:
                flag = True
                print(line_trim)
            else:
                flag = False
        else:
            if flag:
                print(line_trim)

