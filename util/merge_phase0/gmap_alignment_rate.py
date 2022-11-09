import sys

args = sys.argv

if len(args) != 2:
    print("Usage:")
    print("python gmap_alignment_rate.py [gmap_result.gmap]")
    sys.exit()

total_dict = {}

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        tmp_id = ""
        if line_trim:
            if line_trim[0] == ">":
                tmp_ID = line_trim[1:]
                total_dict[tmp_ID] = []
            else:
                line_trim_2 = line_trim.strip(" ")
                if line_trim_2[0] == "+" or line_trim_2[0] == "-":
                    line_table = line_trim_2.split()
                    total_dict[tmp_ID].append(line_table[2][:-1])

for key, val in total_dict.items():
    for i in range(len(val)):
        print(key + "\t" + str(i+1) + "\t" + val[i])
