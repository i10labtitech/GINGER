import sys
from operator import itemgetter

args = sys.argv

def calc_len(data):
    data.sort(key=itemgetter(0))
    for idx, item in enumerate(data):
        if idx != 0:
            print(data[idx][0] - data[idx-1][1])

if len(args) != 4:
    print("Usage:")
    print("python " + args[0] + " [gff] [target] [delimiter]")
    sys.exit()

total_dict = {}
pos_set = []

tgt = args[2]
delim = args[3]

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8:
            if line_table[2] == tgt:
                if pos_set != []:
                    calc_len(pos_set)
                pos_set = []
            elif line_table[2] == delim:
                s_pos = int(line_table[3])
                e_pos = int(line_table[4])
                pos_set.append([s_pos, e_pos])
if pos_set != []:
    calc_len(pos_set)

