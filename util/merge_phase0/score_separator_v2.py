import sys

args = sys.argv

score_set = []

with open (args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8 and line_table[0][0] != "#":
            if line_table[2] == "CDS":
                score_set.append(float(line_table[5]))

score_set.sort(reverse=True)

#sep_list = score_set[::len(score_set)/10]
#print(sep_list)

parser = int(len(score_set) / 10)
sep_list = []
for i in range(9):
    sep_list.append(score_set[(i+1) * parser])
map_result = map(str, sep_list)
print("#" + " ".join(map_result))

with open (args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8 and line_table[0][0] != '#':
            if line_table[2] == "mRNA":
                print(line_trim)
            elif line_table[2] == "CDS":
                tmp_score = float(line_table[5])
                counter = 9
                for item in sep_list:
                    if tmp_score >= item:
                        line_table[5] = str(counter)
                        print("\t".join(line_table))
                        break
                    counter -= 1
                if counter == 0:
                    line_table[5] = "0"
                    print("\t".join(line_table))


