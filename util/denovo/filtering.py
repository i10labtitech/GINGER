import sys

args = sys.argv
total_data = []

threshold = float(args[2])

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        if line_trim and line_trim[0] == ">":
            tmp_data = [line_trim]
            total_data.append(tmp_data)
        else:
            total_data[-1].append(line_trim)
for entry in total_data:
    id = 0.0
    for item in entry:
        item_table = item.strip(" ").split(" ")
        if item_table[0] == "Percent" and item_table[1] == "identity:":
            id = float(item_table[2])
    if id > threshold:
        for item in entry:
            print(item)

#print(total_data)

