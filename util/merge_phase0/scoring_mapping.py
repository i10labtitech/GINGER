import sys
import math

args = sys.argv

with open(args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8 and line_table[0][0] != "#":
            if line_table[2] == "CDS":
                tmp_score = float(line_table[5])
                if tmp_score >= 100:
                    line_table[5] = "1.0"
                elif tmp_score <= 1:
                    line_table[5] = "0.7"
                else:
                    line_table[5] = "{0:.3f}".format(math.log10(tmp_score) * 0.15 + 0.7)
            print("\t".join(line_table))
        else:
            print(line_trim)
