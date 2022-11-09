import sys

args = sys.argv

total_dict = {}

#print(args[1], file=sys.stderr); # tany
#print(args[2], file=sys.stderr); # tany

#coverage_stats
with open (args[2], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
#        print(line_table, file=sys.stderr); # tany 
        total_dict.setdefault(line_table[0], []).append(line_table[2])

#GFF
with open (args[1], "r") as data:
    exon_num_counter = 0
    for line in data:
        line_trim = line.strip("\n")
        line_table = line.split("\t")
        if len(line_table) > 8 and line_table[0][0] != "#":
            if line_table[2] == "mRNA":
                print(line_trim)
                exon_num_counter = -1
            if line_table[2] == "exon":
                exon_num_counter += 1
            if line_table[2] == "CDS":
                tag_table = line_table[8].split(";")
                for tag in tag_table:
                    if "Parent" in tag:
                        tmp_ID = tag.split("=")[1].strip("\n")
                        tmp_ID_alt = ".".join(tmp_ID.split(".")[:-1])
#                        if not tmp_ID_alt in total_dict:
#                            print("BBB0 %s %s" %(tmp_ID_alt, args[1]), file=sys.stderr) # tany 20220522
#                        if len(total_dict[tmp_ID_alt]) <= exon_num_counter:
#                            print("BBB1 %s %d %d %s" %(tmp_ID_alt, len(total_dict[tmp_ID_alt]), exon_num_counter, args[1]), file=sys.stderr) # tany 20220522
                        line_table[5] = total_dict[tmp_ID_alt][exon_num_counter]
                print("\t".join(line_table))
#                print("\t".join(line_table), file=sys.stderr) # tany


