import sys

args = sys.argv

with open (args[1], "r") as data:
    for line in data:
        line_trim = line.strip("\n")
        line_table = line_trim.split("\t")
        if len(line_table) > 8 and line_table[0][0] != "#":
            if line_table[2] == "exon":
                tmp_tag = line_table[8].split(";")
                tmp_id = ""
                tmp_cov = ""
                tmp_exon = ""
                for item in tmp_tag:
                    if "exon_number" in item:
                        tmp_exon = item.split('\"')[1]
                    elif "cov" in item:
                        tmp_cov = item.split('\"')[1]
                    elif "transcript_id" in item:
                        tmp_id = item.split('\"')[1]
                print(tmp_id + "\t" + tmp_exon + "\t" + tmp_cov)
            
            
