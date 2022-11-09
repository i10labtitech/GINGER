#coding:utf-8

import sys

infile = open(sys.argv[1],"r")
errlog = open("gff2gtf.err.log","w")
#outfile_name = sys.argv[1].split(".")[0] + ".gtf"
#outfile = open(outfile_name,"w")
s = 0
t0 = ""
t1 = ""
t2 = ""
t3 = ""
t4 = ""

for i in infile:
	i = i.strip()
	if len(i) > 0:
		if i[0] == "#":
			#outfile.write(i+"\n")
			print(i)
		elif i[0] != "#":
			ii = i.split("\t")
			if len(ii) > 8:
				scaffold_name = ii[0]
				tool_name = ii[1]
				data_type = ii[2]
				start_pos = ii[3]
				end_pos = ii[4]
				u1 = ii[5]
				strand = ii[6]
				u2 = ii[7]
				tag_name = ii[8]
				if data_type == "mRNA":
					t0 = tag_name.split(";")[0].split("=")[1]
					#t0 = tag_name.split(";")[0].split(":")[1]
					#t1 = tag_name.split(";")[1].split(":")[1]
					t2 = 'gene_id "'+t0+'"; transcript_id  "'+t0+'"'
					#outfile.write(ii[0]+"\t"+ii[1]+"\t"+ii[2]+"\t"+ii[3]+"\t"+ii[4]+"\t"+ii[5]+"\t"+ii[6]+"\t"+ii[7]+"\t"+t2+"\n")
					print(ii[0]+"\t"+ii[1]+"\t"+ii[2]+"\t"+ii[3]+"\t"+ii[4]+"\t"+ii[5]+"\t"+ii[6]+"\t"+ii[7]+"\t"+t2)
				if data_type != "mRNA":
					t3 = tag_name.split(";")[1].split("=")[1]
					t4 = 'gene_id "'+t3+'"; transcript_id  "'+t3+'"'
					#outfile.write(ii[0]+"\t"+ii[1]+"\t"+ii[2]+"\t"+ii[3]+"\t"+ii[4]+"\t"+ii[5]+"\t"+ii[6]+"\t"+ii[7]+"\t"+t4+"\n")
					print(ii[0]+"\t"+ii[1]+"\t"+ii[2]+"\t"+ii[3]+"\t"+ii[4]+"\t"+ii[5]+"\t"+ii[6]+"\t"+ii[7]+"\t"+t4)
			if len(ii) <= 8:
				s += 1
				errlog.write("Error"+str(s)+"：列数少ない行があるみたい\n")
if s == 0:
	errlog.write("無事にError無く終わりましたよ\n")
elif s > 0:
	errlog.write("合計"+str(s)+"のErrorが見つかりましたよ\n")


