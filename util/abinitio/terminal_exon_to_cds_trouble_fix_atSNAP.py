#coding:utf-8
#SNAPはexonを予測するため、CDSとしてそのまま翻訳すると"*X"が見られることがある
#このスクリプトはSNAPのgffをCDSベースに変えるプログラムである。
#入力は塩基配列、アミノ酸配列、gff配列である

import sys
import os

#塩基配列
fna=open(sys.argv[1],"r")
#アミノ酸配列(SNAPが出力したアミノ酸配列)
faa=open(sys.argv[2],"r")
#gff
gff=open(sys.argv[3],"r")

#SANP予測のアミノ酸の配列と長さを記録
faa_dict={}
data1=""
for a in faa:
	a = a.strip()
	if len(a) > 0:
		if a[0] == ">":
			data1=a.strip(">")
			faa_dict[data1]=["",0]
		else:
			faa_dict[data1][0]=faa_dict[data1][0]+a
			faa_dict[data1][1]=faa_dict[data1][1]+len(a)
faa.close()
#塩基配列の配列と長さを記録
fna_dict={}
data2=""
for n in fna:
	n=n.strip()
	if len(n) > 0:
		if n[0] == ">":
			data2=n.strip(">")
			fna_dict[data2]=["",0]
		else:
			fna_dict[data2][0]=fna_dict[data2][0]+n
			fna_dict[data2][1]=fna_dict[data2][1]+len(n)
fna.close()
#アミノ酸配列と塩基配列を比較
length_fna=len(fna_dict)
length_faa=len(faa_dict)
rm_dict={}
if length_fna != length_faa:
	print("The number of nucleotide sequences and amino acid sequences is different !!")
	sys.exit()
for i in faa_dict:
	nucl_len=fna_dict[i][1]
	amin_len=faa_dict[i][1]*3
	if nucl_len != amin_len:
		dif=nucl_len-amin_len
		if "*" in faa_dict[i][0]:
			rm_dict[i]=dif
		#print(i+"\t"+str(dif))
		#print("#nucl:"+fna_dict[i][0])
		#print("#amin:"+faa_dict[i][0])
#GFFの処理
n=0
tmp=open(sys.argv[3]+".tmp","w")
for g in gff:
	g=g.strip()
	g_type=g.split("\t")[2]
	if g_type == "mRNA" and n==0:
		n=1
		tmp.write(g)
	elif g_type == "mRNA" and n!=0:
		tmp.write("\n"+g)
	else:
		tmp.write("@@_@@"+g)
tmp.write("\n")
gff.close()
tmp.close()
tmp2=open(sys.argv[3]+".tmp","r")
for t in tmp2:
	t=t.strip()
	tt=t.split("@@_@@")
	word=tt[0].split("\t")[8].split("ID=")[1]
	if word in rm_dict:
		strand=tt[0].split("\t")[6]
		if strand=="+":
			#mrnaの値を変更
			mrna=tt[0].split("\t")
			num=int(mrna[4]) - rm_dict[word]
			mrna[4]=str(num)
			tt[0]="\t".join(mrna)
			#CDSの値変更
			cds=tt[-1].split("\t")
			num=int(cds[4]) - rm_dict[word]
			cds[4]=str(num)
			tt[-1]="\t".join(cds)
		elif strand=="-":
			#mrnaの値を変更
			mrna=tt[0].split("\t")
			num=int(mrna[3]) + rm_dict[word]
			mrna[3]=str(num)
			tt[0]="\t".join(mrna)
			#CDSの値を変更
			cds=tt[-1].split("\t")
			num=int(cds[3]) + rm_dict[word]
			cds[3]=str(num)
			tt[-1]="\t".join(cds)
		sys.stdout.write("\n".join(tt)+"\n")
	else:
		sys.stdout.write("\n".join(tt)+"\n")
tmp2.close()
os.system("rm "+sys.argv[3]+".tmp")
