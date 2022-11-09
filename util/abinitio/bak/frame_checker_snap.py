#coding:utf-8

import sys
import random
import os
import string

genome=open(sys.argv[1],"r") #snapの予測に用いたゲノム配列(データの名前はシンプルに)
faa=open(sys.argv[2],"r") #snapが出力するアミノ酸配列
gff=open(sys.argv[3],"r") #.evm.gffのファイル

#1.ゲノム配列から1データ1行のファイルを作成
g1_tag=str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))
g1=open("genome."+g1_tag+".g1","w")
n=0
for g in genome:
	g = g.strip()
	if len(g) > 0:
		if g[0] == ">":
			if n == 0:
				n+=1
				g1.write(g.strip(">")+"\t")
			else:
				n+=1
                                g1.write("\n"+g.strip(">")+"\t")
		else:
			g1.write(g)
g1.write("\n")
genome.close()
g1.close()
##########################################
#2.step1のファイルから辞書を作成
dict_genome={}
g2=open("genome."+g1_tag+".g1","r")
for e in g2:
	e = e.strip()
	ee = e.split("\t")
	dict_genome[ee[0]]=ee[1]
g2.close()
os.system("rm genome."+g1_tag+".g1")
#########################################
#3.遺伝子配列から1データ1行のファイルを作成
f1_tag=str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))
f1=open("gene."+f1_tag+".f1","w")
m=0
#box=[]
for f in faa:
        f = f.strip()
        if len(f) > 0:
                if f[0] == ">":
                        if m == 0:
                                m+=1
                                f1.write(f.strip(">")+"\t")
                        else:
                                m+=1
                                f1.write("\n"+f.strip(">")+"\t")
                else:
                        f1.write(f)
f1.write("\n")
faa.close()
f1.close()
##########################################
#4.step3のファイルから辞書を作成
dict_gene={}
f2=open("gene."+f1_tag+".f1","r")
for a in f2:
        a = a.strip()
        aa = a.split("\t")
	#box.append(aa[0])
        dict_gene[aa[0]]=aa[1]
f2.close()
os.system("rm gene."+f1_tag+".f1")
#########################################
#for t in box:
#	print(">"+t)
#	print(dict_gene[t])
########################################
#########################################
#5.gffを参照し、frameを辞書に入れる
tmp_tag=str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))+"_"+str(random.randrange(1000000))
tmp=open("tmp."+tmp_tag+".tmp","w")
dict_gs={}
box=[]
for s in gff:
	s = s.strip()
	if len(s) > 0:
		ss = s.split("\t")
		sca = ss[0]
		typ = ss[2]
		st = int(ss[3]) - 1
		en = int(ss[4])
		di = ss[6]
		if typ == "mRNA":
			tmp.write(s+"\n")
		elif typ == "CDS":
			tmp.write(s+"\n")
			iD = ss[8].split(";")[0].split("=cds.model.")[1]
			if di=="+":
				if iD in box:
					dict_gs[iD]+=dict_genome[sca][st:en]
				elif iD not in box:
					box.append(iD)
					dict_gs[iD]=dict_genome[sca][st:en]
			elif di=="-":
				if iD in box:
					pre=dict_genome[sca][st:en]
					rev=pre[::-1]
					tran= rev.translate(string.maketrans("ATGCatgc","TACGtacg"))
                                        dict_gs[iD]+=tran
                                elif iD not in box:
                                        box.append(iD)
                                        pre=dict_genome[sca][st:en]
                                        rev=pre[::-1]
                                        tran= rev.translate(string.maketrans("ATGCatgc","TACGtacg"))
                                        dict_gs[iD]=tran
tmp.close()
gff.close()
#########################################
#for t in box:
#       print(">"+t)
#       print(dict_gs[t])
########################################
#6.コドン表を作成
DNA2Protein = {
        'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
        'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
        'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
        'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

        'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R',
        'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R',
        'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R',
        'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R',

        'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S',
        'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S',
        'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R',
        'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R',

        'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
        'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
        'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
        'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'
}
#######################################
#7.frameを特定
dict_tran={}
for l in box:
	dict_tran[l]={}
	seq0=dict_gs[l][0:]
	seq1=dict_gs[l][1:]
	seq2=dict_gs[l][2:]
	p1=""
	for i in range(0, int(len(seq0) / 3)):
		codon = seq0[3 * i : 3 * i + 3]
		if codon in DNA2Protein:
			p1 += DNA2Protein[codon]
		else:
			p1 += 'X'
	p2=""
        for x in range(0, int(len(seq1) / 3)):
                codon1 = seq1[3 * x : 3 * x + 3]
                if codon1 in DNA2Protein:
                        p2 += DNA2Protein[codon1]
                else:
                        p2 += 'X'
	p3=""
        for y in range(0, int(len(seq2) / 3)):
                codon2 = seq2[3 * y : 3 * y + 3]
                if codon2 in DNA2Protein:
                        p3 += DNA2Protein[codon2]
                else:
                        p3 += 'X'
	if p1 == dict_gene[l]:
		dict_tran[l]= 0
	elif p2 == dict_gene[l]:
                dict_tran[l]= 1
	elif p3 == dict_gene[l]:
                dict_tran[l]= 2
	else:
		print("ERROR:マッチする読み枠がない遺伝子を見つけました！["+l+"]")
	#print(l+"\t"+str(dict_tran[l]))
#######################################
#8.特定したframeをgffの最初のCDSに書き込む
gff_2nd=open("tmp."+tmp_tag+".tmp","r") #.evm.gffのファイル
check_box=[]
trans_table={2:1,1:2,0:0}
amari=0
s_w = 0
for k in gff_2nd:
	k = k.strip()
	kk = k.split("\t")
	Tp=kk[2]
	St=int(kk[3])
	En=int(kk[4])
	if Tp == "mRNA" and s_w == 0:
		s_w += 1
		print(k)
	elif Tp == "mRNA" and s_w > 0:
                print("\n"+k)
	elif Tp == "CDS":
		ID=kk[8].split(";")[0].split("=cds.model.")[1]
		if ID not in check_box:
			check_box.append(ID)
			kk[7]=str(dict_tran[ID])
			TSt=St+int(dict_tran[ID])
			size=En-TSt+1
			amari=size%3
			print("\t".join(kk))
		elif ID in check_box:
			slid=trans_table[amari]
			kk[7]=str(slid)
			TSt=St+slid
			size=En-TSt+1
			amari=size%3
			print("\t".join(kk))
gff_2nd.close()
os.system("rm tmp."+tmp_tag+".tmp")
####################################
