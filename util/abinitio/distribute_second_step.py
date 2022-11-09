#coding:utf-8

import os
import sys

Num=int(sys.argv[2])

#ファイル生成
rang=range(Num)
for i in rang:
	nUm=i+1
	TF=open("sPlitfile"+str(nUm)+".tMp","w")
	TF.close()
	del TF

del i

#出力ファイルに分配
m=1
x=Num
lIst=open(sys.argv[1],"r")
for i in lIst:
	i = i.strip()
	ii = i.split("@")[1:]
	iii="".join(ii)	
	if m <= Num:
		TF=open("sPlitfile"+str(m)+".tMp","a")
		TF.write(iii+"\n")
		m+=1
		TF.close()
	elif m > Num:
		if x > 0:
			TF=open("sPlitfile"+str(x)+".tMp","a")
			TF.write(iii+"\n")
			x=x-1
			TF.close()
		elif x == 0:
			m = 1
			x = Num
			TF=open("sPlitfile"+str(m)+".tMp","a")
			TF.write(iii+"\n")
			TF.close()
			m+=1
		#m=1
		#TF=open("sPlitfile"+str(m)+".tMp","a")
		#TF.write(iii+"\n")
	del TF
lIst.close()

##fasta生成
#NUM=range(Num)

#for i in NUM:
#	ii=str(i+1)
#	os.system("python /data/yuasa/script/augustus_crafty/tab_separate.py sPlit"+str(ii)+".tMp > sPlit"+str(ii)+".fasta")
#
#os.system("rm *.tMp")
