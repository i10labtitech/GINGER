#coding:utf-8
import sys

infile = open(sys.argv[1],'r')
rmkfile = open('rmk.fasta','w')
switch = 0

for i in infile:
	i = i.strip()
	if len(i) > 0:
		if i[0] == '>' and switch == 0:
			rmkfile.write(i)
			rmkfile.write('\n')
			switch = 1
		elif i[0] == '>' and switch == 1:
			rmkfile.write('\n' + i)
			rmkfile.write('\n')
		else:
			rmkfile.write(i)
rmkfile.close()
