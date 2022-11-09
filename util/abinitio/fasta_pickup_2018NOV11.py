#coding:utf-8
import sys

infile = open(sys.argv[1],'r')
name = sys.argv[2] #チェックする文字
s = 0

rmk = open('rmk.fasta','r')
for o in rmk:
	o = o.strip()
	oo = o.strip(">")
	if o[0] == '>' and name == oo:
		print(o)
		s = 1
	elif '>' not in o and s == 1:
		print(o)
		s = 0
rmk.close()
