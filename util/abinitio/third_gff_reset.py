#coding:utf-8

import sys
import os

file = open(sys.argv[1],"r")

for i in file:
	i = i.strip()
	ii = i.replace("KIREME","\n")
	print(ii)
