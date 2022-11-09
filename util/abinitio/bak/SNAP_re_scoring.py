#coding:utf-8

import sys

#############
snap=open(sys.argv[1],"r")
score_box=[]
for i in snap:
	ii = i.split("\t")
	if ii[5]!=".":
		score_box.append(float(ii[5]))
snap.close()
sorted_box=sorted(set(score_box), reverse=True)
length=len(sorted_box)
unit=length//10
#############
re_snap=open(sys.argv[1],"r")
for l in re_snap:
	l=l.strip()
	ll = l.split("\t")
	if ll[5]!=".":
		S=float(ll[5])
		if S>=sorted_box[unit]:
			ll[5]="9"
			print("\t".join(ll))
		elif sorted_box[unit]>S>=sorted_box[unit*2]:
			ll[5]="8"
			print("\t".join(ll))
		elif sorted_box[unit*2]>S>=sorted_box[unit*3]:
			ll[5]="7"
			print("\t".join(ll))
		elif sorted_box[unit*3]>S>=sorted_box[unit*4]:
			ll[5]="6"
			print("\t".join(ll))
		elif sorted_box[unit*4]>S>=sorted_box[unit*5]:
			ll[5]="5"
			print("\t".join(ll))
		elif sorted_box[unit*5]>S>=sorted_box[unit*6]:
			ll[5]="4"
			print("\t".join(ll))
		elif sorted_box[unit*6]>S>=sorted_box[unit*7]:
			ll[5]="3"
			print("\t".join(ll))
		elif sorted_box[unit*7]>S>=sorted_box[unit*8]:
			ll[5]="2"
			print("\t".join(ll))
		elif sorted_box[unit*8]>S>=sorted_box[unit*9]:
			ll[5]="1"
			print("\t".join(ll))
		elif sorted_box[unit*9]>S:
			ll[5]="0"
			print("\t".join(ll))
	elif ll[5]==".":
		print(l)
re_snap.close()

