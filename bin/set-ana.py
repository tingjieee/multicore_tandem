#!/usr/bin/env python
# -*- coding: utf-8 -*-

f1 = open('Hydrolyze1','r');
lines1 = f1.readlines();
l1 = [];
for line in lines1:
	line = line[:-1]
	elems = tuple(line.split(','))
	l1.append(elems);
#print(len(l1))

f2 = open('Hydrolyze0','r');
#f2 = open('../../tandem-bak-serial/bin/Hydrolyze','r');
lines2 = f2.readlines();
l2 = [];
for line in lines2:
	line = line[:-1]
	elems = tuple(line.split(','))
	l2.append(elems);
#print(len(l2))

#i  = 0
k = 0;
for item in l1:
	s1 = set()
	s2 = set()
	j = 0
	ss = ""
	for ele in item:
		if(j == 2) :
#			print ele
			s1.add(ele)
		else :
			j = j+1;
#	for p in s1 :
#		print p,
#	print "\n"

	i = 0
	for ele in l2[k] :
		if(i == 2) :
			s2.add(ele)
		else :
			i = i+1
	k = k + 1
	print s1 < s2

#print len(s1);
'''	
#k  = 0
s2 = set()
for item in l1:
	j = 0
	for ele in item:
		if(j == 2) :
			s2.add(ele)
		else :
			j = j+1;
#print len(s2);
'''

