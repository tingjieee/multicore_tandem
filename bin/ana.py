#!/usr/bin/env python
# -*- coding: utf-8 -*-

f = open('Protein_cnt.csv','r');
lines = f.readlines();

s1 = set();

for line in lines:
	line = line[:-1];
	elems = tuple(line.split('|'));
#	print elems[0];
	s1.add(elems[0]);
print len(s1);
#s1.add(1);
f1 = open('../../tandem-bak-serial/bin/protein_cnt.csv','r');
#f1 = open('../../PTandem-pipline/bin/Protein_cnt.csv','r');

lines = f1.readlines();

s2 = set();

for line in lines:
	line = line[:-1];
	elems = tuple(line.split('|'));
#	print elems[0];
	s2.add(elems[0]);
print len(s2);

s3 = s1 & s2;
#for i in s3:
#	print i;
#	j = j+1;
print len(s3);
