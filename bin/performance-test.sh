#!/bin/bash
i=0;
for i in 1 
do
#	time ./ptandem input.xml
	for t in {2..8..2}
	do
		echo $t;
		./ptandem input.xml
	done
done
