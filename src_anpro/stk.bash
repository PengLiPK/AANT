#!/bin/bash

mkdir aft-stack
cd aft-cut
cp ../stksac .
for stapr in $(ls | awk -F. '{print $2}' | sort -u)
do
	ls *$stapr*.sac > stack.inp
	echo $stapr begin
	echo "stack.inp 1" | ./stksac
	mv sum.sac ../aft-stack/sum$stapr.BHZ.sac
	rm stack.inp
done
cd ..
