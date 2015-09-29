#!/bin/bash
# Find small size data, move them to trash.

for dir in SCdata*
do
	echo $dir
	cd $dir
	mkdir badfile
	ls -l | awk '{if($5<6912636)print $0}' >> ../prbmdat.log
	bdfile=$( ls -l | awk '{if($5<6912636)print $NF}' | grep "sac")
	echo $bdfile
	mv $bdfile badfile
	cd ..
done
