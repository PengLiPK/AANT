#!/bin/bash

for crsfile in $(ls sum*.sac)
do
	stapr=$(echo "$crsfile" | awk -F. '{print $1}' | cut -c4- )
	side=$(echo "$crsfile" | awk -F. '{print $3}')
	echo "$stapr"
	grep "$stapr" stadist.txt
	dist=$(grep "$stapr" stadist.txt | awk '{print $6}')
	
	echo "$dist" > InputPrm.txt
	echo "0" >> InputPrm.txt
	echo "1" >> InputPrm.txt
	echo "3 80 1" >> InputPrm.txt
	echo "0.005" >> InputPrm.txt
	echo "1.5 6" >> InputPrm.txt
	echo "$crsfile" >> InputPrm.txt
	echo "$stapr"_"$side".out >> InputPrm.txt
	echo "$stapr"_"$side".gvel >> InputPrm.txt
	./ftan

done
