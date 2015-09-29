#!/bin/bash


for period in 20.000000 19.000000 18.000000 17.000000 16.000000 15.000000 14.000000 13.000000 12.000000 11.000000 10.000000 9.0000000 8.0000000 7.0000000 6.0000000 5.0000000 4.0000000
do
	for file in *.gvel
	do
		pname=$(echo $period | awk -F. '{print $1}')
		pair=$(echo $file | awk -F_ '{print $1}')
		vel=$(grep "^   $period" $file | awk '{print $2}')
		dist=$(grep "$pair" stadist.txt | awk '{print $6}')
		trvtime=$(echo " scale=5; $dist / $vel " | bc )
		echo $pair $trvtime $vel $dist >> gvel_"$pname".txt
	done
	echo $period is finished!
done
