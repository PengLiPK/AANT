#!/bin/bash

rm del.log
for dfile in $(ls data_op_coor_*.txt)
do
	echo $dfile >>del.log
	awk '{if($7>3.5 || $7<2.005) print $0}' $dfile >>del.log
	for stapr in $(awk '{if($7>3.5 || $7<2.005) print $1}' $dfile)
	do
		echo $stapr
		sed -i '/^'"$stapr"'/d' $dfile
	done
done

