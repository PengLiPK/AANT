#!/bin/bash

file=data_op_coor.txt

awk '{print $6,$4,$5,$2,$3}' $file > data.txt

awk '{print $2+360,$3}' $file | uniq > source.txt

i=0
rm receiver.txt
for sta in $(awk '{print substr($1,1,4)}' $file | uniq)
do
	i=$(( $i + 1 ))
	echo rfile"$i".txt >> receiver.txt
	grep "$sta" $file | wc -l > rfile"$i".txt
	grep "$sta" $file | awk '{print $4+360,$5}' >> rfile"$i".txt
done
