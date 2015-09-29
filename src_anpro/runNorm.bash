#!/bin/bash

mkdir aft-norm

cd aft-prcs
for dir in SCdata20040[1-7]
do
echo $dir
cd $dir
cp ../../norm .
for inpf in $(ls *.sac)
do
	echo $inpf
	tm=$(echo $inpf | cut -c1-8)
	net=$(echo $inpf | awk -F "." '{print $2}')
	stn=$(echo $inpf | awk -F "." '{print $3}')
	otf=$(echo nm$tm.$net.$stn.BHZ.sac)
  	echo "$inpf outf 100 1" | ./norm
	mv outf ../../aft-norm/$otf
done
cd ..
done
cd ..


#echo "$inpf 10 0"|./rdsac
#echo "$outf 10 1"|./rdsac
