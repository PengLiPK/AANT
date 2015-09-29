#!/bin/bash

#for sacfile in $(ls sum*.BHZ.sac)
#do
#sac <<EOF
#r $sacfile
#bp co 0.005 0.5
#w over
#quit
#EOF
#done

mkdir aft-cut-lr

cd aft-stack

cp ../cutsac .

for stapr in $(ls sum*.sac | awk -F. '{print $1}')
do
	echo 1 > cutsac.inp
	echo "$stapr".BHZ.sac >> cutsac.inp
	echo "$stapr".BHZ.sac
	./cutsac
	mv left.sac ../aft-cut-lr/"$stapr".BHZ.left.sac
	mv right.sac ../aft-cut-lr/"$stapr".BHZ.right.sac
done

cd ..
