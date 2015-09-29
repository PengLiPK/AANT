#!/bin/bash

cd aft-cut

cp ../cacuns .
for file in *.sac
do
	echo $file
	echo 1 > cacuns.inp
	echo $file >> cacuns.inp
	./cacuns
	saclst DEPMAX DEPMIN f $file > signal.out
	paste signal.out noise.out > snsingle.out
	cat snsingle.out >> ../sn.txt
done

cd ..
