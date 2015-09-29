#!/bin/bash


mkdir aft-norm_cvt
cd aft-norm


for file in *.sac
do
	saccvt -I < $file > temp.sac
	mv temp.sac ../aft-norm_cvt/$file
done

