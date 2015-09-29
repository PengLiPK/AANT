#!/bin/bash

mkdir lowsnr

for file in $(cat snr_lt_10.txt | awk '{print $1}')
do
	mv aft-cut/"$file" lowsnr/
done
