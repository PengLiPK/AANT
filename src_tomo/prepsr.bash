#!/bin/bash

sourcenum=30
outf=testsource.txt



echo $outf > prepsr.inp
echo $sourcenum >> prepsr.inp
echo "100 110" >> prepsr.inp
echo "30 40" >> prepsr.inp

./prepsr

echo "Random points generating finished!"

num1=$(wc -l $outf | awk '{print $1}')
num2=$(($num1 - 1))
echo $num2
tail -$num2 $outf | psxy -JM4i -R100/110/30/40 -Ba5f2g1 -Sa0.3 -Gred > testsource.ps

echo "Drawing figure finished!"
