#!/bin/bash


paste data.txt traveltime.txt | awk '{print $6-$1}' > resabs.txt
paste data.txt traveltime.txt | awk '{print ($6-$1)/$1}' > resperc.txt

num=$(wc -l resperc.txt | awk '{print $1}')

aveabs=$(awk '{a+=$1}END{print a/"'$num'"}' resabs.txt)
aveperc=$(awk '{a+=$1}END{print a/"'$num'"}' resperc.txt)

stdperc=$(awk '{a+=($1-"'$aveperc'")*($1-"'$aveperc'")}END{print sqrt(a/"'$num'")}' resperc.txt)
stdabs=$(awk '{a+=($1-"'$aveabs'")*($1-"'$aveabs'")}END{print sqrt(a/"'$num'")}' resabs.txt)


abave=$(awk '{a+=sqrt($1*$1)}END{print a/"'$num'"}' resperc.txt)

echo $aveabs ave > resabs_stat.txt
echo $stdabs std >> resabs_stat.txt

echo $aveperc ave > resperc_stat.txt
echo $stdperc std >> resperc_stat.txt
echo $abave abave >> resperc_stat.txt

cat resperc_stat.txt
