#!/bin/bash

mkdir aft-crscrltn
cd aft-whiten

for ftime in $( ls *.sac | awk -F. '{print $1}' | sort -u )
do

echo $ftime begin
size=$(ls $ftime*.sac | wc | awk '{print $1}')
fmaxn=$(echo "$size - 1" | bc)
sacf1day=($(ls $ftime*.sac))

for((i=0;i<$fmaxn;i++))
do
echo i $i
sta1=$(echo ${sacf1day[i]} | awk -F. '{print $3}')
for((j=i+1;j<=$fmaxn;j++))
do
echo j $j
sta2=$(echo ${sacf1day[j]} | awk -F. '{print $3}')

sac <<EOF
r $ftime*$sta1*.sac $ftime*$sta2*.sac
cor m 1
w ../aft-crscrltn/$ftime.$sta1-$sta1.BHZ.sac ../aft-crscrltn/$ftime.$sta1-$sta2.BHZ.sac
quit
EOF

done
done

done
cd ..
