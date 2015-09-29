#!/bin/bash

mkdir part2_aftprcs
cd Rawdata

for file in *.SAC
do

sta=$(echo $file | awk -F. '{print $8}')
respf=$(echo TA_"$sta"_BHZ.resp)
echo $respf
echo $file

dlt=$(saclst DELTA f $file | awk '{print $2}')
if [ "$dlt" = 0.025 ]
then
sac <<EOF
r $file
trans from polezero s /nethome/pli/resp/$respf to none freq 0.005 0.006 39 40
dec 5
dec 4
dec 2
rmean
rtrend
w ../part2_aftprcs/$file
quit
EOF
else
echo $file >> prb.log
fi
echo $file
done
cd ..
