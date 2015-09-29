#!/bin/bash

sta1=ALP
sta2=CLC

mkdir aft-whiten

cd aft-norm_cvt

for file in *.sac
do

gsac << EOF
r $file
whiten freqlimits 0.001 0.005 0.13 0.3 absolute
w ../aft-whiten/wh$file
quit
EOF
done
cd ..

