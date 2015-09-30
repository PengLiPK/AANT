#!/bin/bash


rcs=30
srcs=10
dlon=0.1
dlan=0.1
vconst=4.0


cat edge_delgrid.txt newnode.txt > tempvel
wc -l tempvel | awk '{print $1}' > vel.txt
wc -l edge_delgrid.txt | awk '{print $1}' >> vel.txt
cat tempvel >>vel.txt

echo 2 > node.txt
nodenum=$(head -1 vel.txt)
echo $nodenum >>node.txt
tail -$nodenum vel.txt | awk '{print $1,$2}' >> node.txt

qdelaunay Qt i < node.txt >tri.txt


echo "source.txt" > fmm.inp
echo "receiver.txt" >> fmm.inp
echo "vel.txt" >> fmm.inp
echo "tri.txt" >> fmm.inp
echo "$srcs" >> fmm.inp
echo "$rcs" >> fmm.inp
echo "$dlon $dlan" >> fmm.inp
echo "5.0 5.0" >> fmm.inp
echo "1.0 1.0" >> fmm.inp
echo "100.00 110.0" >> fmm.inp
echo "30.00 40.00" >> fmm.inp
echo "2" >> fmm.inp
echo "1" >> fmm.inp
echo "2" >> fmm.inp
echo "4.0" >> fmm.inp

./fmm_forward



