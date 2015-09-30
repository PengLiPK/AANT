#!/bin/bash


srcs=$(wc -l source.txt | awk '{print $1}')
data=$(wc -l data.txt | awk '{print $1}')
minlon=80.5
maxlon=82.6
minlan=34.3
maxlan=35.9
dlon=0.02
dlan=0.02
vconst=2.9

# Delete node with small ray weight
echo 2 > node.txt
nodenum=$(head -1 vel.txt)
echo $nodenum >>node.txt
tail -$nodenum vel.txt | awk '{print $1,$2}' >> node.txt

qdelaunay Qt i < node.txt >tri.txt


echo "source.txt" > fmm_adpt_del.inp
echo "receiver.txt" >> fmm_adpt_del.inp
echo "vel.txt" >> fmm_adpt_del.inp
echo "tri.txt" >> fmm_adpt_del.inp
echo "$srcs" >> fmm_adpt_del.inp
echo "$data" >> fmm_adpt_del.inp
echo "$dlon $dlan" >> fmm_adpt_del.inp
echo "5.0 5.0" >> fmm_adpt_del.inp
echo "0.2 0.2" >> fmm_adpt_del.inp
echo "$minlon $maxlon" >> fmm_adpt_del.inp
echo "$minlan $maxlan" >> fmm_adpt_del.inp
echo "2" >> fmm_adpt_del.inp
echo "1" >> fmm_adpt_del.inp
echo "2" >> fmm_adpt_del.inp
echo "$vconst" >> fmm_adpt_del.inp

./fmm_adpt_del


cat edge_delgrid.txt newnode.txt > tempvel
wc -l tempvel | awk '{print $1}' > vel.txt
wc -l edge_delgrid.txt | awk '{print $1}' >> vel.txt
cat tempvel >>vel.txt
cp vel.txt vel_5.txt

cp G_1st_norm.txt G_1stnorm_5.txt
cp edge_delgrid.txt edge_delgrid_5.txt
cp fd.txt fd_5.txt
cp newnode.txt newnode_5.txt
cp delgrid.txt delgrid_5.txt
