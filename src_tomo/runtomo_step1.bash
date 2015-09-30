#!/bin/bash


srcs=$(wc -l source.txt | awk '{print $1}')
ndata=$(wc -l data.txt | awk '{print $1}')
minlon=240.5
maxlon=242.6
minlan=34.3
maxlan=35.9
dlon=0.02
dlan=0.02
vconst=$(awk '{total+=$7}END{print total/'"$ndata"'}' data_op_coor.txt)
vout=1
nvtx=5
vtxfile=studyarea_vtx.txt
period=5


#cat edge_delgrid.txt newnode.txt > tempvel
#wc -l tempvel | awk '{print $1}' > vel.txt
#wc -l edge_delgrid.txt | awk '{print $1}' >> vel.txt
#cat tempvel >>vel.txt

echo 2 > node.txt
nodenum=$(head -1 vel.txt)
echo $nodenum >>node.txt
tail -$nodenum vel.txt | awk '{print $1,$2}' >> node.txt

qdelaunay Qt i < node.txt >tri.txt


echo "source.txt" > fmm_forward.inp
echo "receiver.txt" >> fmm_forward.inp
echo "vel.txt" >> fmm_forward.inp
echo "tri.txt" >> fmm_forward.inp
echo "$srcs" >> fmm_forward.inp
echo "$dlon $dlan" >> fmm_forward.inp
echo "5.0 5.0" >> fmm_forward.inp
echo "0.2 0.2" >> fmm_forward.inp
echo "$minlon $maxlon" >> fmm_forward.inp
echo "$minlan $maxlan" >> fmm_forward.inp
echo "2" >> fmm_forward.inp
echo "2" >> fmm_forward.inp
echo "2" >> fmm_forward.inp
echo "$vconst" >> fmm_forward.inp
echo "$nvtx" >> fmm_forward.inp
echo "$vtxfile" >> fmm_forward.inp
echo "$vout" >> fmm_forward.inp
