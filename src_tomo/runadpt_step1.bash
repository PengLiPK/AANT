#/bin/bash


srcs=$(wc -l source.txt | awk '{print $1}')
data=$(wc -l data.txt | awk '{print $1}')
minlon=240.5
maxlon=242.6
minlan=34.3
maxlan=35.9
dlon=0.02
dlan=0.02
vconst=$(awk '{total+=$7}END{print total/'"$data"'}' data_op_coor.txt)
thrshd0=0.25
thrshd1=4
thrshd2=10
vout=1
nvtx=5
vtxfile=studyarea_vtx.txt
period=5

echo 2 > node.txt
nodenum=$(head -1 veltemp.txt)
echo $nodenum >>node.txt
tail -$nodenum veltemp.txt | awk '{print $1,$2}' >> node.txt

qdelaunay Qt i < node.txt
qdelaunay Qt i < node.txt > tri.txt

echo "$period" > fmm_adpt.inp
echo "source.txt" >> fmm_adpt.inp
echo "receiver.txt" >> fmm_adpt.inp
echo "veltemp.txt" >> fmm_adpt.inp
echo "tri.txt" >> fmm_adpt.inp
echo "$srcs" >> fmm_adpt.inp
echo "$data" >> fmm_adpt.inp
echo "$dlon $dlan" >> fmm_adpt.inp
echo "5.0 5.0" >> fmm_adpt.inp
echo "0.2 0.2" >> fmm_adpt.inp
echo "$minlon $maxlon" >> fmm_adpt.inp
echo "$minlan $maxlan" >> fmm_adpt.inp
echo "$thrshd0 $thrshd1 $thrshd2" >> fmm_adpt.inp
echo "2" >> fmm_adpt.inp
echo "2" >> fmm_adpt.inp
echo "2" >> fmm_adpt.inp
echo "$vconst" >> fmm_adpt.inp
echo "$nvtx" >> fmm_adpt.inp
echo "$vtxfile" >> fmm_adpt.inp
echo "$vout" >> fmm_adpt.inp
