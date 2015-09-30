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
ndamp=700


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

./fmm_forward



paste data.txt traveltime.txt | awk '{print $6-$1}'
echo "Average residual!" 
paste data.txt traveltime.txt | awk '{a+=sqrt(($6-$1)*($6-$1))}END{print a/NF}'

echo "fd.txt" > inverseprb.inp
echo "data.txt" >> inverseprb.inp
echo "outmodel.txt" >> inverseprb.inp
echo "velmodel.txt" >> inverseprb.inp
echo "vel.txt" >> inverseprb.inp
echo "$ndata" >> inverseprb.inp
echo "$minlon $maxlon" >> inverseprb.inp
echo "$minlan $maxlan" >> inverseprb.inp
echo "$vconst" >> inverseprb.inp
echo "$ndamp" >> inverseprb.inp
echo "1 1" >> inverseprb.inp

./inverseprb

cp velresult.txt velresult0.txt
#cp velmodel.txt velmodel0.txt
for((i=1;i<4;i++))
do
	echo "source.txt" > fmm_forward.inp
	echo "receiver.txt" >> fmm_forward.inp
	echo "velresult.txt" >> fmm_forward.inp
	echo "tri.txt" >> fmm_forward.inp
	echo "$srcs" >> fmm_forward.inp
	echo "$dlon $dlan" >> fmm_forward.inp
	echo "5.0 5.0" >> fmm_forward.inp
	echo "0.2 0.2" >> fmm_forward.inp
	echo "$minlon $maxlon" >> fmm_forward.inp
	echo "$minlan $maxlan" >> fmm_forward.inp
	echo "2" >> fmm_forward.inp
	echo "1" >> fmm_forward.inp
	echo "2" >> fmm_forward.inp
	echo "$vconst" >> fmm_forward.inp
	echo "$nvtx" >> fmm_forward.inp
	echo "$vtxfile" >> fmm_forward.inp
	echo "$vout" >> fmm_forward.inp

	./fmm_forward



	paste data.txt traveltime.txt | awk '{print $6-$1}'
	echo "Average residual!" 
	paste data.txt traveltime.txt | awk '{a+=sqrt(($6-$1)*($6-$1))}END{print a/NF}'

	echo "fd.txt" > inverseprb.inp
	echo "data.txt" >> inverseprb.inp
	echo "outmodel.txt" >> inverseprb.inp
	echo "velmodel.txt" >> inverseprb.inp
	echo "vel.txt" >> inverseprb.inp
	echo "$ndata" >> inverseprb.inp
	echo "$minlon $maxlon" >> inverseprb.inp
	echo "$minlan $maxlan" >> inverseprb.inp
	echo "$vconst" >> inverseprb.inp
	echo "$ndamp" >> inverseprb.inp
	echo "1 1" >> inverseprb.inp

	./inverseprb
	cp velresult.txt velresult"$i".txt
done
