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
thrshd0=0.4
thrshd1=4
thrshd2=10
vout=1
nvtx=5
vtxfile=studyarea_vtx.txt
period=$(pwd | awk -F_ '{print $3}')

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

./fmm_adpt

cat edge_delgrid.txt newnode.txt > tempvel
wc -l tempvel | awk '{print $1}' > vel.txt
wc -l edge_delgrid.txt | awk '{print $1}' >> vel.txt
cat tempvel >>vel.txt
cp vel.txt vel_1.txt

cp G_1st_norm.txt G_1stnorm_1.txt
cp edge_delgrid.txt edge_delgrid_1.txt
cp fd.txt fd_1.txt
cp newnode.txt newnode_1.txt
cp delgrid.txt delgrid_1.txt

for((i=1;i<4;i++))
do
	j=$(( $i + 1 ))
	echo 2 > node.txt
	nodenum=$(head -1 vel.txt)
	echo $nodenum >>node.txt
	tail -$nodenum vel.txt | awk '{print $1,$2}' >> node.txt

	qdelaunay Qt i < node.txt >tri.txt


	echo "$period" > fmm_adpt.inp
	echo "source.txt" >> fmm_adpt.inp
	echo "receiver.txt" >> fmm_adpt.inp
	echo "vel.txt" >> fmm_adpt.inp
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

	./fmm_adpt

	
	cat edge_delgrid.txt newnode.txt > tempvel
	wc -l tempvel | awk '{print $1}' > vel.txt
	wc -l edge_delgrid.txt | awk '{print $1}' >> vel.txt
	cat tempvel >>vel.txt
	cp vel.txt vel_$j.txt

	cp G_1st_norm.txt G_1stnorm_$j.txt
	cp edge_delgrid.txt edge_delgrid_$j.txt
	cp fd.txt fd_$j.txt
	cp newnode.txt newnode_$j.txt
	cp delgrid.txt delgrid_$j.txt
done

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
echo "2" >> fmm_adpt_del.inp
echo "2" >> fmm_adpt_del.inp
echo "$vconst" >> fmm_adpt_del.inp
echo "$nvtx" >> fmm_adpt_del.inp
echo "$vtxfile" >> fmm_adpt_del.inp
echo "$vout" >> fmm_adpt_del.inp

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
