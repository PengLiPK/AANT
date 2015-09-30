#/bin/bash


srcs=$(wc -l source.txt | awk '{print $1}')
data=$(wc -l data.txt | awk '{print $1}')
minlon=66.0
maxlon=67.8
minlan=34.3
maxlan=35.9
dlon=0.02
dlan=0.02
echo $data
vconst=$(awk '{total+=$7}END{print total/'"$data"'}' data_op_coor.txt)
echo $vconst
