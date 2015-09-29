#!/bin/bash



for((pd=4;pd<21;pd++))
do
	paste stadist.txt gvel_$pd.txt | awk '{print $1,$2,$3,$4,$5,$9,$10,$11,$8}' > data_op_coor_$pd.txt
done
