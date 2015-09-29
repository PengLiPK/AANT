#!/bin/bash

###########################################################
# Rename
###########################################################

mkdir aftprcs-rn
cd part2_aftprcs

for file in *.SAC
do
	yr=$(echo $file | awk -F. '{print $1}')
	day=$(echo $file | awk -F. '{print $2}')
	sta=$(echo $file | awk -F. '{print $8}')
	cp $file ../aftprcs-rn/$yr$day.TA.$sta.BHZ.R.sac
done
cd ..

cd aftprcs-rn
mkdir badfile
ls -l | awk '{if($5<346232)print $0}' >> ../prbdata.log
bdfile=$(ls -l | awk '{if($5<346232)print $NF}' | grep "sac")
mv $bdfile badfile/

cd ..

###########################################################
# Norm
###########################################################
mkdir aft-norm

cd aftprcs-rn
cp ../norm .
for file in *.sac
do
	echo "$file outf 100 1" | ./norm
	mv outf ../aft-norm/nm"$file"
done
cd ..


###########################################################
# Whitening
###########################################################
mkdir aft-whiten

cd aft-norm

for file in *.sac
do

gsac << EOF
r $file
whiten absolute
w ../aft-whiten/wh$file
quit
EOF
done
cd ..


############################################################
# Cross-correlation
############################################################

mkdir aft-crscrltn

cd aft-whiten

for ftime in $( ls *.sac | awk -F. '{print $1}' | sort -u )
do

echo $ftime begin
stanum=$(ls $ftime*.sac | wc | awk '{print $1}')
fmaxn=$(echo "$stanum - 1" | bc)
sacf1day=($(ls $ftime*.sac))

for((i=0;i<=$fmaxn;i++))
do
echo i $i
sta1=$(echo ${sacf1day[i]} | awk -F. '{print $3}')
for((j=i+1;j<=$fmaxn;j++))
do
echo j $j
sta2=$(echo ${sacf1day[j]} | awk -F. '{print $3}')

sac <<EOF
r $ftime*$sta1*.sac $ftime*$sta2*.sac
cor m 1
w ../aft-crscrltn/$ftime.$sta1-$sta1.BHZ.sac ../aft-crscrltn/$ftime.$sta1-$sta2.BHZ.sac
quit
EOF

done
done

done
cd ..


############################################################################
# move sta-pairs with same station
############################################################################
mkdir aft-crs-smsta
cd aft-crscrltn
cp ../sta.lst .
for sta in $(cat sta.lst)
do
	echo $sta
	mv *."$sta"-"$sta".*.sac ../aft-crs-smsta/
done
cd ..


############################################################################
# cut
############################################################################
mkdir aft-cut
cd aft-crscrltn
cp ../cutsac .

for sacfile in *.sac
do
	echo 1 > cutsac.inp
	echo "$sacfile" >> cutsac.inp
	echo out.sac >> cutsac.inp
	echo -10000 >> cutsac.inp
	echo 10000 >> cutsac.inp
	echo "$sacfile"
	./cutsac
	mv out.sac ../aft-cut/ct"$sacfile"
done
cd ..

############################################################################
# Delete low SNR files
############################################################################
#cd aft-cut

#cp ../cacuns .
#for file in *.sac
#do
#	echo $file
#	echo 1 > cacuns.inp
#	echo $file >> cacuns.inp
#	./cacuns
#	saclst DEPMAX DEPMIN f $file > signal.out
#	paste signal.out noise.out > snsingle.out
#	cat snsingle.out >> ../sn.txt
#done
#cd ..

#awk '{print $1,($2-$3)/(2*$5),($2-$3)/(2*$6)}' sn.txt > snr.txt
#awk '{if($2<10)print $0}' snr.txt > snr_lt_10.txt
#mkdir lowsnr
#for file in $(cat snr_lt_10.txt | awk '{print $1}')
#do
#	mv aft-cut/"$file" lowsnr/
#done


############################################################################
# Stacking
############################################################################

mkdir aft-stack

cd aft-cut
cp ../stksac .
for stapr in $(ls | awk -F. '{print $2}' | sort -u)
do
	ls *$stapr*.sac > stack.inp
	echo $stapr begin
	echo "stack.inp 1" | ./stksac
	mv sum.sac ../aft-stack/sum$stapr.BHZ.sac
	rm stack.inp
done
cd ..
