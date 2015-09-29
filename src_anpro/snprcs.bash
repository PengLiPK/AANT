#!/bin/bash


awk '{print $1,$2/$5,$2/$6}' sn.txt > snr.txt

awk '{if($2<10) print $1,$2,$3}' snr.txt > snr_lt_10.txt
