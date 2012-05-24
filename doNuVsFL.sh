#!/bin/bash

#echo -n "Enter the lightcurve line number > "
#read lineNumber

#gawk -v var=$lineNumber 'NR==var {for (i = 2; i <= NF; i++) print $i}' 
#lightcurve.dat > temp1.dat

tail -n 1 lightcurve.dat | awk '{for (i = 2; i <= NF; i++) print $i}' > temp1.dat
#awk -f 'transpose.awk' temp1.dat > temp2.dat

paste nuRange.dat temp1.dat > nuVsFL.dat
rm temp1.dat
