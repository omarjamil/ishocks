#!/bin/bash

echo -n "Enter the shellFile name: "
read shF
echo -n "Enter the total run time: "
read totalT
echo -n "Enter the time gap dt: "
read dt
echo -n "Enter the lc column number: "
read column

gawk -v col=$column '{print $1, $col, $col*0.001}' lightcurve.dat > lc_1_"$column"_err.dat

echo "File lc_1_"$column"_err.dat written with made up errors!"

gawk -v col=$column '{print $1, $col}' lightcurve.dat > lc_1_"$column".dat

echo "File lc_1_"$column".dat written with no error column"

gawk '{print $1, $2}' $shF > m.dat
gawk '{print $1, $3}' $shF > g.dat
gawk '{print $1, $2*$3}' $shF > gm.dat
echo "m.dat. g.dat and gm.dat written"

lcTstep=$(wc -l lc_1_"$column"_err.dat)
set -- $lcTstep
lcBin=$(echo "scale=4;$totalT/$1" | bc)
oldBin=$(echo "scale=6;$dt-0.000001" | bc)
echo "lc_1_"$column"_err.dat
$lcBin $oldBin 
irlc.dat" > rebin.par

mv m.dat g.dat gm.dat rebin.par lc_1_"$column".dat lc_1_"$column"_err.dat crossCorr/
echo "Rebinning lc_1_"$column"_err.dat in crossCorr"
cd crossCorr
rm irlc.dat
./rebin

wait
echo "Rebinning in crossCorr complete: irlc.dat writte"
cd ..