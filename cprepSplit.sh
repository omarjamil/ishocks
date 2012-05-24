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

gawk -v col=$column '{print $1, $col, $col*0.001}' lightcurve1.dat > lc1_1_"$column"_err.dat
echo "File lc1_1_"$column"_err.dat written with made up errors!"
gawk -v col=$column '{print $1, $col}' lightcurve.dat > lc1_1_"$column".dat
echo "File lc1_1_"$column".dat written with no error column"

gawk -v col=$column '{print $1, $col, $col*0.001}' lightcurve2.dat > lc2_1_"$column"_err.dat
echo "File lc2_1_"$column"_err.dat written with made up errors!"
gawk -v col=$column '{print $1, $col}' lightcurve.dat > lc2_1_"$column".dat
echo "File lc2_1_"$column".dat written with no error column"

gawk -v col=$column '{print $1, $col, $col*0.001}' lightcurve3.dat > lc3_1_"$column"_err.dat
echo "File lc3_1_"$column"_err.dat written with made up errors!"
gawk -v col=$column '{print $1, $col}' lightcurve3.dat > lc3_1_"$column".dat
echo "File lc3_1_"$column".dat written with no error column"

gawk -v col=$column '{print $1, $col, $col*0.001}' lightcurve4.dat > lc4_1_"$column"_err.dat
echo "File lc4_1_"$column"_err.dat written with made up errors!"
gawk -v col=$column '{print $1, $col}' lightcurve4.dat > lc4_1_"$column".dat
echo "File lc4_1_"$column".dat written with no error column"

gawk -v col=$column '{print $1, $col, $col*0.001}' lightcurve5.dat > lc5_1_"$column"_err.dat
echo "File lc5_1_"$column"_err.dat written with made up errors!"
gawk -v col=$column '{print $1, $col}' lightcurve5.dat > lc5_1_"$column".dat
echo "File lc5_1_"$column".dat written with no error column"

gawk -v col=$column '{print $1, $col, $col*0.001}' lightcurve6.dat > lc6_1_"$column"_err.dat
echo "File lc6_1_"$column"_err.dat written with made up errors!"
gawk -v col=$column '{print $1, $col}' lightcurve6.dat > lc6_1_"$column".dat
echo "File lc6_1_"$column".dat written with no error column"

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

lcTstep=$(wc -l lc1_1_"$column"_err.dat)
set -- $lcTstep
lcBin=$(echo "scale=4;$totalT/$1" | bc)
oldBin=$(echo "scale=6;$dt-0.000001" | bc)
echo "lc1_1_"$column"_err.dat
$lcBin $oldBin 
irlc1.dat" > rebin1.par

lcTstep=$(wc -l lc2_1_"$column"_err.dat)
set -- $lcTstep
lcBin=$(echo "scale=4;$totalT/$1" | bc)
oldBin=$(echo "scale=6;$dt-0.000001" | bc)
echo "lc2_1_"$column"_err.dat
$lcBin $oldBin 
irlc2.dat" > rebin2.par

lcTstep=$(wc -l lc3_1_"$column"_err.dat)
set -- $lcTstep
lcBin=$(echo "scale=4;$totalT/$1" | bc)
oldBin=$(echo "scale=6;$dt-0.000001" | bc)
echo "lc3_1_"$column"_err.dat
$lcBin $oldBin 
irlc3.dat" > rebin3.par

lcTstep=$(wc -l lc4_1_"$column"_err.dat)
set -- $lcTstep
lcBin=$(echo "scale=4;$totalT/$1" | bc)
oldBin=$(echo "scale=6;$dt-0.000001" | bc)
echo "lc4_1_"$column"_err.dat
$lcBin $oldBin 
irlc4.dat" > rebin4.par

lcTstep=$(wc -l lc5_1_"$column"_err.dat)
set -- $lcTstep
lcBin=$(echo "scale=4;$totalT/$1" | bc)
oldBin=$(echo "scale=6;$dt-0.000001" | bc)
echo "lc5_1_"$column"_err.dat
$lcBin $oldBin 
irlc5.dat" > rebin5.par

lcTstep=$(wc -l lc6_1_"$column"_err.dat)
set -- $lcTstep
lcBin=$(echo "scale=4;$totalT/$1" | bc)
oldBin=$(echo "scale=6;$dt-0.000001" | bc)
echo "lc6_1_"$column"_err.dat
$lcBin $oldBin 
irlc6.dat" > rebin6.par

mv m.dat g.dat gm.dat lc_1_"$column".dat lc_1_"$column"_err.dat lc1_1_"$column".dat lc1_1_"$column"_err.dat lc2_1_"$column".dat lc2_1_"$column"_err.dat lc3_1_"$column".dat lc3_1_"$column"_err.dat lc4_1_"$column".dat lc4_1_"$column"_err.dat lc5_1_"$column".dat lc5_1_"$column"_err.dat lc6_1_"$column".dat lc6_1_"$column"_err.dat  rebin* crossCorr/

cd crossCorr
rm irlc*

echo "Rebinning lc_1_"$column"_err.dat in crossCorr"
./rebin
wait
echo "irlc.dat written"

echo "Rebinning lc1_1_"$column"_err.dat in crossCorr"
mv rebin1.par rebin.par
./rebin
wait
echo "irlc1.dat written"

echo "Rebinning lc2_1_"$column"_err.dat in crossCorr"
mv rebin2.par rebin.par
./rebin
wait
echo "irlc2.dat written"

echo "Rebinning lc3_1_"$column"_err.dat in crossCorr"
mv rebin3.par rebin.par
./rebin
wait
echo "irlc3.dat written"

echo "Rebinning lc4_1_"$column"_err.dat in crossCorr"
mv rebin4.par rebin.par
./rebin
wait
echo "irlc4.dat written"

echo "Rebinning lc5_1_"$column"_err.dat in crossCorr"
mv rebin5.par rebin.par
./rebin
wait
echo "irlc5.dat written"

echo "Rebinning lc6_1_"$column"_err.dat in crossCorr"
mv rebin6.par rebin.par
./rebin
wait
echo "irlc6.dat written"

cd ..