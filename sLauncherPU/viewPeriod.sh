#!/bin/bash

echo -n "Enter the shell file name: "
read shFile

awk '{print $1,$2}' $shFile >> tempM.dat
awk '{print $1,$3}' $shFile >> tempG.dat

mv tempM.dat tempG.dat ../lombPeriodFFT/

cd ../lombPeriodFFT/
rm parameters.par
echo "input_file = tempM.dat
output_file = mScargle.dat
analysis_choice = 1
oversampling_factor = 4
frequency_factor = 1." > parameters.par
echo "Lomb-Scargle on Mass"
./periodogram
wait
echo " "
rm parameters.par
echo "input_file = tempG.dat
output_file = gScargle.dat
analysis_choice = 1
oversampling_factor = 4
frequency_factor = 1." > parameters.par                
echo "Lomb-Scargle on Gamma"
./periodogram
wait

gnuplot -persist "plot1.gp"
gnuplot -persist "plot2.gp"


