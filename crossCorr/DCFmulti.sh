#!/bin/bash
echo -n "Enter tau_min: "
read tau
echo -n "Enter the Delta_tau: "
read delta

echo "First_lc: g.dat 
Second_lc: irlc.dat
Outfile: corrg.dat
Tau_min: "$tau"
Delta_tau: "$delta"
"> read.par
echo "writing corrg.dat" 
./runDCF
wait
echo "corrg.dat done" 

echo "First_lc: m.dat 
Second_lc: irlc.dat
Outfile: corrm.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrm.dat" 
./runDCF
wait
echo "corrm.dat done" 

echo "First_lc: g.dat 
Second_lc: irlc1.dat
Outfile: corrg1.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrg1.dat" 
./runDCF
wait
echo "corrg1.dat done" 

echo "First_lc: m.dat 
Second_lc: irlc1.dat
Outfile: corrm1.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrm1.dat" 
./runDCF
wait
echo "corrm1.dat done" 

echo "First_lc: g.dat 
Second_lc: irlc2.dat
Outfile: corrg2.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrg2.dat" 
./runDCF
wait
echo "corrg2.dat done" 

echo "First_lc: m.dat 
Second_lc: irlc2.dat
Outfile: corrm2.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrm2.dat" 
./runDCF
wait
echo "corrm2.dat done" 

echo "First_lc: g.dat 
Second_lc: irlc3.dat
Outfile: corrg3.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrg3.dat" 
./runDCF
wait
echo "corrg3.dat done" 

echo "First_lc: m.dat 
Second_lc: irlc3.dat
Outfile: corrm3.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrm3.dat" 
./runDCF
wait
echo "corrm3.dat done" 

echo "First_lc: g.dat 
Second_lc: irlc4.dat
Outfile: corrg4.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
./runDCF
wait
echo "corrg4.dat done" 

echo "First_lc: m.dat 
Second_lc: irlc4.dat
Outfile: corrm4.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrm4.dat" 
./runDCF
wait
echo "corrm4.dat done" 

echo "First_lc: g.dat 
Second_lc: irlc5.dat
Outfile: corrg5.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrg5.dat" 
./runDCF
wait
echo "corrg5.dat done" 

echo "First_lc: m.dat 
Second_lc: irlc5.dat
Outfile: corrm5.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrm5.dat" 
./runDCF
wait
echo "corrm5.dat done" 

echo "First_lc: g.dat 
Second_lc: irlc6.dat
Outfile: corrg6.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrg6.dat" 
./runDCF
wait
echo "corrg6.dat done" 

echo "First_lc: m.dat 
Second_lc: irlc6.dat
Outfile: corrm6.dat
Tau_min: "$tau"
Delta_tau: "$delta"
" > read.par
echo "writing corrm6.dat" 
./runDCF
wait
echo "corrm6.dat done" 








