#!/bin/bash

echo "First_lc: g.dat 
Second_lc: irlc.dat
Outfile: corrg.dat
Tau_min: 10
Delta_tau: 0.1
" > read.par
echo "writing corrg.dat" 
./runDCF
wait
echo "corrg.dat done" 

echo "First_lc: m.dat 
Second_lc: irlc.dat
Outfile: corrm.dat
Tau_min: 10
Delta_tau: 0.1
" > read.par
echo "writing corrm.dat" 
./runDCF
wait
echo "corrm.dat done" 

echo "First_lc: gm.dat 
Second_lc: irlc.dat
Outfile: corrgm.dat
Tau_min: 10
Delta_tau: 0.1
" > read.par
echo "writing corrgm.dat" 
./runDCF
wait
echo "corrgm.dat done" 