#!/bin/bash

##Obviously a very simple script to take the appropriate
##columns from lightcurve and stick them in the file
##to be read by Lomb-Scargle program.
echo -n "Enter the read file name: "
read file

echo -n "Enter the second column number: "
read column 

echo -n "Write file name: "
read out

gawk -v col=$column '{print $1,$col, $col*0.001}' $file > $out

echo "file "$out" written with arbitrary errors"
                
