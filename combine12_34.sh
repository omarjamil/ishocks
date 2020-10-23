#!/bin/bash

echo -n "File name with columns 1, 2: "
read file1
echo -n "File name with columns 3, 4: "
read file2

gawk '{print $1, $2}' $file1 > temp12.dat
gawk '{print $3, $4}' $file2 > temp34.dat

paste temp12.dat temp34.dat > temp1234.dat

mv temp1234.dat $file1

rm temp12.dat temp34.dat
