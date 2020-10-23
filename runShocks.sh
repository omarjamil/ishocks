#!/bin/bash

./cleanResFile.sh
rm locCount.dat
./ishocks
wait
if [ -e "locCount.dat" ];then
    gawk -f sumall.awk locCount.dat    
fi

echo "All done. Have a nice day!"