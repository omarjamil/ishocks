#!/bin/bash

if [ -e "lightcurve.dat" ];then
    rm lightcurve*
    echo "deleted: all lightcurve*"
fi

if [ -e "resultsFile.dat" ];then
    rm resultsFile.dat
    echo "deleted: resultsFile.dat"
fi

if [ -e "tauResFile.dat" ];then
    rm tauResFile.dat
    echo "deleted: tauResFile.dat"
fi

if [ -e "finalTS.dat" ];then
    rm finalTS.dat
    echo "deleted: finalTS.dat"
fi

if [ -e "finalTStau.dat" ];then
    rm finalTStau.dat
    echo "deleted: finalTStau.dat"
fi

if [ -e "nuRange.dat" ];then
    rm nuRange.dat
    echo "deleted: nuRange.dat"
fi

