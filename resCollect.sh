#!/bin/bash

folder=$1
echo $folder

mkdir $folder

cp parameters.par splitLC.par $folder
cp *.dat $folder
cp sLauncherPU/launch.in $folder
cp crossCorr/*.dat $folder

mv jets.log $folder
