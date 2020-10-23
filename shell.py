#!/usr/bin/env python

from pylab import *
import os

datafile=open('sh.dat','wr')

tStart=1.
tMax=86400.
tIncr=1000.
shMass=2.e+12
shGamma=2.1
shWidth=1.e+3

x=arange(tStart,tMax,tIncr)

for t in x:
    s1=str(t)
    s2=str(shMass)
    s3=str(shGamma)
    s4=str(shWidth)
    datafile.write(s1+"\t"+s2+"\t"+s3+"\t"+s4+"\n")



