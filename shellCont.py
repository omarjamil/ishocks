#!/usr/bin/env python

from pylab import *
import os

datafile=open('shCont.dat','wr')

tStart=1.
tMax=10000.
tIncr=1.
jetLum=1.e31
lScale=0.6
c=3.e8
shGMin=2.9
shGMax=3.0
shGRange=(shGMax-shGMin)
#shMass=2.e+12
dtInjMin=1.
dtInjMax=1000.
dtInjRange=tInjMax-tInjMin

x=arange(tStart,tMax,tIncr)

for t in x:
    sinu=sin(0.0023*t)**2
    #shGamma=2.1
    shGamma=(sinu*shGRange)+shGMin
    dtInj=(sin*dtInjRange)+dtInjMin
    shBeta=sqrt(1.-1./(shGamma**2))
    shLength=lScale*tIncr*shBeta*c
    shMass=(jetLum*tIncr)/(shGamma*(c**2))
    s1=str(t)
    s2=str(shMass)
    s3=str(shGamma)
    s4=str(shLength)
    datafile.write(s1+"\t"+s2+"\t"+s3+"\t"+s4+"\n")

