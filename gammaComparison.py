#!/usr/bin/env python

from pylab import *
import os

datafile=open('gamma.dat','wr')

gMin=1.5
gMax=2.0
incr=0.001
mass=1.0
eta=0.0
mu=mass+eta

g1=arange(gMin,gMax,incr)
g2=arange(gMin,gMax,incr)
for x1 in g1:
    beta1=sqrt(1.-(1./x1**2))
    for x2 in g2:
        beta2=sqrt(1.-(1./x2**2))
#exact merger beta value (squared)
        beta_m_sqr = ((x1*mu*beta1 + x2*mu*beta2)/(x1*mu+x2*mu))**2
#3:exact merger BLF
        g_m = 1./sqrt((1.-beta_m_sqr))
#4:approximate merger BLF
        g_m_app = sqrt((mu*x1+mu*x2)/((mu/x1)+(mu/x2)))
#5:dynamical efficieny with exact BLF
        e_d = (((x1-g_m)*mu)+((x2-g_m)*mu)+(g_m*eta)+(g_m*eta))/((x1*mu)+(x2*mu))
#6:dynamical efficieny with approximate BLF
        e_d_2=(((x1-g_m_app)*mu)+((x2-g_m_app)*mu)+(g_m_app*eta)+(g_m_app*eta))/((x1*mu)+(x2*mu))
#7:dynamical efficiency with exact BLF; incorrect internal Energy formula
        e_d_app = (((x1-g_m)*mu)+((x2-g_m)*mu)+(eta)+(eta))/((x1*mu)+(x2*mu))
#8:dynamical efficieny with approximate BLF; incorrect int Ene formula
        e_d_app_2 = (((x1-g_m_app)*mu)+((x2-g_m_app)*mu)+(eta)+(eta))/((x1*mu)+(x2*mu))
        s1=str(x1)
        s2=str(x2)
        s3=str(g_m)
        s4=str(g_m_app)
        s5=str(e_d)
        s6=str(e_d_2)
        s7=str(e_d_app)
        s8=str(e_d_app_2)
        datafile.write(s1+"\t"+s2+"\t"+s3+"\t"+s4+"\t"+s5+"\t"+s6+"\t"+s7+"\t"+s8+"\n")
    datafile.write("\n")
datafile.close()


