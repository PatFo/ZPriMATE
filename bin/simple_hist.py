#!/usr/bin/python
import sys
import numpy as np
from matplotlib import pyplot as plt


OX = []
OY = []

f=sys.argv[1]

with open(f, 'r') as openedFile :
    for line in openedFile:
        tab = line.split()
        OX.append(float(tab[0]))
        OY.append(float(tab[2]))

        
yma=max(OY)
ymi=min(OY)
xma=max(OX)
xmi=min(OX)



fig=plt.figure()
plt.semilogy(OX, OY, 'r',label='Signal',drawstyle='steps')
#plt.minorticks_on()
plt.grid(True)
plt.ylim((ymi,yma*1.2))
plt.xlim((xmi,xma))
plt.ylabel('Cross Section')
plt.xlabel('M_inv [GeV]')
plt.title('Events')
plt.legend()

f=f[:-3]+"pdf"
plt.savefig(f) 