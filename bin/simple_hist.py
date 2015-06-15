#!/usr/bin/python
import os
import sys
import numpy as np
from matplotlib import pyplot as plt

ZPMSYS =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))  

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
plt.ylabel('Events')
plt.xlabel(r'$m_{inv}$ [GeV]')
plt.title('Predicted number of events')
plt.legend()

#Include ZPriMATE logo on plot
im = plt.imread(ZPMSYS + "/icons/logotype.png")
ax = plt.axes([0.13,0.1, 0.2, 0.2], frameon=False)  # Change the numbers in this array to position your image [left, bottom, width, height])
ax.imshow(im)
ax.axis('off')  # get rid of the ticks and ticklabels


f=f[:-3]+"pdf"
plt.savefig(f,dpi=300) 