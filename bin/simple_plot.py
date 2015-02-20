#!/usr/bin/python
import sys
import numpy as np
from matplotlib import pyplot as plt

#       This scrip needs the absolute path of the datafile as command line arguments 
#       and then plots the three column data.


f= sys.argv[1]
  
x, y1, y2 = np.genfromtxt(f, delimiter='\t\t').transpose()

fig=plt.figure()
plt.semilogy(x, y1, 'r',label='Total')
plt.semilogy(x, y2, 'b',label='SM')
#plt.semilogy(x, y3, 'g',label='SM + Signal')
plt.grid(True)
plt.ylabel('CrossX [pb]')
plt.xlabel('Ecm [GeV]')
plt.title('Total CrossX')
plt.legend()
    
f=f[:-3]+"pdf"
plt.savefig(f)