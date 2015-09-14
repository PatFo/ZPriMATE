#!/usr/bin/python
import os
import sys
import numpy as np
from matplotlib import pyplot as plt



# Switches to enable plotting of cross sections
plot_tot = 1
plot_sm  = 1
plot_int = 0
plot_sig = 1



#       This scrip needs the absolute path of the datafile as command line arguments 
#       and then plots the three column data.
f= sys.argv[1]
ZPMSYS =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))  

  
#x, y1, y2 = np.genfromtxt(f, delimiter='\t\t').transpose()

# Empty containers for data
x  = []
y1 = []
y2 = [] 
y3 = [] 
y4 = []

# Labels for plotting
lx  = ''
ly1 = ''
ly2 = ''
ly3 = '' 
ly4 = ''

# Read data file
with open(f, 'r') as openedFile :
    lno=0
    for line in openedFile:
        tab = line.split()
        #First line are labels 
        if lno==0:
            lx =tab[0]
            ly1=tab[1]
            ly2=tab[2]
            ly3=tab[3]
            ly4=tab[4]
        else:
            x.append( np.double(tab[0]))
            y1.append(np.double(tab[1]))
            y2.append(np.double(tab[2]))
            y3.append(np.double(tab[3]))
            y4.append(np.double(tab[4]))
        lno +=1

# Create plots
fig=plt.figure()
if plot_tot:
    plt.semilogy(x, y1, 'r',label=ly1)
if plot_int:
    plt.semilogy(x, y2, 'y',label=ly2)
if plot_sig:
    plt.semilogy(x, y3, 'g',label=ly3)
if plot_sm:
    plt.semilogy(x, y4, 'b',label=ly4)

plt.grid(True)
plt.ylabel(r'd$\sigma$/d$m^2$ [fb/GeV$^2$]')
plt.xlabel(r'$m_{inv}$ [GeV]')
plt.title('Differential Partonic Cross Section')
plt.legend()
#Include ZPriMATE logo on plot
im = plt.imread(ZPMSYS + "/icons/logo_small.png")
ax = plt.axes([0.13,0.1, 0.17, 0.17], frameon=False)  # Change the numbers in this array to position your image [left, bottom, width, height])
ax.imshow(im)
ax.axis('off')  # get rid of the ticks and ticklabels
    
f=f[:-3]+"pdf"
plt.savefig(f,dpi=300)