#!/usr/bin/python
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter

def contourLogPlot(lobs,masses,mixings,outFile,ZPMSYS):
    #Contour plot of limits

    #levels = np.linspace(0, np.amax(lobs) , 40)
    loglevels = np.logspace(np.log(np.amin(lobs)), np.log(np.amax(lobs)) , 20)
    l_f = LogFormatter(10, labelOnlyBase=False)

    fig, axs = plt.subplots(1,1)
    cs = axs.contourf(masses, mixings, lobs, 
                      norm = LogNorm(),
                    levels=loglevels)
    cs2 = axs.contour(masses, mixings, lobs, [1], colors='k')
    #bar = fig.colorbar(cs, ax=axs, format="%.2f")
    bar = fig.colorbar(cs, ticks=loglevels, format=l_f)
    bar.set_label('R-Value')
    
    plt.title('Contour plot for kinetic mixing')
    axs.set_yscale("log") 
    plt.xlabel(r"$M_{Z^\prime}$ [GeV]")
    plt.ylabel(r'Kin. mixing $\chi$')
    #plt.legend()
    ##Include ZPriMATE logo on plot
    im = plt.imread(ZPMSYS + "/icons/logotype.png")
    ax = plt.axes([0.53,0.1, 0.2, 0.2], frameon=False)  # Change the numbers in this array to position your image [left, bottom, width, height])
    ax.imshow(im)
    ax.axis('off')  # get rid of the ticks and ticklabels
    plt.savefig(outFile, dpi=300) 
