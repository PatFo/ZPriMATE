#!/usr/bin/python
import numpy as np
import os
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter

ZSYS =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

def contourLogPlot(lobs,masses,mixings,outFile,ZPMSYS=ZSYS):
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
    plt.close()

def plotBisectResult(lobs,whid,outFile,ZPMSYS=ZSYS):
    lobs = OrderedDict(sorted(lobs.items(), key=lambda t: t[0]))
    data=OrderedDict()
    massCount=0
    for mass in lobs:
        massCount+=1
        rValues=lobs[mass]
        chiOld=-1
        for chi in rValues:
            Robs = rValues[chi]
            # rValues is Ordered dict. Once I pass 1.0 I can extract
            # relevant mixing angle
            if Robs>1.0:
                tmp=(chiOld+chi)/2.0
                data[mass]=tmp
                break
            chiOld=chi
        else:
            # Model not excluded
            data[mass]=None
    masses = np.zeros(massCount)
    mixings = np.zeros(massCount)
    index=0
    for mass in data:
        masses[index]=mass
        mixings[index]=data[mass]
        index+=1
        
        fig, axs = plt.subplots()
    axs.plot(masses,mixings,"^-",label=str(whid))
    plt.xlabel(r"$M_{Z^\prime}$ [GeV]")
    plt.ylabel(r'Kin. mixing $\chi$')
    axs.set_ylim([1e-3,1.0])
    axs.set_yscale("log")
    axs.legend(loc='best')
    im = plt.imread(ZPMSYS + "/icons/logotype.png")
    ax = plt.axes([0.7,0.1, 0.2, 0.2], frameon=False)  # Change the numbers in this array to position your image [left, bottom, width, height])
    ax.imshow(im)
    ax.axis('off')

    if outFile:
        plt.savefig(outFile, dpi=300)
        plt.close()
    return plt


def plotBisectComplete(fileStem,widths,outFile,ZPMSYS=ZSYS):

    for wHid in widths:
        lobs = parseOutput(fileStem+str(int(100*wHid))+".dat")
        plt=plotBisectResult(lobs,wHid,None,ZPMSYS)
    plt=savefig(outFile,dpi=300)
    plt.close()
    
def parseBisectOutput(fileName):
    lobs=dict()
    rValues=dict()
    with open(fileName,'r') as inputFile:
        for line in inputFile:
            if line.startswith('$'):
                parse=True
                mass = float(line.split()[1])
                continue
            if parse:
                if len(line.strip())==0:
                    parse=False
                    rValues=OrderedDict(sorted(rValues.items(), key=lambda t: t[0]))
                    lobs[mass]=rValues
                    break
                chi=float(line.split()[0])
                Robs=float(line.split()[1])
                rValues[chi]=Robs
                
    return OrderedDict(sorted(lobs.items(), key=lambda t: t[0]))

