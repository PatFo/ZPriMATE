#!/usr/bin/python
import numpy as np
import os
import sys
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



def floatToString(number):
    return str(float(number)).replace('.','d')

"""
Plot the results for bisection run. 
'inp' argument is either file stem for width loop 
or directly 'lops' for a single plot
"""
def plotBisectResult(inp,whid,outFile,ZPMSYS=ZSYS,logo=True):
    fig, axs = plt.subplots()
    if type(whid) in [int,float]:
        whid=[whid]
    elif whid==None:
        wdir=os.path.dirname(inp)
        files = [f for f in os.listdir(wdir) if os.path.isfile(os.path.join(wdir,f))]
        whid=[]
        for f in files:
            fileName, fileExtension = os.path.splitext(f)
            if fileExtension=='.dat':
                try:
                    whid.append(float(fileName.replace('limits','').replace('d','.'))/100.0)
                except:
                    print "Automatic extraction of widths not possible."
                    raise

    whid = sorted(whid)
    maxMass=-1
    minMass=1e10
    for width in whid:
        if type(inp)!=str:
            lobs = OrderedDict(sorted(inp.items(), key=lambda t: t[0]))
        else:
            lobs = parseBisectOutput(inp+floatToString(width*100)+".dat")
        data=OrderedDict()
        massCount=0
        for mass in lobs:
            # Extract boundaries if plots are requested for samples
            # with a different number of data points
            if mass>maxMass:
                maxMass=mass
            if mass< minMass:
                minMass=mass
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
                # Set to None s.t. matlab doesn't plot this value
                data[mass]=None
        masses = np.zeros(massCount)
        mixings = np.zeros(massCount)
        index=0
        for mass in data:
            masses[index]=mass
            mixings[index]=data[mass]
            index+=1

        axs.plot(masses,mixings,"^-",label=str(width*100)+'%')
        
    plt.xlabel(r"$M_{Z^\prime}$ [GeV]")
    plt.ylabel(r'Kin. mixing $\chi$')
    
    axs.set_xlim([None,maxMass*1.02])
    axs.set_ylim([1e-3,2.0])
    axs.set_yscale("log")

    # Print logo
    if logo:
        im = plt.imread(ZPMSYS + "/icons/logotype.png")
        # Change the numbers in this array to position your image [left, bottom, width, height])
        ax = plt.axes([0.7,0.1, 0.2, 0.2], frameon=False)  
        ax.imshow(im)
        ax.axis('off')

    
    axs.legend(loc='best',title='Width in percent of mass',frameon=True)
    
    plt.savefig(outFile, dpi=300)

    plt.close()
    return plt


    
def parseBisectOutput(fileName):
    lobs=dict()
    rValues=dict()
    parse=False
    with open(fileName,'r') as inputFile:
        for line in inputFile:
            if line.startswith('$'):
                parse=True
                mass = float(line.split()[1])
                continue
            if parse:
                if len(line.split())<2:
                    rValues=OrderedDict(sorted(rValues.items(), key=lambda t: t[0]))
                    lobs[mass]=rValues
                    rValues=dict()
                    parse=False
                    continue
                chi=float(line.split()[0])
                
                Robs=float(line.split()[1])
                
                rValues[chi]=Robs
                
    return OrderedDict(sorted(lobs.items(), key=lambda t: t[0]))


def progressBar(progress,barLength=10, status = ""):
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
