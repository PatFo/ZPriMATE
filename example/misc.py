#!/usr/bin/python
import numpy as np
import os
from os.path import join as pjoin
import sys
from collections import OrderedDict
import matplotlib as mpl
# Set mpl backend. This way no running X-server is needed
# Has to be done before pyplot is imported!
mpl.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


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

def getWidthsDirectory(directory):
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory,f))]
    whid=[]
    for f in files:
        fileName, fileExtension = os.path.splitext(f)
        if fileExtension=='.dat':
            try:
                whid.append(float(fileName.replace('limits','').replace('d','.'))/100.0)
            except:
                print "Automatic extraction of widths not possible."
                raise
    return sorted(whid)

def plotMixVsWidth(masses,fileStem,outFile):
    if type(masses) in [int,float]:
        masses=[masses]
    # A lot of duplicate code incoming...
    wdir = os.path.dirname(fileStem)
    widths=sorted(getWidthsDirectory(wdir))
    
    fig, axs = plt.subplots()
    for mass in masses:
        mixings=[]
        for iWidth,width in enumerate(widths):
            lobs = parseBisectOutput(fileStem+floatToString(width*100)+".dat")

            if mass in lobs:
                mixings.append(getBestChi(lobs[mass]))
            else:
                # loop until we pass 'mass'
                # take take average between this and the last chi
                for m in lobs:
                    if m > mass:
                        chi = getBestChi(lobs[m])
                        mixings.append((chi+chiOld)/2.0)
                        break
                    else:
                        chiOld = getBestChi(lobs[m])
        axs.plot(np.multiply(widths,100.0),mixings,label='Mass %s'%str(mass))
    
    axs.legend(loc='best')
    
    plt.xlabel(r"Hidden width [%]")
    plt.ylabel(r'Kin. mixing $\chi$')
    plt.savefig(outFile,dpi=300)
    plt.close()

def plotBisectContour(directory,outFile):
    widths = getWidthsDirectory(directory)
    mixings = None
    Initialized = False
    masses=[]
    for iWidth, width in enumerate(widths):
        lobs=parseBisectOutput(pjoin(directory,"limits"+floatToString(width*100.0)+".dat"))
        # Can initilize array only after knowing nubmer of masses
        if not Initialized:
            mixings = np.zeros((len(widths),len(lobs)))
            Initialized=True
        iMass = 0
        for mass in lobs:
            mixings[iWidth,iMass]=getBestChi(lobs[mass])
            iMass+=1

        masses=[mass[0] for mass in lobs.items()]
    #levels = np.linspace(0, np.amax(lobs) , 40)
    loglevels = np.logspace(np.log(np.amin(mixings)), np.log(np.amax(mixings)) , 20)
    l_f = LogFormatter(10, labelOnlyBase=False)

    fig, axs = plt.subplots(1,1)

    cmap=mpl.pyplot.cm.jet
    widthsWoutZero=widths[1:]
    mixingsWoutZero=mixings[1:,:]
    cs = axs.contourf(masses,widthsWoutZero, mixingsWoutZero,
                      locator=ticker.LogLocator(),
                      levels=loglevels
                  )

    plt.title('Contour plot for kinetic mixing')
    axs.set_yscale("log")
    axs.set_xscale("log") 
    plt.xlabel(r"$M_{Z^\prime}$ [GeV]")
    plt.ylabel(r"Hidden width [\%]")
    inset=False
    if inset:
        axins = plt.axes([0.2,0.3,0.2,0.2])#zoomed_inset_axes(axs,1,loc=2)
        axins.contourf(masses,widths, mixings,
                       locator=ticker.LogLocator(),
                       levels=loglevels
                   )
        axins.set_xlabel('')
        axins.set_ylabel('')
        axins.set_xticklabels([])
        axins.set_yticklabels([])
        x1, x2, y1, y2 = 100.0, 200.0, 0.0, 0.01
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        #mark_inset(axs, axins, loc1=2, loc2=4
         #      , fc="none", ec="0.5"
          #     )
    #cs2 = axs.contour(masses, widths, mixings, [1e-3,1e-2,1e-1,1e-0], colors='k')
    #bar = fig.colorbar(cs, ax=axs, format="%.2f")

    bar = fig.colorbar(cs,
                       ax=axs,
                       ticks=loglevels,
                       format=l_f
                   )
    bar.set_label(r'Kin. mixing $\chi$')
    

    #plt.legend()
    ##Include ZPriMATE logo on plot
    
    plt.savefig(outFile, dpi=300)
    plt.close()

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
        whid = getWidthsDirectory(wdir)


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
            data[mass]=getBestChi(rValues)

        masses = np.zeros(massCount)
        mixings = np.zeros(massCount)
        index=0
        for mass in data:
            masses[index]=mass
            mixings[index]=data[mass]
            index+=1

        axs.plot(masses,mixings,"-",label=str(width*100)+'%')

    plt.xlabel(r"$M_{Z^\prime}$ [GeV]")
    plt.ylabel(r'Kin. mixing $\chi$')
    
    axs.set_xlim([100.0,maxMass*1.02])
    axs.set_ylim([1e-3,2.0])
    axs.set_yscale("log")
    axs.set_xscale("log")

    axs.set_xticks([100,200,500,1000,2000,3000])
    axs.set_xticklabels([100,200,500,1000,2000,3000])
    axs.minorticks_off()
    
    # Print logo
    if logo:
        im = plt.imread(ZPMSYS + "/icons/logotype.png")
        # Change the numbers in this array to position your image [left, bottom, width, height])
        ax = plt.axes([0.7,0.1, 0.2, 0.2], frameon=False)  
        ax.imshow(im)
        ax.axis('off')

    
    axs.legend(
        loc=4, # lower right 'best',
        title='Hidden width in percent of boson mass',
        frameon=True,
        ncol=2
    )
    
    plt.savefig(outFile, dpi=300)

    plt.close()
    return plt

def getBestChi(rValues):
    chiOld=-1
    for chi in rValues:
        Robs = rValues[chi]
        # rValues is Ordered dict. Once I pass 1.0 I can extract
        # relevant mixing angle
        if Robs>1.0:
            tmp=(chiOld+chi)/2.0
            return tmp
        chiOld=chi
    else:
        # Model not excluded
        # Set to None s.t. matlab doesn't plot this value
        return None

    
def parseBisectOutput(fileName):
    lobs=dict()
    rValues=dict()
    parse=False
    with open(fileName,'r') as inputFile:
        for line in inputFile:
            if line.startswith('$$'):
                continue
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
