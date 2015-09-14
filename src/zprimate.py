#!/usr/bin/python
import os
import sys
import subprocess
import tempfile
import time
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import fcntl


logoSmall="""
 _______     _ __  __   _ _____ ___ 
|_  / _ \_ _(_)  \/  | /_\_   _| __|
 / /|  _/ '_| | |\/| |/ _ \| | | _| 
/___|_| |_| |_|_|  |_/_/ \_\_| |___|            
"""
logo="""
 _________       _ __  __    _  _____ _____ 
|__  /  _ \ _ __(_)  \/  |  / \|_   _| ____|
  / /| |_) | '__| | |\/| | / _ \ | | |  _|  
 / /_|  __/| |  | | |  | |/ ___ \| | | |___ 
/____|_|   |_|  |_|_|  |_/_/   \_\_| |_____|
"""
print logo

ZPMSYS = os.path.dirname(os.path.realpath(__file__))
#print ZPMSYS
COREX = "./src/core"
PLTX = "./bin/plot_signal.py"
LIMX = "./bin/calc_limit.py"


#Get input file
fullfile = sys.argv[1]

#Check whether inputfile was given with absolute path
if fullfile[0] !='/':
  fullfile = os.getcwd() + "/" +   fullfile

# Create unique temporary file for communication
# Filedescriptor and full path to file
fdTMP, inp = tempfile.mkstemp()


#Cross section calculation
#os.chdir(ZPMSYS)

print "\nCalculating cross section ...\n"
core = subprocess.Popen([COREX,fullfile,inp],stderr=subprocess.PIPE)

# Catch and wait for output and error if send to PIPE
(core_out, core_err) = core.communicate() 

       
    
# Write Core output
# Keep file open for later writing access
# File is closed at the end

logname="zprimate.log"
xl   = []
xr   = []
nev  = []
odir = ""
limdir = ""
eff = ""

#Check for plot file and generate plots
#print "Generating plots from temporary file %s:"%inp

with open(inp,'r') as openedFile:
  fno = 0
  for events in openedFile:
    #Get output directory
    if fno ==0 :
      odir = events.split()[0]
    #Get limits directory
    elif fno == 1 :
      limdir = events.split()[0]
    #Get efficiencies
    elif fno == 2 :
      eff = events.split()[0]
    #Operate on event files
    else:
      events = events.split()[0]
      with open(os.path.join(odir,logname), 'w') as log:
        log.write("Plotting %s\n"%events)
      plot = subprocess.Popen(["%s %s" %(PLTX, events)], stdout=None, shell=True)
      #Combine events for global analysis
      with open(events, 'r') as eventfile:
        lno = 0
        for line in eventfile:
          data  = line.split()
          assert len(data)==3
          #First read of data
          if fno == 3:
            xl.append( np.double(data[0]))
            xr.append( np.double(data[1]))
            nev.append( np.double(data[2]))
          #Append further data
          else:
            # print xl[lno], data[0]
            assert xl[lno] == np.double(data[0])
            assert xr[lno] == np.double(data[1])
            nev[lno] += np.double(data[2])
          lno += 1
    fno += 1          
           
# Remove temporary file
os.remove(inp)

with open(os.path.join(odir,logname), 'a') as log:
  log.write("%s"%core_err)


print """Calculation finished.

You can find the results in:
%s"""%odir

#FIT EFFICIENCIES  
#Allocate empty data arrays
x_eff = []
y_eff = []
with open(eff, 'r') as openedFile :
    for line in openedFile:
        tab = line.split()
        x_eff.append( np.double(tab[0]) )
        y_eff.append(np.double(tab[1]))
#Interpolate efficiencies (Spline fit)
eff_fit = InterpolatedUnivariateSpline(x_eff , y_eff , k=1)

  

            
#Combine events and multiply with efficiencies            
comfile = odir + "/combined.dat"
combined = open(comfile, 'w')            
for i in range(len(xl)):
  midpoint = xl[i] + (xr[i]-xl[i])/2
  combined.write("%s\t%s\t%s\n"%(xl[i], xr[i], nev[i]*eff_fit(midpoint)))
combined.close()
with open(os.path.join(odir,logname), 'a') as log:
  log.write("Plotting %s"%comfile)
plot = subprocess.Popen(["%s %s" %(PLTX, comfile)], stdout=None, shell=True)
  
#Calculate best exclusion limit
limit = subprocess.Popen(["%s %s %s" %(LIMX, limdir, comfile) ], stdout=None, shell=True)
limit.wait()
