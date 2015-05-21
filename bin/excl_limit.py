#!/usr/bin/python
import numpy as np
import os
import subprocess


#Specify the data file and the output directory as input paramters 
datafile  = sys.argv[1]
outputdir = sys.argv[2]

#Check whether output directory exists; if not create
if not os.path.isdir(outputdir):
  os.mkdir(outputdir)
  print "Created output directory %s"%outputdir
  

# Input parameters
python_exec = "mypython"
s95script = os.path.dirname(__file__) + "/get_s95.py"
outfile = outputdir + "/limits"



xl   = []
xr   = []
nobs = []
nexp = []
nerr = []


#Loop over analysis data
with open(datafile, 'r') as openedFile:
  lno=0
  for line in openedFile:
    tab = line.split()
    #First 6 line contain no data 
    if lno>6:
      #Avoid empty lines
      if tab != []:
        # Assert that data has 6 columns
        assert len(tab) == 6
        #print tab
        # Load data
        xl.append(np.double(tab[1]))
        xr.append(np.double(tab[2]))
        nobs.append(np.double(tab[3]))
        nexp.append(np.double(tab[4]))
        nerr.append(np.double(tab[5]))
    lno +=1
    
#Generate lookup tables for combinations of 'nbins' neighbouring bins    
MAX_BINS = len(xl)

for nadd in range(MAX_BINS):
  nbins = 1 + nadd
  print "Combining %s bin(s) to calculate exclusion limit ...\n"%nbins
  lim = open(outfile+"_%s"%nbins, 'w')
  
  #Loop over data and calculate exclusion limits  
  for i in range( len(xl) - nadd  ):
    #Construct new bin size
    lo = xl[i]
    hi = xr[i + nadd]
    #Get bin entries
    obs=0
    exp=0
    #Iterate of nbins consecutive bins and add up the values
    for n in range(nbins):
      obs += nobs[i+n]
      exp += nexp[i+n]
    #Calculate the error as the relative error of the middle bins
    mid = 0
    if nadd%2 == 0:
      mid = i + nadd/2
    else:
      mid = i + (nadd + 1)/2
    rel = nerr[mid]/nexp[mid]
    err = rel * exp
    print lo, hi, obs, exp, err
    
    #Calculate line by line the 95% exclusion limit on the signal (observed and expected)
    proc = subprocess.Popen(["%s %s %s %s %s" %(python_exec, s95script, obs, exp, err) ], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    #print "program output:\n", out     ################################ DEBUG STATEMENT
    o95 = 0 
    e95 = 0
    for row in out.split("\n"):
      items = row.split()
      #print items
      if items != []:
          if items[0]=='S95_obs:':
              #print "observed %s" %items[1]    ################################ DEBUG STATEMENT
              o95 = items[1]
          elif items[0]=='S95_exp:':
              #print "expected %s" %items[1]    ################################ DEBUG STATEMENT
              e95 = items[1]
    lim.write("%s\t%s\t%s\t%s\n"%(lo, hi, e95, o95))
  lim.close()