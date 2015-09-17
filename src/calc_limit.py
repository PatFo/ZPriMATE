#!/usr/bin/python
import os
import sys
import glob
import numpy as np
import subprocess

# Input parameters
s95script = "/remote/pi104a/foldenauer/code/get_s95.py"



def table_entry(label, value):
  return "%-15s | %-5s"%(label, value)
  

""" 
Pass the directory where the lookup tables lie as 1st argument. Directory may *ONLY* contain lookup tables
Pass the model prediction file for the corresponding analysis as 2nd argument
"""
limdir    = sys.argv[1]
modeldata = sys.argv[2]


xlpred  = []
xrpred  = []
npred   = []



#Load model prediction
with open(modeldata, 'r') as openedFile:
  for line in openedFile:
    tab = line.split()
    #Avoid empty lines
    if tab != []:
      # Assert that data has 3 columns
      assert len(tab) == 3
      #print tab     ################################ DEBUG
      # Load data
      xlpred.append(np.double(tab[0]))
      xrpred.append(np.double(tab[1]))
      npred.append(np.double(tab[2]))


# Get lookup table file names
os.chdir(limdir)
limfiles= sorted(glob.glob("*"))
#print limfiles     ################################ DEBUG STATEMENT
#Setup dcitionary of  number of combined bins and filenames
limdict = []
for f in limfiles:
  words = f.split('_')
  i = np.int(words[1])
  limdict.append([i, f])
 
  
#Allocate empty list to store best exclusion
#[xl, xr, rat=nprad/exp95, exp95, obs95, npred, nbins]
best_exclusion = [0] * 7

#Loop over lookup tables
for pair in limdict:
  nbins = pair[0]
  nadd = nbins - 1
  limits = open(limdir+"/%s"%pair[1], 'r')

  #Loop over data and calculate exclusion limits  
  i =0
  for line in limits:
    #Read the lookup table line by line
    tab = line.split()
    #print tab   ################################ DEBUG STATEMENT
    
    lo = xlpred[i]
    hi = xrpred[i + nadd]
    #Iterate of nbins consecutive bins and add up the values
    pred = 0
    for n in range(nbins):
      pred += npred[i+n]
    #print lo, hi   ################################ DEBUG
    #Check that the bins are the same
    assert lo == np.double(tab[0])
    assert hi == np.double(tab[1])
    
    #Calculate exclusion ratio
    exp95 = np.double(tab[2])
    obs95 = np.double(tab[3])
    ratio = pred/exp95
    #If you find a better exclusion ratio store the values
    if ratio > best_exclusion[2]:
      best_exclusion[0] = lo
      best_exclusion[1] = hi
      best_exclusion[2] = ratio
      best_exclusion[3] = exp95
      best_exclusion[4] = obs95
      best_exclusion[5] = pred
      best_exclusion[6] = nbins
    #Increment the file line number
    i +=1
    

    
#Print result table
try:
  ratio = best_exclusion[5]/best_exclusion[4]
except:
  ratio=1e6
table  = "SUMMARY TABLE:\n%s"%('-'*25)
table += "\n" + table_entry("No. of bins", best_exclusion[6])
table += "\n" + table_entry("Best bin", [best_exclusion[0], best_exclusion[1]])
table += "\n" + table_entry("Observed S95", best_exclusion[4])
table += "\n" + table_entry("Model signal", best_exclusion[5])
table += "\n" + table_entry("Exp R-Value", best_exclusion[2])
table += "\n" + table_entry("Obs R-Value", ratio)

print "\n" + table
resultf = os.path.dirname(modeldata) + "/results"
result = open(resultf, 'w')
result.write(table)
result.close()


#If model is excluded in one bin
if best_exclusion[2] > 1:
  print "\nRESULT: Model excluded with R-value %s!\n"%ratio  
else:
  print "\nRESULT: Model NOT excluded!\n"
