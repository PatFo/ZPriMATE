#!/usr/bin/python
import os
import sys
import subprocess
import tempfile
import time
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import getopt

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
odir = ""
limFile = ""
effFile = ""
logname="zprimate.log"
eventFiles=[]
tmpFile="" # If nothing is set the file is generated automatically
plot=False

debug=False
force=False

def startZPriMATE(
    settingsFile,
    ZPMSYS=os.path.split(os.path.dirname(os.path.realpath(__file__)))[0], # ZPDIR
    CORE="./src/core"
):
  global tmpFile
  global force
  # If temporary file is not explicitly set, open a random one
  if tmpFile=="":
    # First, file descriptor is not used
    fdTMP, tmpFile = tempfile.mkstemp()
  #Cross section calculation

  os.chdir(ZPMSYS)
  cmd=[CORE,settingsFile,tmpFile]
  if force:
    cmd.append("--force")

  print "\nCalculating cross section ...\n"
  if not debug:
    core = subprocess.Popen(cmd,stderr=subprocess.PIPE)
  else:
    core = subprocess.Popen(cmd)
  (core_out, core_err) = core.communicate()
  return core_out, core_err, core.returncode


def getDirs():
  global tmpFile
  global odir
  global limFile
  global effFile
  global eventFiles
  with open(tmpFile,'r') as openedFile:
    fno = 0
    for line in openedFile:
      #Get output directory
      if fno ==0 :
        odir = line.split()[0]
        #Get limits directory
      elif fno == 1 :
        limFile = line.split()[0]
      #Get efficiencies
      elif fno == 2 :
        effFile = line.split()[0]
      #Operate on event files
      else:
        eventFiles.append(line.split()[0])
      fno += 1

def plotEvents(
    eventFile,
    logname,
    plotExec="./bin/plot_signal"
):
  
  global plot
  if not plot:
    return
  logFile = os.path.join(odir,logname)
  if os.path.exists(logFile):
    with open(logFile, 'a') as log:
      log.write("Plotting %s\n"%eventFile)
  plot = subprocess.Popen(["%s %s" %(plotExec, eventFile)], stdout=None, shell=True)


def combineEvents(combinedFile,eventFiles=eventFiles):
  xl   = []
  xr   = []
  nev  = []
  #Check for plot file and generate plots
  eventNum=0
  for event in eventFiles:
    eventNum+=1
    plotEvents(event,logname)
    #Combine events for global analysis
    with open(event, 'r') as eventfile:
      lno = 0
      for line in eventfile:
        data  = line.split()
        assert len(data)==3
        #First read of data
        if eventNum == 1:
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
        
  eff_fit=fitEfficiencies(effFile)

  #Combine events and multiply with efficiencies     
  combined = open(combinedFile, 'w')            
  for i in range(len(xl)):
    midpoint = xl[i] + (xr[i]-xl[i])/2
    combined.write("%s\t%s\t%s\n"%(xl[i], xr[i], nev[i]*eff_fit(midpoint)))
  combined.close()
  
  plotEvents(combinedFile,logname)


  
def fitEfficiencies(effFile):
  x_eff = []
  y_eff = []

  with open(effFile, 'r') as openedFile :
      for line in openedFile:
          tab = line.split()
          x_eff.append( np.double(tab[0]) )
          y_eff.append(np.double(tab[1]))
  #Interpolate efficiencies (Spline fit)
  eff_fit = InterpolatedUnivariateSpline(x_eff , y_eff , k=1)
  return eff_fit

def main(settingsFile, options):
  global tmpFile
  global odir
  global limFile
  global effFile
  global eventFiles
  print logo

  #Check whether inputfile was given with absolute path
  settingsFile = os.path.abspath(settingsFile)

  (core_out,core_err,core_returnval) = startZPriMATE(settingsFile)
  
  # Read out directories from temporary file
  getDirs()

  with open(os.path.join(odir,logname), 'w') as log:
    log.write("%s"%core_err)

  if not core_returnval==0:
    print """
Calculation terminated unsuccessful. Please check the logfile for more details.
"""
    return 1

  comfile = os.path.join(odir,"combined.dat")
  combineEvents(comfile)

  #Calculate best exclusion limit
  limit = subprocess.Popen(["%s %s %s" %("./bin/calc_limit", limFile, comfile) ], stdout=None, shell=True)
  limit.wait()

  print """All results can be found in:
  %s

"""%odir


def usage():
  print """
Usage:
zprimate [options] settingsFile 

Available options are:
  -h (--help)    : Display this usage message
  -f (--force)   : Force deletion of files in output directory
  -v (--verbose) : Start ZPriMATE in verbose mode
  -d (--debug)   : Display additional debugging information
"""

def getOptions(argv):

  # At this point the returning of the options is redundant but keep mechanics for now...
  global force
  global debug
  global plot

  options = dict()  
  try:
    shortOptions="phvdf"
    longOptions=["plot","help","verbose","debug","force"]
    if(len(argv)==1):
      return argv[0],[]
    
    if(len(argv)==0):
      print """
You didn't supply any arguments!"""
      opts=[("-h","")]
      args=[]
    else:

      opts,args = getopt.getopt(argv,shortOptions,longOptions)

    # If there are too many arguments recognized, check if input was given in the wrong order
    if len(args)>1:
      newInput = args[1:]
      newInput.append(args[0])
      opts,args = getopt.getopt(newInput,shortOptions,longOptions)

    # If there are still too many arguments raise exception
    if len(args)>1:
      
      print """

You supllied %s arguments without operands. Please check your input.

"""%len(args)
      print args
      opts=[("-h","")]
      
  except getopt.GetoptError as err:
    print """
Input error encountered:"""
    print err.args[0]
    opts=[("-h","")]
    
  except:
    print "Error: ", sys.exc_info()
    opts=[("-h","")]


  for opt,arg in opts:
    if opt in ("-h","--help"):
      usage()
      return "",[]
    elif opt in ("-d","--debug","-v","--verbose"):
      print "Starting in debugging mode"
      debug=True
      options["debug"]=True
      options["verbose"]=True
    elif opt in ("-f","--force"):
      options["force"]=True
      force=True
    elif opt in ("-p","--plot"):
      plot=True

  return args[0],options




if __name__=="__main__":
  settingsFile, options = getOptions(sys.argv[1:])

  if not settingsFile=="":  
    main(settingsFile,options)
