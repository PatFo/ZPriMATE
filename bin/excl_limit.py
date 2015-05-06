#!/usr/bin/python
import subprocess

# Input parameters
python_exec = "mypython"
s95script = "/remote/pi104a/foldenauer/code/get_s95.py"
analysis = "/remote/pi104a/foldenauer/data/xscan/dimuon.dat"
outfile=  "/remote/pi104a/foldenauer/data/xscan/exclusion_limits_test.dat"




# Outputfile
lim = open(outfile, 'w')

#Loop over analysis data
with open(analysis, 'r') as openedFile:
  lno=0
  for line in openedFile:
    tab = line.split()
    #First 6 line contain no data 
    if lno>6:
      #Avoid empty lines
      if tab != []:
        #print tab[3], tab[4], tab[5]
        #Calculate line by line the 95% exclusion limit on the signal (observed and expected)
        proc = subprocess.Popen(["%s %s %s %s %s" %(python_exec, s95script, tab[3], tab[4], tab[5]) ], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        print "program output:\n", out
        o95 = 0 
        e95 = 0
        for row in out.split("\n"):
          items = row.split()
          #print items
          if items != []:
              if items[0]=='S95_obs:':
                  print "observed %s" %items[1]
                  o95 = items[1]
              elif items[0]=='S95_exp:':
                  print "expected %s" %items[1]
                  e95 = items[1]
        lim.write("%s\t%s\t%s\t%s\n"%(tab[1], tab[2], e95, o95))
    lno +=1
    print lno
  
  
lim.close()