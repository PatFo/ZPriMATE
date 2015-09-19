import misc
import sys
from os.path import join as pjoin
from os.path import expanduser
argc = len(sys.argv)

if argc==2:
    try:
        procID=int(sys.argv[1])
    except:
        sys.exit("Didn't recognize procID %s"%sys.argv[1])
else:
    sys.exit("Please supply procID.")

    
if procID==1:
    suffix="EE"
elif procID==2:
    suffix="MM"
else:
    sys.exit("Didn't recognize procID %s"%str(procID))
    
directory=pjoin(expanduser("~"),"KINETIC_MIXING_"+suffix)

misc.plotBisectResult(
    pjoin(directory,"limits"),
    None,
    pjoin(directory,"limitsTotBisection.pdf"),
    logo=False # Logo is in conflict with legend
)

misc.plotBisectContour(directory,pjoin(directory,"CountourBisection.pdf"))
