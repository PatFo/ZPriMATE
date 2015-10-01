import misc
import sys
from os.path import join as pjoin
from os.path import expanduser

def main(suffix):
    
    directory=pjoin(expanduser("~"),"KINETIC_MIXING_"+suffix)

    misc.plotBisectResult(
        pjoin(directory,"limits"),
        None,
        pjoin(directory,"limitsTotBisection"+suffix+".pdf"),
        logo=False, # Logo is in conflict with legend,
        exclude=[0.005,0.05,0.25]
    )
    
    misc.plotMixVsWidth(
        [200.0,500.0,1000.0,2000.0],
        pjoin(directory,"limits"),
        pjoin(directory,"MixingVsWidth"+suffix+".pdf"),
        fit=True
    )
    #misc.plotBisectContour(directory,pjoin(directory,"CountourBisection.pdf"))

if __name__=="__main__":
    argc = len(sys.argv)

    if argc==2:
        try:
            procID=int(sys.argv[1])
        except:
            sys.exit("Didn't recognize procID %s"%sys.argv[1])
    else:
        sys.exit("Please supply procID.")

    
    if procID==1:
        suffix=["EE"]
    elif procID==2:
        suffix=["MM"]
    elif procID==0:
        suffix=["EE","MM"]
    else:
        sys.exit("Didn't recognize procID %s"%str(procID))
    for suf in suffix:
        main(suf)
