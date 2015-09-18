import misc
import sys
from os.path import join as pjoin

directory="/home/jetter/KINETIC_MIXING_MM"

misc.plotBisectResult(
    pjoin(directory,"limits"),
    None,
    pjoin(directory,"limitsTest.pdf"),
    logo=False
)

misc.plotBisectContour(directory,pjoin(directory,"CountourTest.pdf"))
