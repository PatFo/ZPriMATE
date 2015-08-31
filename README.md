# <img src=icons/logotype.png width=325 height=187 /> 

ZPriMATE -- Z Prime Models At Terascale Energies


#Obtaining ZPriMATE
-----------------------

First and foremost, it has to be pointed out that ZPriMATE is work in progress and so far only a minimal working 
pre-release version for dilepton final states exists. At this stage, sources are available only from a public 
git
repository at https://github.com/PatFo/ZPriMATE.git


##Prerequisites
-------------------------

ZPriMATE has been developed for the UNIX systems  Linux and OS X and hast been tested on the following systems:

- Debian 8.1 (Jessie)
- Ubuntu 14.04 LTS (Trusty Tahr)
- Ubuntu 15.04 (Vivid Vervet)
- OS X 10.9 (Mavericks)
- OS X 10.10 (Yosemite)

The program package was designed to be as self-contained as possible. However, there are some 
dependencies 
that are required for installing and running ZPriMATE:


- **C/C++ Compiler**: The package needs an installed version of the g++ or clang compiler.
- **Python 2.7**: Python 2.7.3 is the minimum required version.  Additionally the following Python packages are 
needed
  - Numpy 
  - Scipy
  - Matplotlib



During configuration of ZPriMATE, the presence of all dependencies is inquired. If 
any dependencies are missing, an error is thrown and the user is informed of the absence of the  missing prerequisite.

#Installation

##Obtaining sources with git
-------------------------

The most convenient and platform independent way to obtain a copy of ZPriMATE is to clone the git repository. 
If git is installed on the machine, it is sufficient to run the following command in a terminal
```
$ git clone https://github.com/PatFo/ZPriMATE.git
```
This will create a subdirectory called ZPriMATE in the directory where you have issued the command. 

##Obtaining sources without git
-------------------------
If git is not installed, one can obtain a copy of the master branch from the github mirror.

On Linux:
```
$ wget https://github.com/PatFo/ZPriMATE/archive/master.tar.gz -O - | tar xz
```
On OS X:
```
$ curl -Lk https://github.com/PatFo/ZPriMATE/archive/master.tar.gz | tar xz
```
Issuing one of those commands create a subdirectory called ZPriMATE-master. 
Building

Once the source files have been obtained, installation works  as usual
```
$ cd ZPriMATE(-master)
$ ./configure [--prefix=<install-path>]
$ make 
$ make install
```
The optional '--prefix' variable allows to specify an installation path. The default path is 
'/usr/local'. Depending on where you install ZPriMATE to, you might be asked to set up the program properly by 
running
```
$ source setup.sh
```

#The program structure
---------------------------------

ZPriMATE was designed as a modular program package implemented in a hybrid Python/C++ approach. The input the user has to provide essentially consists
of the model characterized by a set of model parameters and a LHC analysis that the model shall be tested 
against. The input is passed to the  C++ core application (*Core*) that calculates the semi-analytical 
cross section and turns it into a prediction of the signal events s_i. The signal prediction and the analysis are 
then passed on to a Python routine  (*Limit Calculator*) responsible for the statistical evaluation and 
the determination of the *R*-value. Based on the determined *R*-value, a model can be excluded or not. 
A second Python routine  (*Plotter*)  plots the  calculated signal prediction s_i. 
In the following section, the individual parts of ZPriMATE are described in more detail.


##Input
-------------------

ZPriMATE is a command line tool that is invoked after successful installation  by running 
```
$ zprimate <settings/file/path>
```

```
      //Minimal settings file example
      //-----------------------------

      //PROGRAM PARAMETERS
      // $VERBOSE       = 1
      // $ACC           = 1e-2
      // $ODIR          = <outdir>
      // $PDF           = mstw/grid/mstw2008lo.00.dat


      //SPECIFY MODEL FILE OR SET MASS OF Z_SSM
      // $MODEL         = example/example.conf
      $MODEL         = 2000


      //ANALYSIS PARAMETERS
      $PROC          = 1
      $EBEAM         = 8000
      $LUM           = 20.3  
      $BINS          = analyses/arXiv_1405_4123/bins.dat
      $EFFICIENCIES  = analyses/arXiv_1405_4123/eff_el.dat
      $LIMITS        = analyses/arXiv_1405_4123/el_lims
```
It needs as input parameter the path to a valid settings file, which contains the main parameters needed for running. A minimum example file is shipped with the ZPriMATE package and is 
found under\\
($ZPriMATE)/example/settings. In essence, the parameters can be grouped into three blocks. The first block 
of program parameters  concerning the 
run behavior of ZPriMATE is optional and will not be discussed here.

###Model parametrization
 
 The second block consists of a single variable called $MODEL. This variable can take two types of values:
 
- A floating point number: If a single number is specified the program automatically chooses the SSM as the 
tested model and interprets the number as the mass M_{Z'_{SSM}} of the SSM Z'. 
- A string: Alternatively the variable can be assigned a path to a full model file defining all model parameters.

```
---------------------------------------
---------------------------------------

MODEL PARAMETERS:

---------------------------------------
---------------------------------------

$GENERAL
mzp =1000 # Mass 
gx = 0.1 # Gauge coupling
chi = 1.1 # kinetic mixing
whid = 100 # hidden width
dm = 200 # bare mass mixing (not yet supported)
$END

$DOWNL
0.1 # dd~ coupling
0.1 # ss~
0.1 # bb~
0.0 # ds~ / sd~
0.0 # db~ / bd~
0.0 # sb~ / bs~
$END
```
###LHC Analysis

The last block of parameters in the settings file are the analysis parameters. The information that ZPriMATE needs for 
running is

- $PROC: The final state id (1=dielectron, 2=dimuon)
- $EBEAM: The center of mass energy $\sqrt{s}$ in GeV
- $LUM: The integrated luminosity $L$ at which data has been taken in fb$^{-1}$
- $BINS: The file containing the bins of invariant mass used in the anlaysis
- $EFFICIENCIES: The file containing acceptance $\times$ efficiency of selected final state
- $LIMITS: The directory containing the limitfiles

So far, one ATLAS analysis for 
dilepton resonances at $\sqrt{s}=8$ TeV ([arXiv](http://arxiv.org/abs/1405.4123)) has been implemented. The corresponding files are part of 
the ZPriMATE package and are stored under ($ZPriMATE)/analyses/arXiv_1405_4123}. 
