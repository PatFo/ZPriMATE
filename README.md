# <img src=icons/logotype.png width=325 height=187 /> 

ZPriMATE -- Z Prime Models At Terascale Energies


Obtaining ZPriMATE
-----------------------

First and foremost, it has to be pointed out that ZPriMATE is work in progress and so far only a minimal working 
pre-release version for dilepton final states exists. At this stage, sources are available only from a public 
\texttt{git} 
repository at \url{https://github.com/PatFo/ZPriMATE.git}


Prerequisites
-------------------------

ZPriMATE has been developed for the \texttt{UNIX} systems  Linux\footnote{Tested on Debian 8.1 
(Jessie), Ubuntu 14.04 LTS (Trusty Tahr) and Ubuntu 15.04 (Vivid Vervet)} and OS X\footnote{Tested on OS X 10.9 
(Mavericks) and OS X 10.10 
(Yosemite)}. The program package was designed to be as self-contained as possible. However, there are some 
dependencies 
that are required for installing and running ZPriMATE:
\begin{itemize}
 \item \emph{C/C++ Compiler}: The package needs an installed version of the \texttt{g++} or \texttt{clang} compiler.
 \item \emph{Python 2.7}: Python 2.7.3 is the minimum required version.  Additionally the following Python packages are 
needed
\begin{itemize}
 \item Numpy \cite{scipy}
 \item Scipy \cite{scipy}
 \item Matplotlib  \cite{Hunter:2007}
\end{itemize}
\end{itemize}

During configuration of ZPriMATE, the presence of all dependencies is inquired. If 
any dependencies are missing, an error is thrown and the user is informed of the absence of the  missing prerequisite.

\subsection{Installation}
\subsubsection{Obtaining sources with \texttt{git}}

The most convenient and platform independent way to obtain a copy of ZPriMATE is to clone the \texttt{git} repository. 
If \texttt{git} is installed on the machine, it is sufficient to run the following command in a terminal
\begin{lstlisting}[%backgroundcolor = \color{lightgray},
		  language = bash,
		  basicstyle=\footnotesize\ttfamily,
		  xleftmargin = 20pt,
		  framexleftmargin = 0em]
$ git clone https://github.com/PatFo/ZPriMATE.git
\end{lstlisting}
This will create a subdirectory called \texttt{ZPriMATE} in the directory where you have issued the command. 

\subsubsection{Obtaining sources without \texttt{git}}
If \texttt{git} is not installed, one can obtain a copy of the master branch from the github mirror.
\newline
\\ \underline{On Linux}:
\begin{lstlisting}[%backgroundcolor = \color{lightgray},
		  language = bash,
		  basicstyle=\footnotesize\ttfamily,
		  xleftmargin = 0pt,
		  framexleftmargin = 0em]
$ wget https://github.com/PatFo/ZPriMATE/archive/master.tar.gz -O - | tar xz
\end{lstlisting}
\underline{On OS X}:
\begin{lstlisting}[%backgroundcolor = \color{lightgray},
		  language = bash,
		  basicstyle=\footnotesize\ttfamily,
		  xleftmargin = 0pt,
		  framexleftmargin = 0em]
$ curl -Lk https://github.com/PatFo/ZPriMATE/archive/master.tar.gz | tar xz
\end{lstlisting}
Issuing one of those commands create a subdirectory called \texttt{ZPriMATE-master}. 
\subsubsection{Building}
Once the source files have been obtained, installation works  as usual
    \begin{lstlisting}[%backgroundcolor = \color{lightgray},
		  language = bash,
		  basicstyle=\footnotesize\ttfamily,
		  xleftmargin = 20pt,
		  framexleftmargin = 0em]
$ cd ZPriMATE(-master)
$ ./configure [--prefix=<install-path>]
$ make 
$ make install
\end{lstlisting}
The optional -\texttt{-prefix} variable allows to specify an installation path. The default path is 
\texttt{/usr/local}. Depending on where you install ZPriMATE to, you might be asked to set up the program properly by 
running
    \begin{lstlisting}[%backgroundcolor = \color{lightgray},
		  language = bash,
		  basicstyle=\footnotesize\ttfamily,
		  xleftmargin = 20pt,
		  framexleftmargin = 0em]
$ source setup.sh
\end{lstlisting}

\section{The program structure}

\begin{figure}[h]
\begin{tikzpicture}[%
  auto,
  block/.style={
    rectangle,
    draw=black,
    thick,
    text width=6em,
    align=center,
    rounded corners,
    minimum height=2em
  },
  block1/.style={
    ellipse,
    draw=blue,
    thick,
    fill=blue!20,
    text width=6em,
    align=center,
%     rounded corners,
    minimum height=2em
  },
  block2/.style={
    rectangle,
    draw=black,
    thick,
    text width=5em,
    align=center,
    minimum height=5em
  },
  block3/.style={
    rectangle,
    draw=black,
    thick,
    fill=blue!20,
    text width=6em,
    align=center,
    rounded corners,
    minimum height=2em
  },
  line/.style={
    draw,thick,
    -latex',
    shorten >=2pt
  }
]
\draw (0,0) node[block1] (in1) {LHC Analysis};
\draw (0,-2) node[block1] (in2) {Model parameters $\boldsymbol\theta$};
\draw (4,-1) node[block2] (c) {Core};
\draw (9,-0.5) node[block] (limit){Limit Calculator} ;
\draw (9,-2.5) node[block] (plt) {Plotter};
\draw (limit.west)+(0.1,0.2) node[] (lim1) {};
\draw (limit.west)+(0.1,-0.2) node[] (lim2) {};
\draw (13,-0.5) node[block3] (r) {$R$-value};
\draw (13,-1.5) node[block3] (pred) {signal $\{s_i\}$};
\draw (13,-2.5) node[block3] (pl) {plots};

\draw (c.west) -- ++(0,0.25) coordinate (cuin);
\draw (c.west) -- ++(0,-0.25) coordinate (cdin);
\draw (c.east)  ++(0,0.25) coordinate (cuout);
\draw (c.east)  ++(0,-0.25) coordinate (cdout);
\draw (cuout)  ++(1,0) coordinate (intu);
\draw (cuout) -- ++(0.5,0) coordinate (intd);

\draw[-latex,thick] (in1.east) |- (cuin);
\draw[-latex,thick] (in2.east) |- (cdin);
\draw[-latex,thick] (intd) |- (lim1);
\draw[-latex,thick] (intu) |- (lim2);
\draw[-latex,thick] (intu) |- (plt.west);
\draw[-latex,thick] (intu) |- (pred.west);
\draw[-latex,thick] (cdout) |- (pred.west);
\draw[-latex,thick] (plt.east) |- (pl.west);
\draw[-latex,thick] (limit.east) |- (r.west);

\draw  (limit.west)+(-1.5,0.5) node (label) {analysis};
\draw  (plt.west)+(-1.5,-0.5) node {signal $\{s_i\}$};
\draw  (label)+(0,1.5) node {\large ZPriMATE};
\draw  (in1)+(0,1.5) node {\large Input};
\draw  (r)+(0,2) node {\large Output};

\draw[red,thick,rounded corners,fill=red!40, opacity=.3] ($(c.north west)+(-0.5,0.6)$)  rectangle ($(plt.south 
east)+(0.5,-0.6)$);
\end{tikzpicture}
\caption{Working structure of ZPriMATE}
\label{zpschematic}
\end{figure}

ZPriMATE was designed as a modular program package implemented in a hybrid \texttt{Python/C++} approach. Its 
conceptional structure is depicted in \cref{zpschematic}. The input the user has to provide essentially consists
of the model characterized by a set of parameters $\boldsymbol\theta$ and a LHC analysis that the model shall be tested 
against. The input is passed to the  \texttt{C++} core application (\textit{Core}) that calculates the semi-analytical 
cross section and turns it into a prediction of the signal events $\{s_i\}$. The signal prediction and the analysis are 
then passed on to a \texttt{Python} routine  (\textit{Limit Calculator}) responsible for the statistical evaluation and 
the determination of the $R$-value. Based on the determined $R$-value, a model can be excluded or not. 
A second \texttt{Python} routine  (\textit{Plotter})  plots the  calculated signal prediction $\{s_i\}$. \\
In the following section, the individual parts of ZPriMATE are described in more detail.


\subsection{Input}

ZPriMATE is a command line tool that is invoked after successful installation  by running 
\begin{lstlisting}[%backgroundcolor = \color{lightgray},
	      language = bash,
	      basicstyle=\footnotesize\ttfamily,
	      xleftmargin = 20pt,
	      framexleftmargin = 0em]
$ zprimate <settings/file/path>
\end{lstlisting}

\begin{figure}[h]
     \begin{lstlisting}[%backgroundcolor = \color{lightgray},
		      basicstyle=\fontsize{9}{9}\ttfamily,
		      frame = single,
		      xleftmargin = -25pt,
		      framexleftmargin = -3em]
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
    \end{lstlisting}
    \caption{Minimum example of settings file}
    \label{settings}
\end{figure}
It needs as input parameter the path to a valid settings file, which contains the main parameters needed for running. In 
\cref{settings} a minimum example of a settings file is shown, which is shipped with the ZPriMATE package and is 
found under\\
\texttt{(\$ZPriMATE)/example/settings}. In essence, the parameters can be grouped into three blocks. The first block 
of program parameters  concerning the 
run behavior of ZPriMATE is optional and will not be discussed here.

\subsubsection{Model parametrization}
 
 The second block consists of a single variable called \texttt{\$MODEL}. This variable can take two types of values:
 \begin{itemize}
  \item A floating point number: If a single number is specified the program automatically chooses the SSM as the 
tested model and interprets the number as the mass $M_{Z'_{SSM}}$ of the SSM $Z'$. 
\item A string: Alternatively the variable can be assigned a path to a full model file defining all model parameters.
 \end{itemize}



\begin{figure}[h]
\begin{lstlisting}[%backgroundcolor = \color{lightgray},
		  basicstyle=\fontsize{9}{9}\ttfamily,
		  frame = single,
		  xleftmargin = 0pt,
		  framexleftmargin = 0em]
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
\end{lstlisting}
  \caption{Excerpt from \texttt{example.conf} model file.}
  \label{modelfile}
\end{figure}


\subsubsection{LHC Analysis}

The last block of parameters in the settings file are the analysis parameters. The information that ZPriMATE needs for 
running is
\begin{itemize}
 \item \texttt{\$PROC}: The final state id (1=dielectron, 2=dimuon)
 \item \texttt{\$EBEAM}: The center of mass energy $\sqrt{s}$ in GeV
 \item \texttt{\$LUM}: The integrated luminosity $L$ at which data has been taken in fb$^{-1}$
 \item \texttt{\$BINS}: The file containing the bins of invariant mass used in the anlaysis
 \item \texttt{\$EFFICIENCIES}: The file containing acceptance $\times$ efficiency of selected final state
 \item \texttt{\$LIMITS}: The directory containing the limitfiles generated as described in \cref{sec_limits}
\end{itemize}
So far, one ATLAS analysis for 
dilepton resonances at $\sqrt{s}=8$ TeV \cite{Aad:2014cka} has been implemented. The corresponding files are part of 
the ZPriMATE package and are stored under \texttt{(\$ZPriMATE)/analyses/arXiv\_1405\_4123}. 
