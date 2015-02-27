#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "pheno.h"
#include "xsec.h"
#include "spectrum_analysis.h"


using namespace std;

int main(int argc, char** argv){
  
  //Initialize model with model configuration file
  pheno::ZpModel m(argv[1]);
  
  
  //Construct sampling scheme for Cross section plotting
  pheno::sampling_scheme intervals;
  
  
  double stepsize=m.wz_();
  double min = 5;
  double max = 1.5*m.mzp_();
  if(m.wzp_()>m.wz_()){
    intervals.push_back(new double[3]);
    intervals.back()[0]=min;      intervals.back()[1]=max;     intervals.back()[2]=stepsize; 
  }else{
    double offset = m.wzp_()*5;
    intervals.push_back(new double[3]);
    intervals.back()[0]=min;      intervals.back()[1]=m.mzp_()-offset;     intervals.back()[2]=stepsize; 
    intervals.push_back(new double[3]);
    intervals.back()[0]=m.mzp_()-offset;      intervals.back()[1]=m.mzp_()+offset;     intervals.back()[2]=m.wzp_()/5; 
    intervals.push_back(new double[3]);
    intervals.back()[0]=m.mzp_()+offset;      intervals.back()[1]=max;     intervals.back()[2]=stepsize; 
  }
  
  
  //Setup SpectrumScanner
  char pdfset[] = "/remote/pi104a/foldenauer/local/MSTW/Grids/mstw2008nnlo.00.dat";
  
  pheno::HadronXSec hsec(&m.mu, &m, pdfset);
  pheno::SpectrumScanner<pheno::HadronXSec> scanner(&m, &hsec, 1);    
  for(pheno::sampling_scheme::iterator it = intervals.begin(); it != intervals.end(); ++it)
  {
    //Set the sampling intervals
    scanner.add_interval(*it); 
    //Free the set interval array
    delete[] *it;
  }
  
  //Generate the Cross Section Spectrum
  char outfile[] = "/scratch/foldenauer/data/xscan/hadron_scan.dat";
  scanner.scan(outfile);
  
    
  return 0;
}