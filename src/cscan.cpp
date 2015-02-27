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
  
  
    
  //Setup SpectrumScanner
  char pdfset[] = "/remote/pi104a/foldenauer/local/MSTW/Grids/mstw2008nnlo.00.dat";
  
  pheno::HadronXSec hsec(&m.mu, &m, pdfset);
  pheno::SpectrumScanner<pheno::HadronXSec> scanner(&m, &hsec, 1);    
 
  
  
  //Construct sampling scheme for Cross section plotting
  double stepsize=m.wz_();
  double min = 5;
  double max = 1.5*m.mzp_();
  if(m.wzp_()>m.wz_())
  {
    scanner.add_interval(min, max, stepsize);
  }
  else
  {
    double offset = m.wzp_()*5;
    scanner.add_interval( min            , m.mzp_()-offset, stepsize   );
    scanner.add_interval( m.mzp_()-offset, m.mzp_()+offset, m.wzp_()/5 );
    scanner.add_interval( m.mzp_()+offset, max            , stepsize   );
  }
  
  
  
  
  //Generate the Cross Section Spectrum
  char outfile[] = "/scratch/foldenauer/data/xscan/hadron_scan.dat";
  scanner.scan(outfile);
  
    
  return 0;
}