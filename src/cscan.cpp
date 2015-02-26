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
  double s1[3]={5, 200, 200./20};
  double s2[3]={200., 3700, 3300./40};
  double s3[3]={3700, 4300, 600./50};
  double s4[3]={4300, 6000, 1500./20};
  
  
  pheno::PartonXSec xsec(&m.u, &m.mu, &m);
  
  pheno::SpectrumScanner<pheno::PartonXSec> pscan(&m, &xsec, 1);
  pscan.add_interval(s1);
  pscan.add_interval(s2);
  pscan.add_interval(s3);
  pscan.add_interval(s4);
  char pout[] = "/scratch/foldenauer/data/xscan/parton_scan.dat";
  pscan.scan(pout);
  
  
  
  char pdfset[] = "/remote/pi104a/foldenauer/local/MSTW/Grids/mstw2008nnlo.00.dat";
  
  pheno::HadronXSec hsec(&m.mu, &m, pdfset);
  pheno::SpectrumScanner<pheno::HadronXSec> scanner(&m, &hsec, 1);    
  scanner.add_interval(s1);
  scanner.add_interval(s2);
  scanner.add_interval(s3);
  scanner.add_interval(s4);
  char outfile[] = "/scratch/foldenauer/data/xscan/hadron_scan.dat";
  scanner.scan(outfile);
  
    
  return 0;
}