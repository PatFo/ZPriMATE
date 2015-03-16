#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "pheno.h"
#include "xsec.h"
#include "spectrum_analysis.h"


using namespace std;


//Linear binning scheme
struct LinBin{
  double operator() (double low)
  {
    return low + 5;
  }
};



struct LogBin{
  ///Binning from Atlas paper: arXiv 1405.4123v
  double offset=log(4208)-log(3934);
  double operator() (double low)
  {
    return exp(log(low)+offset);
  }
};





int main(int argc, char** argv){
  
  //Initialize SSM
//   pheno::ZpModel ssm(91.1876);  
  pheno::ZpModel ssm(1000);
  printf("Relative width of Zp in SSM: %g%%\n",ssm.wzp_()/ssm.mzp_()*100);
  
  char pdfset[] = "/remote/pi104a/foldenauer/local/MSTW/Grids/mstw2008lo.00.dat";
  pheno::HadronXSec ssmxsec(&ssm.mu, &ssm, pdfset);
  ssmxsec.set_monte_calls(1000);
  printf("Total cross section: %g\n", ssmxsec.zpXsec(5, 8000, 0.01));
  
  
  //Output Histogram
  pheno::HistWriter<LogBin, pheno::HadronXSec> hist(&pheno::HadronXSec::totXsec, &ssmxsec);
  hist.writeHist( 40, 1500, 0.05,"/scratch/foldenauer/data/xscan/hist.dat");
  
  
  pheno::PartonXSec emu(&ssm.el, &ssm.mu, &ssm);
  double mee =   800;
  printf("Leptonic cross section for e+ e- -> A -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigGam(mee));
  printf("Leptonic cross section for e+ e- -> Z -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigZ(mee));
  printf("Leptonic cross section for e+ e- -> Zp -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigZp(mee));
  printf("Leptonic cross section for e+ e- -> (A,Z) -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigGamZ(mee));
  printf("Leptonic cross section for e+ e- -> (A,Zp) -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigGamZp(mee));
  printf("Leptonic cross section for e+ e- -> (Z,Zp) -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigZZp(mee));
  
  printf("\nTotal cross section for e+ e- -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigTot(mee));
  
  //Initialize model with model configuration file
  pheno::ZpModel m(argv[1]);    
  
  //Setup SpectrumScanner
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
    double offset = m.wzp_()*10;
    scanner.add_interval( min            , m.mzp_()-offset, stepsize   );
    scanner.add_interval( m.mzp_()-offset, m.mzp_()+offset, m.wzp_()/5 );
    scanner.add_interval( m.mzp_()+offset, max            , stepsize   );
  }
  
  
  
  
  //Generate the Cross Section Spectrum
  char outfile[] = "/scratch/foldenauer/data/xscan/hadron_scan.dat";
  scanner.scan(outfile);
  
    
  return 0;
}