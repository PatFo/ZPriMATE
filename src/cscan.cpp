#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <ctime>


#include "pheno.h"
#include "xsec.h"
#include "spectrum_analysis.h"


using namespace std;


//Linear binning scheme
struct LinBin{
  double operator() (double low)
  {
    return low + 20;
  }
};



struct LogBin{
  ///Binning from Atlas paper: arXiv 1405.4123v
  const double offset = log(4208)-log(3934);
  double operator() (double low)
  {
    return exp(log(low)+offset);
  }
};



double gaussian(double mu, double x)
{
  double sigma = 0.02*mu;
  return 1./(sqrt(2 *M_PI)*sigma) * exp(-1*pow(mu-x,2)/(2*sigma*sigma));
}




int main(int argc, char** argv){
  
    char pdfset[] = "/remote/pi104a/foldenauer/local/MSTW/Grids/mstw2008lo.68cl.21.dat";
  //Initialize SSM
//   pheno::ZpModel ssm(91.1876);  
  for(double mzp =1000; mzp<=1000; mzp+=500)
  {
    pheno::ZpModel ssm(mzp);
    printf("Relative width of Zp in SSM: %g%%\n",ssm.wzp_()/ssm.mzp_()*100);
    
    pheno::HadronXSec ssmxsec(&ssm.mu, &ssm, pdfset);
//     ssmxsec.set_monte_calls(100000);
    printf("Total cross section: %g\n", ssmxsec.zpXsec(5, 8000, 0.01));
    
    double mqq = 90;
    printf("Partonic cross section for d d~ -> A -> mu+ mu- @ %g GeV: %g\n", mqq, ssmxsec.dxsec->sigGam(mqq) );
    printf("Partonic cross section for d d~ -> Z -> mu+ mu- @ %g GeV: %g\n", mqq, ssmxsec.dxsec->sigZ(mqq) );
    printf("Partonic cross section for d d~ -> Zp -> mu+ mu- @ %g GeV: %g\n", mqq, ssmxsec.dxsec->sigZp(mqq) );
    printf("Partonic cross section for d d~ -> (A,Z) -> mu+ mu- @ %g GeV: %g\n", mqq, ssmxsec.dxsec->sigGamZ(mqq) );
    printf("Partonic cross section for d d~ -> (A,Zp) -> mu+ mu- @ %g GeV: %g\n", mqq, ssmxsec.dxsec->sigGamZp(mqq) );
    printf("Partonic cross section for d d~ -> (Z,Zp) -> mu+ mu- @ %g GeV: %g\n", mqq, ssmxsec.dxsec->sigZZp(mqq) );
    printf("\nTotal cross section for d d~ -> mu+ mu- @ %g GeV: %g\n", mqq, ssmxsec.dxsec->sigTot(mqq));
    
    //Output Histogram
    struct timeval tv;
    double fb2pb = 0.001;
    pheno::HistWriter<LogBin, pheno::HadronXSec> hist(&ssmxsec, &pheno::HadronXSec::zpXsec, NULL);
    char histf[] = "/scratch/foldenauer/data/xscan/hist0.dat";
    gettimeofday(&tv, NULL);  
    double t0=tv.tv_sec+(tv.tv_usec/1000000.0);     
    hist.writeHist( 40, 4500, 0.05, histf, fb2pb);    
    
    pheno::HistWriter<LogBin, pheno::HadronXSec> hist2(&ssmxsec, &pheno::HadronXSec::zpXsec, &gaussian);
    char histf2[] = "/scratch/foldenauer/data/xscan/hist_smear.dat";
    gettimeofday(&tv, NULL);
    double t0b=tv.tv_sec+(tv.tv_usec/1000000.0);     
    printf("Writing histogram took: %g s\n", t0b-t0);
    hist2.writeHist( 40, 4500, 0.05, histf2, fb2pb); 
    
    gettimeofday(&tv, NULL);  
    double t1=tv.tv_sec+(tv.tv_usec/1000000.0);   
    printf("Writing histogram took: %g s\n", t1-t0b);
    printf("\nIntegral yields: %g\n", ssmxsec.zpXsec(50, 2*mzp, 0.1));
    gettimeofday(&tv, NULL);  
    double t2=tv.tv_sec+(tv.tv_usec/1000000.0); 
    printf("Execution time: %g s\n", t2-t1);
    printf("\nSmeared integral yields: %g\n", ssmxsec.zpXsec(50, 2*mzp, 0.1, &gaussian));
    gettimeofday(&tv, NULL);  
    double t3=tv.tv_sec+(tv.tv_usec/1000000.0); 
    printf("Execution time: %g s\n", t3-t2);
    
    //Cross check with mzp=mz --> Should yield the SM partonic cross sections!
    pheno::PartonXSec emu(&ssm.el, &ssm.mu, &ssm);
    double mee =   91;
    printf("Leptonic cross section for e+ e- -> A -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigGam(mee));
    printf("Leptonic cross section for e+ e- -> Z -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigZ(mee));
    printf("Leptonic cross section for e+ e- -> Zp -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigZp(mee));
    printf("Leptonic cross section for e+ e- -> (A,Z) -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigGamZ(mee));
    printf("Leptonic cross section for e+ e- -> (A,Zp) -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigGamZp(mee));
    printf("Leptonic cross section for e+ e- -> (Z,Zp) -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigZZp(mee));
    
    printf("\nTotal cross section for e+ e- -> mu+ mu- @ %g GeV: %g\n", mee, emu.sigTot(mee));
  }
  
  
  
  
  
  
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