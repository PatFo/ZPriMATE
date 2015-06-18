// #include <stdlib.h>
// #include <ctime>
#include <fstream>
#include <math.h>
#include <stdexcept>
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/time.h>
// #include <ctime>

//ZPriMATE headers
#include "input.h"
#include "pheno.h"
#include "spectrum_analysis.h"
#include "xsec.h"


using namespace std;





//SMEARING FUNCTIONS FOR DIFFERENT PARTICLE TYPES
//-------------------------------------------------------------

//Gaussian distribution
double gaussian(double mu, double x, double sigma) 
{  
  return 1./(sqrt(2 *M_PI)*sigma) * exp(-1*pow(mu-x,2)/(2*sigma*sigma)); 
}

//Smearing function for muons
double mu_smear(double mu, double x)
{
  double sigma = 6.41632793049e-05*mu*mu + 0.0224623026794*mu + 3.75054782287;
  return gaussian(mu, x, sigma);
}


//Smearing function for electrons
double el_smear(double mu, double x)
{
  double a = 0.651869775669; double b = 0.155008630639;  double c = 0.00545431923191; 
  double sigma = mu * sqrt( pow(a/mu, 2) + pow(b/sqrt(mu), 2) + c*c );
//   double sig_lin = 0.00561922585506*mu + 1.62532781953;
  return gaussian(mu, x, sigma);
}




//MAIN
//##########################


int main(int argc, char** argv){
  
  
  //Allocate pointer to store model
  pheno::ZpModel * model;
  struct timeval tv;

  //Start timing
  gettimeofday(&tv, NULL);  
  double t0=tv.tv_sec+(tv.tv_usec/1000000.0);     
   
  
  //Extract settings for run given in start file 
  settings input(argv[1]);
  
  
  //Check whether SSM is to be used or if a model file was given
  if(input.use_ssm())
  {
    if(input.verbose()) printf("Constructing Sequential Standard Model...\n\n");
    model = new pheno::ZpModel(input.mzssm());
  }
  else
  {
    if(input.verbose()) printf("Constructing Z' Model...\n\n");
    model = new pheno::ZpModel(input.model().c_str());
  }
  
  
  //Vector containg the appropriate cross sections
  std::vector< pheno::HadronXSec*> cs;
  double (* fptr)(double, double);
  char * pdfset = (char *) input.pdfset().c_str();
  //Construct the desired final states
  if(input.proc_id() == 1) 
  {
    fptr = &el_smear;
//     fptr = NULL;
    cs.push_back( new pheno::HadronXSec(&(model->el),  model, pdfset, input.ebeam()) );
  }
  else if (input.proc_id() == 2)
  {
    fptr = &mu_smear;
    cs.push_back( new pheno::HadronXSec(&(model->mu),  model, pdfset, input.ebeam()) );
  }
  else if (input.proc_id() == 3)
  {
    fptr = &mu_smear;
    cs.push_back( new pheno::HadronXSec(&(model->tau), model, pdfset, input.ebeam()) );
  }
  else if (input.proc_id() == 0)        // JETS   ->  vector containing cross section objects for each quark
  {
    fptr = &mu_smear;
    cs.push_back( new pheno::HadronXSec(&(model->u), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->d), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->c), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->s), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->b), model, pdfset, input.ebeam()) );
  }
  else throw std::runtime_error("ERROR: Invalid process id.\nCurrently available:\n\n\t0 = jets\n\t1 = e+ e-\n\t2 = mu+ mu-\n\t3 = tau+ tau-\n");
  
  
  
//   //################################################################## TESTING  
//   pheno::HadronXSec* pcopy =  new pheno::HadronXSec(*cs[0]);
// //   #pragma omp parallel for
//   for(int i=0; i<8; ++i)
//   {
//     Signal sig;
// //     printf("Function call in iteration %i in parralel yields: %g\n",i, sig(copy,100, 1000, 0.01, 2 ));
//     printf("Function call in iteration %i in parralel yields: %g\n",i, (*cs[0]).zpXsec(200, 1000, 0.01, fptr, 2));
//   }
//   pheno::binning testbins = pheno::get_binning((char *)input.binning().c_str());
//   (const pheno::binning) testbins;
//   int len = testbins.size();
//   #pragma omp parallel for ordered firstprivate(pcopy, fptr)//, testbins)
// //   for(pheno::binning::iterator itb = testbins.begin(); itb!=testbins.end(); ++itb)
//   for(int i=0; i<len; ++i)
//   {
//     pheno::HadronXSec* pcopy2 =  new pheno::HadronXSec(*cs[0]);
//     double low = testbins.operator[](i).first;
//     double high = testbins.operator[](i).second;
// //     double low = itb->first;
// //     double high = itb->second;
//     printf("Signal prediction in bin [%g, %g]: %i \n",low, high, pcopy2->dxsec->pdg_in());    
//     printf("Signal prediction in bin [%g, %g]: %g \n",low, high, pcopy2->accuracy_goal);    
//     printf("Signal prediction in bin (iteration %i) [%g, %g] is: %g\n", i, low, high, pcopy2->zpXsec(low, high, 0.01, fptr, 2));
// //     delete pcopy2;
//   }
//   
//   
//    printf("\n ***************** Finished prediction **************\n\n");
//   
//   //################################################################## TESTING
  
  
   
  //Create communication file for python input lying at TMP
  char const * tmp = getenv("TMPDIR");
  if (tmp == 0)
    tmp = "/tmp";
  std::string ofname(tmp);
  ofname.append("/hist");
  std::ofstream of(ofname.c_str());
  of<<input.odir()<<std::endl; //Specify the output directory
  of<<input.limdir()<<std::endl; //Specify the limits directory
  of<<input.efficiencies()<<std::endl; //Specify the effficiency file
  
  //Calculating and writing the cross section
  int count=1;
  double totx =0 ;
  pheno::binning bins = pheno::get_binning((char *)input.binning().c_str());
  for(std::vector< pheno::HadronXSec*>::iterator it = cs.begin(); it!=cs.end(); ++it)
  {
    //Construct filename
    string filename = input.odir();
    filename.append("/events");
    stringstream ss;
    ss<<count;
//     printf("%s\n", ss.str().c_str());
    if(cs.size()>1) { filename.append(ss.str().c_str()); }
    filename.append(".dat");
    //Mark for plotting
    of<<filename<<std::endl;
    
    //Write signal data
    pheno::HistWriter<pheno::HadronXSec> hist( *it, &pheno::HadronXSec::zpXsec, fptr);
    hist.writeHist(&bins, input.int_acc(), (char *)filename.c_str(), input.luminosity());
    
    //Calculate total production cross section for process
    totx += (*it)->zpXsec(5, input.ebeam(), 1e-4);
    ++count;
  }
  of.close();
  
  //Write total cross section
  string crossfile(input.odir());
  crossfile.append("/totx");
  std::ofstream cross(crossfile.c_str());
  cross<<totx<<" fb";
  cross.close();
  
  
  
  //Print total execution time
  gettimeofday(&tv, NULL);
  double t1=tv.tv_sec+(tv.tv_usec/1000000.0);     
  printf("Signal calculation took %g s\n\n", t1-t0);
  
  
  
  //Construct sampling scheme for Cross section plotting
//   double stepsize=model->wz_();
  double stepsize=10;
  double min = 5;
  double max = 1.5*model->mzp_();
  
  pheno::SpectrumScanner<pheno::HadronXSec> scanner(model, cs[0]);    
  if(model->wzp_()>model->wz_())
  {
    scanner.add_interval(min, max, stepsize);
  }
  else
  {
    double offset = model->wzp_()*10;
    scanner.add_interval( min            , model->mzp_()-offset, stepsize   );
    scanner.add_interval( model->mzp_()-offset, model->mzp_()+offset, model->wzp_()/5 );
    scanner.add_interval( model->mzp_()+offset, max            , stepsize   );
  }
  //Generate the Cross Section Spectrum
  string scan(input.odir());
  scan.append("/hadron_scan.dat");
  scanner.scan((char *) scan.c_str());
  
  
  //Print scanning time
  gettimeofday(&tv, NULL);
  double t2=tv.tv_sec+(tv.tv_usec/1000000.0);     
  printf("Cross section scan took %g s\n\n", t2-t1);
  
  
  //Free memory 
  for(unsigned int pos =0 ; pos < cs.size(); ++pos)
  {
    delete cs[pos];
  }
  delete model;
  
  return 0;
}