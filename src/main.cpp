// #include <stdlib.h>
// #include <ctime>
#include <fstream>
#include <math.h>
#include <stdexcept>
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/time.h>
#include <cstring>

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

// Intermediate dijet resolution, i.e. variance of gaussian
// This includes only smearing but nothing else
double dijet_resolution(double x)
{
  double 
    a = 0.00109900233544, 
    b = 1.1507949966, 
    c = 47.2285104911,
    resolution = x * sqrt( a + b / x + c / pow(x,2) )
    ;
  return resolution;
}

// Complete shape function for dijets
// parameter mu: central value, x: position to evaluate function
double dijet_shape(double mu, double x)
{
  // Preliminary values
  double
    alpha = 0.3,
    n = 10.0,
    norm = 1.0,
    t,a,b
    ;
  // Relative distance to central value
  t = ( x - mu )/dijet_resolution(mu);

  // Determine if power law or gaussian
  if ( t >= alpha ) {
    return norm * exp(-0.5 * pow(t,2) );
  } else {
    a = pow( n/alpha ,n ) * exp(-0.5*pow(alpha,2));
    b = n/alpha - alpha;
    return norm*(a / pow(b-t , n));
  }
}

// Check if file exists

inline bool fileExists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

// 
void usage(char** argv) {
  cout <<std::endl;
  cout << " Usage: " << argv[0] << " inFile [tmpFile]" << std::endl;
  cout << std::endl;
  cout << " inFile:\t\t\tSettings file with run specific parameters and paths" << std::endl;
  cout << " tmpFile (optional):\tTemporary file for communication to python process." <<std::endl;
  cout << "\t\t\tIf non is given, a default is chosen." << std::endl;
  cout <<std::endl;
}


//MAIN
//##########################


int main(int argc, char** argv){
  
  if(argc<2 || argc>3){
    usage(argv);
    return 1;
  }else if(!fileExists(argv[1])) {
    cout << "Input file " << argv[1] << " doesn't exist or cannot be opened." << std::endl;
    return 1;
  }

  //Extract settings for run given in start file 
  settings input(argv[1]);
  

  // Create communication file for python input in TMP
  // If file is given specifically use this one instead
  std::string ofname;
  if (argc == 2) { 
    char const * tmp = getenv("TMPDIR");
    if (tmp == 0)
      tmp = "/tmp";
    ofname = tmp;
    ofname.append("/hist");

  } else if(argc==3) {
    ofname = argv[2];
  }

  std::ofstream of(ofname.c_str());
  of<<input.odir()<<std::endl; //Specify the output directory
  of<<input.limdir()<<std::endl; //Specify the limits directory
  of<<input.efficiencies()<<std::endl; //Specify the effficiency file
  
  //Allocate pointer to store model
  pheno::ZpModel * model;
  struct timeval tv;

  //Start timing
  gettimeofday(&tv, NULL);  
  double t0=tv.tv_sec+(tv.tv_usec/1000000.0);     
  
  
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

  char * pdfset = new char [input.pdfset().length()+1];
  std::strcpy(pdfset, input.pdfset().c_str());
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
    fptr = 0x0;
    cs.push_back( new pheno::HadronXSec(&(model->u), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->d), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->c), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->s), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->b), model, pdfset, input.ebeam()) );
  }
  else throw std::runtime_error("ERROR: Invalid process id.\nCurrently available:\n\n\t0 = jets\n\t1 = e+ e-\n\t2 = mu+ mu-\n\t3 = tau+ tau-\n");
  
  

  
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
