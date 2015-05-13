// #include <stdlib.h>
#include <fstream>
#include <math.h>
#include <stdexcept>
#include <stdio.h>
#include <sstream>
#include <sys/stat.h>
// #include <ctime>

//CSCAN headers
#include "input.h"
#include "pheno.h"
#include "xsec.h"
#include "spectrum_analysis.h"


using namespace std;



//Smearing function
double gaussian(double mu, double x)
{
  double sigma = 6.41632793049e-05*mu*mu + 0.0224623026794*mu + 3.75054782287;
  return 1./(sqrt(2 *M_PI)*sigma) * exp(-1*pow(mu-x,2)/(2*sigma*sigma));
}




//MAIN
//##########################


int main(int argc, char** argv){
  
  //Extract settings for run given in start file 
  settings input(argv[1]);
  
  
  //Allocate pointer to store model
  pheno::ZpModel * model;
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
  char * pdfset = (char *) input.pdfset().c_str();
  //Construct the desired final states
  if(input.proc_id() == 1)              cs.push_back( new pheno::HadronXSec(&(model->el),  model, pdfset, input.ebeam()) );
  else if (input.proc_id() == 2)        cs.push_back( new pheno::HadronXSec(&(model->mu),  model, pdfset, input.ebeam()) );
  else if (input.proc_id() == 3)        cs.push_back( new pheno::HadronXSec(&(model->tau), model, pdfset, input.ebeam()) );
  else if (input.proc_id() == 0)        // JETS   ->  vector containing cross section objects for each quark
  {
    cs.push_back( new pheno::HadronXSec(&(model->u), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->d), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->c), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->s), model, pdfset, input.ebeam()) );
    cs.push_back( new pheno::HadronXSec(&(model->b), model, pdfset, input.ebeam()) );
  }
  else throw std::runtime_error("ERROR: Invalid process id.\nCurrently available:\n\n\t0 = jets\n\t1 = e+ e-\n\t2 = mu+ mu-\n\t3 = tau+ tau-\n");
  
  
   
  //Create communication file for python input lying at TMP
  char const * tmp = getenv("TMPDIR");
  if (tmp == 0)
    tmp = "/tmp";
  std::string ofname(tmp);
  ofname.append("/hist");
  std::ofstream of(ofname.c_str());
  of<<input.odir()<<std::endl; //Specify the output directory
  of<<input.limdir()<<std::endl; //Specify the limits directory
  
  //Calculating and writing the cross section
  int count=1;
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
    pheno::HistWriter<pheno::HadronXSec> hist( *it, &pheno::HadronXSec::zpXsec, &gaussian);
    hist.writeHist(&bins, input.int_acc(), (char *)filename.c_str(), input.luminosity());
    
    ++count;
  }
  of.close();
  
  
  //Free memory 
  for(unsigned int pos =0 ; pos < cs.size(); ++pos)
  {
    delete cs[pos];
  }
  delete model;
  
  return 0;
}