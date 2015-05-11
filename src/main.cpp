// #include <stdlib.h>
#include <stdio.h>
#include <math.h>
// #include <sys/time.h>
// #include <ctime>

#include "input.h"
#include "pheno.h"
// #include "xsec.h"
// #include "spectrum_analysis.h"


using namespace std;



//Smearing function
double gaussian(double mu, double x)
{
  double sigma = 6.41632793049e-05*mu*mu + 0.0224623026794*mu + 3.75054782287;
  return 1./(sqrt(2 *M_PI)*sigma) * exp(-1*pow(mu-x,2)/(2*sigma*sigma));
}


//Functions that have verbosity flag output more information 
bool verbose=true;



//MAIN
//##########################


int main(int argc, char** argv){
  
  //Extract settings for run given in start file 
  settings input(argv[1], verbose);
  
  //Allocate pointer to sotre model
  pheno::ZpModel * pmodel;
  
  //Check whether SSM is to be used or if a model file was given
  if(input.use_ssm())
  {
    pmodel = new pheno::ZpModel(input.mzssm());
  }
  else
  {
    pmodel = new pheno::ZpModel(input.model().c_str());
  }

  
  printf("Main")  ;
  
  //Free memory 
  delete pmodel;
  
  return 0;
}