// #include <stdlib.h>
#include <stdio.h>
#include <math.h>
// #include <sys/time.h>
// #include <ctime>

#include "input.h"
// #include "pheno.h"
// #include "xsec.h"
// #include "spectrum_analysis.h"


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
  
      
  settings input(argv[1]);

  printf("Main")  ;
  
  return 0;
}