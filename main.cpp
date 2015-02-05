#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "model.h"


using namespace std;

int main(int argc, char** argv){
  
  fundamental::fermion up(1, 1, 2.3, 2./3, 1., 1.);
  
  std::cout<<"The em charge is "<< up.get_emcharge();
  
  return 0;
}