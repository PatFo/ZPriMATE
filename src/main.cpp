#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "model.h"
#include "pheno.h"
#include "read_config.h"


using namespace std;

int main(int argc, char** argv){
  
  pheno::up up;
  
  std::cout<<"The pheno-up parameters:"
  <<"\n\tFamily:\t\t"<<up.get_family()
  <<"\n\tIsospin:\t"<<up.get_iso3()
  <<"\n\tMass:\t\t"<<up.get_mass()
  <<"\n\tE.m. charge:\t"<<up.get_emcharge()
  <<"\n\tQxl:\t\t"<<up.get_xlcharge()
  <<"\n\tQxr:\t\t"<<up.get_xrcharge()
  <<endl;
  
//   up.change_mass(9.5e-2);
//   up.update_emcharge(4./7);
//   up.update_xlcharge(1./2);
//   up.update_xrcharge(-1./2);
//   
//   std::cout<<"\n\nUpdated pheno-up parameters:\n"
//   <<"\n\tFamily:\t\t"<<up.get_family()
//   <<"\n\tIsospin:\t"<<up.get_iso3()
//   <<"\n\tMass:\t\t"<<up.get_mass()
//   <<"\n\tE.m. charge:\t"<<up.get_emcharge()
//   <<"\n\tQxl:\t\t"<<up.get_xlcharge()
//   <<"\n\tQxr:\t\t"<<up.get_xrcharge()
//   <<endl;
  
  
  pheno::muon mu;
  cout<<mu.get_iso3()<<endl;
  
  fundamental::bsm_parameters par(0.1, 1500);
  fundamental::vcoeff upco(up, par);
  fundamental::vcoeff muco(mu, par);
  
  cout<<"\nUp-quark vector couplings:\n-------------------"
  <<"\nqzl="<<upco.q_zl
  <<"\nqzr="<<upco.q_zr
  <<"\nqzpl="<<upco.q_zpl
  <<"\nqzpr="<<upco.q_zpr
  <<endl;
  
  cout<<"\nMuon vector couplings:\n-------------------"
  <<"\nqzl="<<muco.q_zl
  <<"\nqzr="<<muco.q_zr
  <<"\nqzpl="<<muco.q_zpl
  <<"\nqzpr="<<muco.q_zpr
  <<endl;
  
  conf_reader reader(argv[1]);
  
  dict parameters=reader.get_config();
  
  cout<<(parameters["up"])["cxl"];
  
  return 0;
}