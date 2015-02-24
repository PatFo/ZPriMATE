#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "model.h"
#include "pheno.h"
#include "xsec.h"
#include "read_config.h"


using namespace std;

int main(int argc, char** argv){
  
  pheno::up up;
  
  std::cout<<"The pheno-up parameters:"
  <<"\n\tPDG:\t\t"<<up.get_pdg()
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
  
  up.set_vecc( new fundamental::vcoeff(up, par));

  
  cout<<"\nUp-quark vector couplings:\n-------------------"
  <<"\nqzl="<<upco.q_zl<<"\t up.qzl="<<up.vecc().q_zl
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
  
  pheno::zpmodel m(argv[1]);
  
  cout<<"\nsw2="<<m.sw2_()
  <<"\ng1="<<m.g1_()
  <<"\ng2="<<m.g2_()
  <<"\naew="<<m.aew_()
  <<"\ne="<<m.e_()
  <<"\nxi="<<m.xi_()
  <<"\nvev="<<m.vev_()
  <<"\nfef="<<m.fef_()
  <<"\ngz="<<m.gz_()
  <<"\ntan chi="<<tan(m.mixing_())
  <<"\nmzp="<<m.mzp_()
  <<endl;
  
  double w =m.wzp_();
  w=2*w;


  cout<<"nu_el parameters\nNc="<<m.ne.Nc()
  <<"\nmass="<<m.ne.get_mass()
  <<"\nqzpl="<<m.ne.vecc().q_zpl
  <<"\nqzpr="<<m.ne.vecc().q_zpr
  <<"\nyl="<<m.ne.vecc().get_hypl()
  <<"\nyr="<<m.ne.vecc().get_hypr()
  <<endl;
  
  pheno::PartonXSec xsec(&m.u, &m.mu, &m);
  
  cout<<"CrossX of pdg="<<xsec.pdg_in()<<" is "<<xsec.sigTot(1000.)<<endl;
  
  
  char pdfset[] = "/remote/pi104a/foldenauer/local/MSTW/Grids/mstw2008nnlo.00.dat";
  pheno::HadronXSec hsec(&m.mu, &m, pdfset);
//   cout<<"Total hadronic cross section for E=1000GeV:  "<<hsec.sigTotal(1000.)<<endl;
//   std::printf("%10s|%10g|%10g|%10g|%10g\n","Cross Sec", hsec.sigSM(1000.), hsec.sigInt(1000.), hsec.sigSignal(1000.), hsec.sigTotal(1000.));
  
  //Plotting the partonic cross section
  std::ofstream outf("sample_data.dat");
  float low(5), high(1.5*m.mzp_());
  float step = (high-low)/500;
  
  for(float E=low; E<high;E+=step)
  {
//     outf<<E<<"\t\t"<<xsec.sigTot(E)<<"\t\t"<<xsec.sigSM(E)<<"\n";
    outf<<E<<"\t\t"<<hsec.sigTotal(E)<<"\t\t"<<hsec.sigSM(E)<<"\n";
//     cout<<E<<"\t\t"<<xsec.sigTot(E)<<"\t\t"<<xsec.sigInt(E)<<"\n";
    
  }
  outf.close();
  
  cout<<"File written.\n";

  return 0;
}