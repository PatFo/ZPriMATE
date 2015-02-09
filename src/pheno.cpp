#define _USE_MATH_DEFINES

#include "pheno.h"
#include <iostream>
#include <cmath>

using namespace pheno;




//************************************//
//            QUARK CLASSES           //
//************************************//

///Implemenation of quark constructors

down::down(float xlc, float xrc, bool massive): fermionExt(massive, 1, -1./2, 4.8e-3, -1./3, xlc, xrc, 3) { }

up::up(float xlc, float xrc, bool massive): fermionExt(massive, 2, 1./2, 2.3e-3, 2./3, xlc, xrc, 3) { }

strange::strange(float xlc, float xrc, bool massive): fermionExt(massive, 3, -1./2, 9.5e-2, -1./3, xlc, xrc, 3) { }

charm::charm(float xlc, float xrc, bool massive): fermionExt(massive, 4, 1./2, 1.275, 2./3, xlc, xrc, 3) { }

bottom::bottom(float xlc, float xrc, bool massive): fermionExt(massive, 5, -1./2, 4.18, -1./3, xlc, xrc, 3) { }

top::top(float xlc, float xrc, bool massive): fermionExt(massive, 6, 1./2, 173.5, 2./3, xlc, xrc, 3) { }



//************************************//
//            LEPTON CLASSES          //
//************************************//


///Implemenation of lepton constructors

electron::electron(float xlc, float xrc, bool massive): fermionExt(massive, 11, -1./2, 5.11e-4, -1, xlc, xrc, 1) { }

muon::muon(float xlc, float xrc, bool massive): fermionExt(massive, 13, -1./2, 0.1057, -1, xlc, xrc, 1) { }

tauon::tauon(float xlc, float xrc, bool massive): fermionExt(massive, 15, -1./2, 1.777, -1, xlc, xrc, 1) { }

nu_el::nu_el(float xlc, float xrc): fermionExt(false, 12, 1./2, 0, 0, xlc, xrc, 1) { }

nu_mu::nu_mu(float xlc, float xrc): fermionExt(false, 14, 1./2, 0, 0, xlc, xrc, 1) { }

nu_tau::nu_tau(float xlc, float xrc): fermionExt(false, 16, 1./2, 0, 0, xlc, xrc, 1) { }




//************************************//
//            ZP-MODEL CLASS          //
//************************************//



//Constructor: The whole model is set up HERE
zpmodel::zpmodel(const char* configfile): bsm_parameters(0.1, 1500, 0) /*partial_widths(),*/  //Default values if no parameters are specified in config file
{
  //Set up list of pointers to fermions for map iteration
  flst["up"]=&u; flst["charm"]=&c; flst["top"]=&t;
  flst["down"]=&d; flst["strange"]=&s; flst["bottom"]=&b;
  
  flst["electron"]=&el; flst["muon"]=&mu; flst["tauon"]=&tau;  
  flst["nu_el"]=&nt; flst["nu_mu"]=&nm; flst["nu_tau"]=&nt;

  //Get model configuration from config file
  conf_reader reader(configfile);
  dict init = reader.get_config();
  
  std::cout<<"Applying configuration:\n";
  
  //Set up MODEL PARAMETERS
  dict::iterator it= init.find("model_parameters");
  if(it == init.end())
  {
    std::cout<<"No model parameters specified. Using default:\n";
  }
  else
  {
    std::cout<<"Found \"$"<<it->first<<"\":\n";
    set_gx( (it->second)["gx"] );
    set_mzp( (it->second)["mzp"] );
    set_mixing( (it->second)["chi"] );
  }
  std::cout<<"\ngx="<<gx_()<<"\nmzp="<<mzp_()<<"\nmixing="<<mixing_()<<"\n\n";
  
  
  //Applying FERMION CONFIGURATION:
  //Iterate over the whole fermion list and check for initialization values passed in config file
  for (fermion_list::iterator ferms=flst.begin(); ferms!=flst.end(); ++ferms)
  {
    it = init.find(ferms->first);  //Fermion label(string)
    if(it == init.end())
    {
      std::cout<<ferms->first<<": default.\n";
    }
    else
    {
      //Set new fermion parameters
      (ferms->second)->update_xlcharge( (it->second)["cxl"] );
      (ferms->second)->update_xrcharge( (it->second)["cxr"] );
      if( !(it->second)["massive"] )
      {
        (ferms->second)->change_mass(0);
      }
      
      //Initialize fermion vector couplings
      (ferms->second)->set_vecc( new fundamental::vcoeff( *(ferms->second), *this) );
      
      //Print fermion parameters after initialization
      std::cout<<ferms->first<<":\t\tcxl="<<(ferms->second)->get_xlcharge()<<"\tcxr="<<(ferms->second)->get_xlcharge()<<"\tmass:"<<(ferms->second)->m()<<"\tqgam:"<<((ferms->second)->vecc()).q_gam/e_()<<"\n";      
    }
  }  
    
  //Initialize widths to -1 ("not yet calculated")
  partial_widths=NULL;
  higgs_width=-1.;
  wzp=-1.;
}





//Calculate partial fermionic width of Zp
double zpmodel::calc_width(fundamental::fermionExt& f)
{
  double ratio = pow(f.get_mass()/mzp_(), 2);
//   std::cout<<f.get_mass();
  double kin = 1 - 4*ratio;
  return  double(f.Nc()) /(24*M_PI) * sqrt(kin) * ( (pow(f.vecc().q_zpl, 2) + pow(f.vecc().q_zpr, 2))*kin + 6*f.vecc().q_zpl*f.vecc().q_zpr*ratio );
}



//Calculate total Zp width
double zpmodel::wzp_()
{
  if(wzp==-1) //-1 means not yet calculated
  {
    wzp=0;
    for(fermion_list::iterator it=flst.begin(); it!=flst.end(); ++it)
    {
      wzp+=calc_width( *(it->second) );
      std::cout<<calc_width( *(it->second) )<<std::endl;
    }
    //Higgs width
    if(mixing_()!=0)
    {
      double mh = 125.36;
      double geff = (gz_()*sin(xi_()) + g1_()*cos(xi_())*tan(mixing_()))*(gz_()*cos(xi_()) - g1_()*sin(xi_())*tan(mixing_()));
      double pref = pow(geff*vev_(),2)/(96*sqrt(2)*M_PI*pow(mzp_(),3));
      double fac1 = 2 + pow( (mh*mh - (mz_()*mz_()+ mzp_()*mzp_()))/(2*mz_()*mzp_()) ,2);
      double fac2 = sqrt( pow(mzp_()*mzp_()+mz_()*mz_() -mh*mh, 2) - pow(2*mz_()*mzp_(),2) );
      higgs_width= pref*fac1*fac2;
      wzp+=higgs_width;
      std::cout<<higgs_width<<std::endl;
    }
  }
  return wzp;  
}

