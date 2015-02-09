#include "pheno.h"
#include <iostream>

using namespace pheno;




//************************************//
//            QUARK CLASSES           //
//************************************//

///Implemenation of quark constructors

down::down(float xlc, float xrc, bool massive): fermionExt(massive, 1, -1./2, 4.8e-3, -1./3, xlc, xrc) { }

up::up(float xlc, float xrc, bool massive): fermionExt(massive, 1, 1./2, 2.3e-3, 2./3, xlc, xrc) { }

strange::strange(float xlc, float xrc, bool massive): fermionExt(massive, 2, -1./2, 9.5e-2, -1./3, xlc, xrc) { }

charm::charm(float xlc, float xrc, bool massive): fermionExt(massive, 2, 1./2, 1.275, 2./3, xlc, xrc) { }

bottom::bottom(float xlc, float xrc, bool massive): fermionExt(massive, 3, -1./2, 4.18, -1./3, xlc, xrc) { }

top::top(float xlc, float xrc, bool massive): fermionExt(massive, 3, 1./2, 173.5, 2./3, xlc, xrc) { }



//************************************//
//            LEPTON CLASSES          //
//************************************//


///Implemenation of lepton constructors

electron::electron(float xlc, float xrc, bool massive): fermionExt(massive, 1, -1./2, 5.11e-4, -1, xlc, xrc) { }

muon::muon(float xlc, float xrc, bool massive): fermionExt(massive, 2, -1./2, 0.1057, -1, xlc, xrc) { }

tauon::tauon(float xlc, float xrc, bool massive): fermionExt(massive, 3, -1./2, 1.777, -1, xlc, xrc) { }

nu_el::nu_el(float xlc, float xrc): fermionExt(false, 1, 1./2, 0, 0, xlc, xrc) { }

nu_mu::nu_mu(float xlc, float xrc): fermionExt(false, 2, 1./2, 0, 0, xlc, xrc) { }

nu_tau::nu_tau(float xlc, float xrc): fermionExt(false, 3, 1./2, 0, 0, xlc, xrc) { }




//************************************//
//            ZP-MODEL CLASS          //
//************************************//

//Constructor: The whole model is set up HERE
zpmodel::zpmodel(const char* configfile): bsm_parameters(0.1, 1500, 0)  //Default values if no parameters are specified in config file
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
      std::cout<<ferms->first<<":\t\tcxl="<<(ferms->second)->get_xlcharge()<<"\tcxr="<<(ferms->second)->get_xlcharge()<<"\tmass:"<<(ferms->second)->get_mass()<<"\tqgam:"<<((ferms->second)->vecc()).q_gam/e_()<<"\n";
    }
  }  
}

