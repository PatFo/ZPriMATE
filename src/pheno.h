#ifndef PHENO_H
#define PHENO_H

#include <map>
#include <string>
#include "model.h"
#include "read_config.h"

namespace pheno{

  
  //************************************//
  //            QUARK CLASSES           //
  //************************************//
    
  class down : public fundamental::fermionExt{
  /// Down-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      down(double xlc=1., double xrc=1., bool massive = true);
  };
  
  class up : public fundamental::fermionExt{
  /// Up-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      up(double xlc=1., double xrc=1., bool massive = true);
  };
  
  class strange : public fundamental::fermionExt{
  /// Strange-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      strange(double xlc=1., double xrc=1., bool massive = true);
  };
  
  class charm : public fundamental::fermionExt{
  /// Charm-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      charm(double xlc=1., double xrc=1., bool massive = true);
  };
  
  class bottom : public fundamental::fermionExt{
  /// Bottom-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      bottom(double xlc=1., double xrc=1., bool massive = true);
  };
  
  class top : public fundamental::fermionExt{
  /// Top-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      top(double xlc=1., double xrc=1., bool massive = true);
  };


  
  //************************************//
  //            LEPTON CLASSES          //
  //************************************//
  
  //Charged leptons
  class electron : public fundamental::fermionExt{
  /// Electron class that has default U(1)_X charges xlc=1, xrc=1
    public:
      electron(double xlc=1., double xrc=1., bool massive = true);
  };
  
  class muon : public fundamental::fermionExt{
  /// Muon class that has default U(1)_X charges xlc=1, xrc=1
    public:
      muon(double xlc=1., double xrc=1., bool massive = true);
  };
  
  class tauon : public fundamental::fermionExt{
  /// Tauon class that has default U(1)_X charges xlc=1, xrc=1
    public:
      tauon(double xlc=1., double xrc=1., bool massive = true);
  };
  
  //Neutrinos
  class nu_el : public fundamental::fermionExt{
  /// Electron neutrino class that has default U(1)_X charges xlc=1, xrc=0 (no right-handed neutrino in SM)
    public:
      nu_el(double xlc=1., double xrc=0);
  };
  
  class nu_mu : public fundamental::fermionExt{
  /// Muon neutrino class that has default U(1)_X charges xlc=1, xrc=0 (no right-handed neutrino in SM)
    public:
      nu_mu(double xlc=1., double xrc=0);
  };
  
  class nu_tau : public fundamental::fermionExt{
  /// Tauon neutrino class that has default U(1)_X charges xlc=1, xrc=0 (no right-handed neutrino in SM)
    public:
      nu_tau(double xlc=1., double xrc=0);
  };
  
  
  //************************************//
  //            HELPER TYPE             //
  //************************************//
  
  
  typedef std::map<std::string, fundamental::fermionExt*> fermion_list;

  
  
  //************************************//
  //            ZP-MODEL CLASS          //
  //************************************//
  
   /// Zp model class consisiting of all parameters, couplings and fermions
  class ZpModel: public fundamental::bsm_parameters{
    public:
      //Quarks
      down 	d;
      up 	u;
      strange 	s;
      charm 	c;
      bottom 	b;
      top 	t;
      //Leptons
      electron 	el;
      muon 	mu;
      tauon 	tau;
      nu_el 	ne;
      nu_mu 	nm;
      nu_tau	nt;
    private:
      fermion_list flst;
      void setup_flst();
      //Zp width related
      double* partial_fwidths;
      double higgs_width;
      double ww_width;
      double wzp;
      double calc_fwidth(fundamental::fermionExt &f);    
      double calc_zhwidth();
      double calc_wwidth();
      void update_width();
    public:
      //Constructors
      ZpModel(const char* configfile);
      ZpModel();
      //Width of Zp
      double wzp_();
  };
  
  
  
}


#endif