#ifndef PHENO_H
#define PHENO_H


#include "model.h"
#include "read_config.h"

namespace pheno{

  
  //************************************//
  //            QUARK CLASSES           //
  //************************************//
    
  class down : public fundamental::fermionExt{
  /// Down-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      down(float xlc=1., float xrc=1., bool massive = false);
  };
  
  class up : public fundamental::fermionExt{
  /// Up-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      up(float xlc=1., float xrc=1., bool massive = false);
  };
  
  class strange : public fundamental::fermionExt{
  /// Strange-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      strange(float xlc=1., float xrc=1., bool massive = false);
  };
  
  class charm : public fundamental::fermionExt{
  /// Charm-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      charm(float xlc=1., float xrc=1., bool massive = false);
  };
  
  class bottom : public fundamental::fermionExt{
  /// Bottom-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      bottom(float xlc=1., float xrc=1., bool massive = false);
  };
  
  class top : public fundamental::fermionExt{
  /// Top-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      top(float xlc=1., float xrc=1., bool massive = false);
  };


  
  //************************************//
  //            LEPTON CLASSES          //
  //************************************//
  
  //Charged leptons
  class electron : public fundamental::fermionExt{
  /// Electron class that has default U(1)_X charges xlc=1, xrc=1
    public:
      electron(float xlc=1., float xrc=1., bool massive = false);
  };
  
  class muon : public fundamental::fermionExt{
  /// Muon class that has default U(1)_X charges xlc=1, xrc=1
    public:
      muon(float xlc=1., float xrc=1., bool massive = false);
  };
  
  class tauon : public fundamental::fermionExt{
  /// Tauon class that has default U(1)_X charges xlc=1, xrc=1
    public:
      tauon(float xlc=1., float xrc=1., bool massive = false);
  };
  
  //Neutrinos
  class nu_el : public fundamental::fermionExt{
  /// Electron neutrino class that has default U(1)_X charges xlc=1, xrc=1
    public:
      nu_el(float xlc=1., float xrc=1.);
  };
  
  class nu_mu : public fundamental::fermionExt{
  /// Muon neutrino class that has default U(1)_X charges xlc=1, xrc=1
    public:
      nu_mu(float xlc=1., float xrc=1.);
  };
  
  class nu_tau : public fundamental::fermionExt{
  /// Tauon neutrino class that has default U(1)_X charges xlc=1, xrc=1
    public:
      nu_tau(float xlc=1., float xrc=1.);
  };
  
  
  
  
  //************************************//
  //            ZP-MODEL CLASS          //
  //************************************//
  
  class zpmodel: public fundamental::bsm_parameters{
   /// Zp model class consisiting of all parameters, couplings and fermions
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
      tauon 	ta;
      nu_el 	ne;
      nu_mu 	nm;
      nu_tau	nt;
      
      //Constructor
      zpmodel(const char* configfile);
  };
  
  
  
}


#endif