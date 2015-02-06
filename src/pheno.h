#ifndef PHENO_H
#define PHENO_H


#include "model.h"

namespace pheno{

  
  //************************************//
  //            QUARK CLASSES           //
  //************************************//
    
  class up : public fundamental::fermion{
  /// Up-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      up(float xlc=1., float xrc=1.);
  };
  
  class down : public fundamental::fermion{
  /// Down-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      down(float xlc=1., float xrc=1.);
  };
  
  class charm : public fundamental::fermion{
  /// Charm-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      charm(float xlc=1., float xrc=1.);
  };
  
  class strange : public fundamental::fermion{
  /// Strange-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      strange(float xlc=1., float xrc=1.);
  };
  
  class top : public fundamental::fermion{
  /// Top-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      top(float xlc=1., float xrc=1.);
  };

  class bottom : public fundamental::fermion{
  /// Bottom-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      bottom(float xlc=1., float xrc=1.);
  };

  
  //************************************//
  //            LEPTON CLASSES          //
  //************************************//
  
  //Charged leptons
  class electron : public fundamental::fermion{
  /// Electron class that has default U(1)_X charges xlc=1, xrc=1
    public:
      electron(bool massless = true, float xlc=1., float xrc=1.);
  };
  
  class muon : public fundamental::fermion{
  /// Muon class that has default U(1)_X charges xlc=1, xrc=1
    public:
      muon(bool massless = true, float xlc=1., float xrc=1.);
  };
  
  class tauon : public fundamental::fermion{
  /// Tauon class that has default U(1)_X charges xlc=1, xrc=1
    public:
      tauon(bool massless = true, float xlc=1., float xrc=1.);
  };
  
  //Neutrinos
  class nu_el : public fundamental::fermion{
  /// Electron neutrino class that has default U(1)_X charges xlc=1, xrc=1
    public:
      nu_el(float xlc=1., float xrc=1.);
  };
  
  class nu_mu : public fundamental::fermion{
  /// Muon neutrino class that has default U(1)_X charges xlc=1, xrc=1
    public:
      nu_mu(float xlc=1., float xrc=1.);
  };
  
  class nu_tau : public fundamental::fermion{
  /// Tauon neutrino class that has default U(1)_X charges xlc=1, xrc=1
    public:
      nu_tau(float xlc=1., float xrc=1.);
  };
  
  
}


#endif