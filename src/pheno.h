#ifndef PHENO_H
#define PHENO_H


#include "model.h"

namespace pheno{

  
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
  
}


#endif