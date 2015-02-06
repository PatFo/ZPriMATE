#ifndef PHENO_H
#define PHENO_H


#include "model.h"

namespace pheno{

  
  class up : public fundamental::fermion{
  /// Up-quark class that has default U(1)_X charges xlc=1, xrc=1
    public:
      up(float xlc=1., float xrc=1.);
  };

  
}


#endif