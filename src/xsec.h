#ifndef XSEC_H
#define XSEC_H

#include "pheno.h"


namespace pheno{
  
  
  //*******************************************************//
  //            Class for partonic cross sections          //
  //*******************************************************//

  
  class PartonXSec{
    ///Class for calculcation of partonic cross sections of f_in f_in~ --> f_out f_out~
    private:
      //
      pheno::zpmodel *  _model;      
      //Numerical factors for partial cross sections
      double numGam;
      double numZ;
      double numZp;
      double numGamZ;
      double numGamZp;
      double numZZp;
      //Partial cross sections (Zp part can be called->public)
      double sigGam(double Ecm);
      double sigZ(double Ecm);
      double sigGamZ(double Ecm);
      double sigGamZp(double Ecm);
      double sigZZp(double Ecm);
    public:
      //Partonic cross sections
      double sigZp(double Ecm);
      double sigSM(double Ecm);
      double sigInt(double Ecm);
      double sigTot(double Ecm);
      //Class Constructor: Give model as reference &model
      PartonXSec(fundamental::fermionExt* f_in, fundamental::fermionExt* f_out, pheno::zpmodel* p_model);
      ///Give the two fermions and the model as reference to the constructor!
  };

  
  
  
  
  
  
  
  //*******************************************************//
  //            Class for hadronic cross sections          //
  //*******************************************************//

  
  class HadronXSec{
    ///Class for calculcation of partonic cross sections of p p --> f_out f_out~
    private:
      //Internal parton cross sections: no top, as pdf negligible
      PartonXSec* dxsec;
      PartonXSec* uxsec;
      PartonXSec* sxsec;
      PartonXSec* cxsec;
      PartonXSec* bxsec;
    public:
      //Hadronic cross sections
      double sigTot(double Ecm);
      //Constructor and destructor to take care of memory allocations
      HadronXSec(fundamental::fermionExt* f_out, pheno::zpmodel* p_model);
      ~HadronXSec();
    
    
  };

  
}

#endif