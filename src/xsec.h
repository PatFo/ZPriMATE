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
      PartonXSec(fundamental::fermionExt f_in, fundamental::fermionExt f_out, pheno::zpmodel* p_model);
  };

}

#endif