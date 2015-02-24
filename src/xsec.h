#ifndef XSEC_H
#define XSEC_H

#include "pheno.h"
#include "../../../local/MSTW/mstwpdf.h"
//PDF package
#include <mstwpdf.h>


namespace pheno{
  
  
  //*******************************************************//
  //            Class for partonic cross sections          //
  //*******************************************************//

  
  class PartonXSec{
    ///Class for calculcation of partonic cross sections of f_in f_in~ --> f_out f_out~
    private:
      int pdgin;
      //Pointer to the model in use
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
      int pdg_in();
      //Partonic cross sections
      double sigZp(double Ecm);
      double sigSM(double Ecm);
      double sigInt(double Ecm);
      double sigTot(double Ecm);
      //Class Constructor: Give model as reference &model
      PartonXSec(fundamental::fermionExt* f_in, fundamental::fermionExt* f_out, pheno::zpmodel* p_model);
      ///Give the two fermions and the model as reference to the constructor!
  };

  
  
  
  
  
  
  
  //*********************************************************//
  //            Classes for hadronic cross sections          //
  //*********************************************************//

  

  
  class HadronXSec{
    ///Class for calculcation of partonic cross sections of p p --> f_out f_out~
    private:
      double Epp;
      //Internal parton cross sections: no top, as pdf negligible
      PartonXSec* dxsec;
      PartonXSec* uxsec;
      PartonXSec* sxsec;
      PartonXSec* cxsec;
      PartonXSec* bxsec;
      c_mstwpdf* pdf;
      //Subroutines for pdf convoluted cross sections: Specify partial cross section in functor object
      template<class PartialCrossX> double pdf_xsec(PartonXSec* pxsec, double q, double x);
      template<class PartialCrossX> double pdfconvoluted(PartonXSec* pxsec, double Ecm);
    public:
      //Hadronic cross sections
      double sigTot(double Ecm);
      //Constructor and destructor to take care of memory allocations
      HadronXSec(fundamental::fermionExt* f_out, pheno::zpmodel* p_model, char* pdf_grid_file, double Ecoll=8000.);
      ~HadronXSec();    
  };

  
  //Implementation of the member templates
  template<class PartialCrossX> 
  double pheno::HadronXSec::pdf_xsec(PartonXSec* pxsec, double q, double x)
  {
    PartialCrossX f;
    return f(pxsec, q) * pheno::HadronXSec::pdf->parton(pxsec->pdg_in(), x, q);
  }
  
  
  template<class PartialCrossX> 
  double pheno::HadronXSec::pdfconvoluted(PartonXSec* pxsec, double Ecm)
  {
    return 0;
  }
  
}

#endif