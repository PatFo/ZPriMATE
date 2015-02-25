#ifndef SPECTRUM_ANALYSIS_H
#define SPECTRUM_ANALYSIS_H


#include "xsec.h"
#include <vector>


namespace pheno {
  
  
  //*******************************************************//
  //            Class for scanning spectrum                //
  //*******************************************************//
  
  class SpectrumScanner {
    private:
      zpmodel* _model;
      HadronXSec* _hsec;
      //Sampling Regions
      bool is_default;
      std::vector<double*> samplingRegions;
      void reset_regions();
      //Sampler function
      void sampler(char* outfile, double low, double high, double step);
    public:
      void set_interval(double low, double high, double step);
      void scan(char* outfile);
      SpectrumScanner(pheno::zpmodel* pmod, pheno::HadronXSec* phsec);
      ~SpectrumScanner();
    
  };
  
  
}



#endif