#ifndef SPECTRUM_ANALYSIS_H
#define SPECTRUM_ANALYSIS_H


#include "xsec.h"
#include <vector>
#include <fstream>
#include <stdio.h> //############################# DEBUG


namespace pheno {
  
  
  
  typedef std::vector<double*> sampling_scheme;
  
  //*******************************************************//
  //            Class for scanning spectrum                //
  //*******************************************************//
  
  template<class CrossSection>
  class SpectrumScanner {
    private:
      unsigned int strategy;
      ZpModel* _model;
      CrossSection* _hsec;
      //Sampling Regions
      bool is_default;
      sampling_scheme samplingRegions;
      //Sampler function
      void sampler(char* outfile, double low, double high, double step);
    public:
      void add_interval(double low, double high, double step);
      void add_interval(double* pinterval); //WARNING: Needs a static array of type double with 3 elements, e.g. double arr[3]={1.,2.,3.};
      void reset_sampling();
      void scan(char* outfile);
      SpectrumScanner(pheno::ZpModel* pmod, CrossSection* phsec, unsigned int int_strategy=1);
      ~SpectrumScanner();
  };
  
  
  
  
  //*******************************************************//
  //            Class for  histrogram output               //
  //*******************************************************//
  
  template<class Binning, class T>
  class HistWriter {
    ///Class for writing histograms 
    ///Pass the binning scheme functor as temnplate argument 'Binning' and the class whose function gets plotted as argument 'T'
    ///Pass the function as a pointer
    private:
      const static double reldiff=1e5;
      T * pobj;
      double (T::* pfunc)(double, double, double, double(*)(double, double), int);
      double(* psmear)(double, double);
    public:
      void writeHist(double ll, double ul, double acc, char* outfile, double factor=1.);
      HistWriter(T* pobject, double (T::* pfunction)(double, double, double, double(*)(double, double), int), double(* psmearing)(double, double));      
  };
  
  
  
  
  template<class Binning, class T>
  pheno::HistWriter<Binning, T>::HistWriter(T* pobject, double (T::* pfunction)(double, double, double, double(*)(double, double), int), double(* psmearing)(double, double))
  {
    pheno::HistWriter<Binning, T>::pobj=pobject;
    pheno::HistWriter<Binning, T>::pfunc=pfunction;
    pheno::HistWriter<Binning, T>::psmear=psmearing; 
  }
  
  
  
  template<class Binning, class T>
  void pheno::HistWriter<Binning, T>::writeHist(double ll, double ul, double acc, char* outfile, double factor)
  {  
    //Create instance of binning functor
    Binning f;
    std::ofstream outf(outfile);
    double prev=0;
    //Write histogram to file in loop
    for(double low = ll; low<ul; )
    {
      double high = f(low); //Calculate upper bound for bin
      double res = ((this->pobj)->* (this->pfunc))(low, high, acc, this->psmear, 2);
      if(prev!=0)
      {
        double ratio1 = abs ((res-prev)/prev);
        double ratio2 = abs((res-prev)/res);
        std::printf("Checking deviations %g %g\n", ratio1, ratio2);      //#############################################v DEBUG
        //Check whether there is a hughe leap in the integral --> wrong convergence
        if( ( ratio1 > (this->reldiff)) || (ratio2 > (this->reldiff)) )
        {
          std::printf("Match jumping criterion. Relative deviations are %g %g\n", ratio1, ratio2);      //#############################################v DEBUG
          res = ((this->pobj)->* (this->pfunc))(low, high, 1e-2,  this->psmear, 1); 
        }
      }
      outf<<low<<"\t"<<res * factor <<"\n"; //Write the bin to file
      low = high; //Set new lower bound
      prev= res;  //Save the last value
      std::printf("Previous result %g\n", prev);      //#############################################v DEBUG
    }
    outf.close();
  }

  
  
  
 
 
}

#endif