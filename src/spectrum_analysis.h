#ifndef SPECTRUM_ANALYSIS_H
#define SPECTRUM_ANALYSIS_H


#include "xsec.h"
#include <vector>
#include <fstream>
#include <stdio.h> //############################# DEBUG


namespace pheno {
  
  
  
  typedef std::vector<double*> sampling_scheme;
  typedef std::vector<std::pair<double, double> > binning;
  
    
  //*******************************************************//
  //            Class for scanning spectrum                //
  //*******************************************************//
  
  template<class CrossSection>
  class SpectrumScanner {
    private:
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
      SpectrumScanner(pheno::ZpModel* pmod, CrossSection* phsec);
      ~SpectrumScanner();
  };
  
  
  
  
  
  
  //*******************************************************//
  //            Function to get binning from file          //
  //*******************************************************//
  
  
  //Construct a histogram scheme from two column text file
  binning get_binning(char* binfile);
  
  
  
  
  //*******************************************************//
  //            Class for  histrogram output               //
  //*******************************************************//
  
  template< class T>
  class HistWriter {
    ///Class for writing histograms 
    ///Pass the binning scheme functor as temnplate argument 'Binning' and the class whose function gets plotted as argument 'T'
    ///Pass the function as a pointer
    private:
      const static double reldiff=1e2; //Maximum rel difference between two consecutive bins before switch
      T * pobj;
      double (T::* pfunc)(double, double, double, double(*)(double, double), int);
      double(* psmear)(double, double);
      double writeHistCore(double lo, double hi, double acc, double prev);
    public:
      void writeHist(binning *pbins, double acc, char* outfile, double factor=1.);
      template<class Binning> void writeHist(double ll, double ul, double acc, char* outfile, double factor=1.);
      HistWriter(T* pobject, double (T::* pfunction)(double, double, double, double(*)(double, double), int), double(* psmearing)(double, double));      
  };
  
  
  
  
  template< class T>
  pheno::HistWriter<T>::HistWriter(T* pobject, double (T::* pfunction)(double, double, double, double(*)(double, double), int), double(* psmearing)(double, double))
  {
    pheno::HistWriter<T>::pobj=pobject;
    pheno::HistWriter<T>::pfunc=pfunction;
    pheno::HistWriter<T>::psmear=psmearing; 
  }
  
     
     
     
  template<class T>
  double pheno::HistWriter<T>::writeHistCore(double lo, double hi, double acc, double prev)
  {
    //Integrate with cubature
    double res = ((this->pobj)->* (this->pfunc))(lo, hi, acc, this->psmear, 2);
    if(prev!=0)
    {
      double ratio1, ratio2, diff= res-prev;
      if(diff>0)
      {
        ratio1 = diff/prev;
        ratio2 = diff/res;
      }else{
        ratio1 = -diff/prev;
        ratio2 = -diff/res;
      }
//       std::printf("Checking deviations %g %g\n", ratio1, ratio2);      //#############################################v DEBUG
      //Check whether there is a hughe leap in the integral --> wrong convergence
      if( ( ratio1 > (this->reldiff)) || ( ratio2 > (this->reldiff)) )
      {
        //Switch to Monte Carlo integration
        std::printf("Match jumping criterion. Relative deviations are %g %g\n", ratio1, ratio2);      //#############################################v DEBUG
        res = ((this->pobj)->* (this->pfunc))(lo, hi, 1e-2,  this->psmear, 1); 
      }
    }
    return res;
  }

  
  
  
  template<class T>  
  void pheno::HistWriter<T>::writeHist(binning *pbins, double acc, char* outfile, double factor)
  {  
    std::ofstream outf(outfile);
    double prev=0;
    int length = pbins->size();
    
    std::printf("Writing data to %s ...\n", outfile); 
    //Write histogram to file in loop
    for(int i=0; i<length; ++i)
    {
      double accuracy = acc;
      //Calculate first point with high accuracy to have a reliable starting point for convergence test
      if(i==0){accuracy= 1e-6;} 
      
      double low = (pbins->operator[](i)).first; //Get lower bound for bin
      double high = (pbins->operator[](i)).second; //Get upper bound for bin

      //Call core integration function; teh cross section 'res' is given in [fb]
      double res = this->writeHistCore(low, high, accuracy, prev);
      
      outf<<low<<"\t"<<high<<"\t"<<res * factor <<"\n"; //Write the bin to file
      prev= res;  //Save the last value
//       std::printf("Previous result %g\n", prev);      //#############################################v DEBUG
    }
    outf.close();
    std::printf("Finished writing %s\n", outfile);      //#############################################v DEBUG
  }  
     
  
  template<class T>
  template<class Binning>
  void pheno::HistWriter<T>::writeHist(double ll, double ul, double acc, char* outfile, double factor)
  {  
    //Create instance of binning functor
    Binning f;
    std::ofstream outf(outfile);
    double prev=0;
    //Write histogram to file in loop
    for(double low = ll; low<ul; )
    {
      double high = f(low); //Calculate upper bound for bin

      double res = this->writeHistCore(low, high, acc, prev);
      
      outf<<low<<"\t"<<res * factor <<"\n"; //Write the bin to file
      low = high; //Set new lower bound
      prev= res;  //Save the last value
      std::printf("Previous result %g\n", prev);      //#############################################v DEBUG
    }
    outf.close();
  }

  
  
  
 
 
}

#endif