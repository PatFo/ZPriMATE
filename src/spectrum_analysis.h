#ifndef SPECTRUM_ANALYSIS_H
#define SPECTRUM_ANALYSIS_H


#include "xsec.h"
#include <vector>
#include <fstream>
#include <stdio.h> //############################# DEBUG


namespace pheno {
    

  inline void printProgBar( int percent ){
    std::string bar;
    int barSize = 40;
    double dpercent = ((double)percent)/100;
    for(int i = 0; i < barSize; i++){
      if( i < (int)(dpercent*barSize)){
	bar.replace(i,1,"=");
      }else if( i == (int)(dpercent*barSize)){
	bar.replace(i,1,">");
      }else{
	bar.replace(i,1," ");
      }
    }
    
    std::cout<< "\r" "[" << bar << "] ";
    std::cout.width( 3 );
    std::cout<< percent << "%" << std::flush;
  } 
  
  
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
      constexpr static double reldiff=1e2; //Maximum rel difference between two consecutive bins before switch
      T * pobj;
      double (T::* pfunc)(double, double, double, double(*)(double, double), int);
      double(* psmear)(double, double);
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
  void pheno::HistWriter<T>::writeHist(binning *pbins, double acc, char* outfile, double factor)
  {  
    double prev=0;
    int length = pbins->size();
    binning bincopy(*pbins);
    //Create a container of binning size to store results
    std::vector<double> prediction(length);
    double progress=0;
    printProgBar(0);
    for(int i1=0; i1<length; ++i1)
      {

        //Allocating private copies of all pointers in the game to be used by each thread
        
        //Get lower and upper bound of bin
        double low = bincopy[i1].first; 
        double high = bincopy[i1].second;
        
        //Call core integration function
        //The cross section is given in [fb]; multiply by factor=luminosity to obtain events
        prediction[i1] = (pobj->* pfunc)(low, high, acc, psmear, 2) * factor;
	progress = ((double)i1+1)/((double)length);
	printProgBar(progress*100);
      }
    // Newline for progressbar
    std::cout << std::endl << std::endl;

    //Write the prediction to file
    std::fprintf(stderr,"Writing data to %s ...\n", outfile); 
    std::ofstream outf(outfile);

    for(int i2=0; i2<length; ++i2)
    { 
      outf<<(pbins->operator[](i2)).first<<"\t"<<(pbins->operator[](i2)).second<<"\t"<<prediction[i2]<<"\n"; 
    }
    outf.close();
    std::fprintf(stderr,"Finished writing to %s\n", outfile);      //#############################################v DEBUG
  }    
}

#endif
