#ifndef SPECTRUM_ANALYSIS_H
#define SPECTRUM_ANALYSIS_H


#include "xsec.h"
#include <vector>
#include <mutex>
#include <algorithm>
#include <thread>
#include <memory>
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
    std::mutex results_lock;
    std::mutex cout_mutex;
    std::shared_ptr< std::vector< std::pair<double,double> > > results;
    constexpr static double reldiff=1e2; //Maximum rel difference between two consecutive bins before switch
    T * pobj;
    double (T::* pfunc)(double, double, double, double(*)(double, double), int);
    double(* psmear)(double, double);
    void calculateBin(double lo, double hi, double acc);
  public:
    void writeHist(binning *pbins, double acc, char* outfile, double factor=1.);
    template<class Binning> void writeHist(double ll, double ul, double acc, char* outfile, double factor=1.);
    HistWriter(T* pobject, double (T::* pfunction)(double, double, double, double(*)(double, double), int), double(* psmearing)(double, double));
    HistWriter(pheno::HistWriter<T> &copy);
    ~HistWriter();
  };
  
    template<class T>
    pheno::HistWriter<T>::~HistWriter(){
      // Since all HistWriter classes share the same results pointer we have to
      // make sure, that it hasn't been deleted before
      // cout << "Destructor is called" << endl;
      // cout << results.use_count() << endl;
      // if(results.unique()){
      // 	cout <<"Delete results at "<< results <<endl;
      // 	//delete results.get();
      // 	cout << "Did it!"<< endl;
      // }
    }


    template<class T>
    pheno::HistWriter<T>::HistWriter(pheno::HistWriter<T> &copy){
      //cout << "Calling copy constructor " << endl;
      // All copies share the same results pointer to write to it
      //cout << "Use count before" << results.use_count() << endl;
      results = copy.results;
      //cout << "Use count after" << results.use_count() << endl;
      pobj=copy.pobj;
      pfunc=copy.pfunc;
      psmear=copy.psmear; 
    }
  
  
  
  template< class T>
  pheno::HistWriter<T>::HistWriter(T* pobject, double (T::* pfunction)(double, double, double, double(*)(double, double), int), double(* psmearing)(double, double))
  {
    std::vector<std::pair<double,double>> *res = new std::vector<std::pair<double,double>>();
    results = std::shared_ptr< std::vector<std::pair<double,double>> > (res);
    pobj=pobject;
    pfunc=pfunction;
    psmear=psmearing; 
  }
  
     
     
  
  template<class T>
  void pheno::HistWriter<T>::calculateBin(double lo, 
					  double hi, 
					  double acc
					    )
  {
    //cout_mutex.lock();
    //cout <<"Adress in calculateBin: " << this << endl;
    //cout_mutex.unlock();
    double res = (this->pobj->*(this->pfunc))(lo, hi, acc, this->psmear, 2);
    // Need guard to avoid deadlock when writing to results
    std::lock_guard<std::mutex> guard(results_lock);    
    
    //cout <<"Results adress in calculateBin: " << results << endl;
    //cout << results.get()->size() << endl;
    results.get()->push_back(std::pair<double,double>(lo,res));
    //cout << results.get()->size() << endl;
    //cout << "Use count" << results.use_count() << endl;
    //cout <<"Result saved" << endl;
  }
  
  
  template<class T>  
  void pheno::HistWriter<T>::writeHist(binning *pbins, double acc, char* outfile, double factor)
  {  
    double prev=0;
    int length = pbins->size();
    //int numCores = std::thread::hardware_concurrency();
    
    std::vector<std::thread> threads;
    std::vector<pheno::HistWriter<T>*> copies;
   
    //Write histogram to file in loop
    
    printf("Start integration...\n");
    //cout << "Adress before loop " << this << endl;
    //cout << "Result adress before loop" << &results << endl;
    //(this->results)->push_back(std::pair<double,double>(1.0,1.0));
    //cout << "Results size before loop: " << results.get()->size() << endl;

    for(int i=0; i<length; ++i)
    {      
      double low = (pbins->operator[](i)).first; //Get lower bound for bin
      double high = (pbins->operator[](i)).second; //Get upper bound for bin
      // // Create copy of HistWriter to avoid deadlock within thread
      //  pheno::HistWriter<T> *copy = new pheno::HistWriter<T>(*this);
      //  // Save pointers in a vector for deletion later on
      //  cout_mutex.lock();
      // cout << "Copy adress " <<copy << endl;
      // cout_mutex.unlock();
      // //cout << "This adress " << this << endl;
      // if(copy==this){
      // 	throw runtime_error("copy didn't allocate new adress");
      // }
      // for (auto it = copies.begin(); it != copies.end(); it++ ){
      // 	if(*it==copy) throw runtime_error("Adress of copy already exists!");
      // }
      // copies.push_back(copy);

      //calculateBin(low,high,acc);

      threads.push_back(
      			std::thread(
      				    &pheno::HistWriter<T>::calculateBin,
      				    this,
      				    //std::ref(copy),
      				    low, high, acc
      				    )
      			);
    }
	    
    // Wait for all threads to finish before writing the result file
    // cout << "##################################################" << endl;
    // cout << "Start joining" <<endl;
    for(auto th=threads.begin();th!=threads.end();th++) th->join();
	   
    
    cout << "##################################################" << endl;
    std::sort(
	      (*(results.get())).begin(),
	      (*(results.get())).end(),
	      [](std::pair<double,double> InA,std::pair<double,double> InB) {return InA.first < InB.first;} );
    std::printf("Writing data to %s ...\n", outfile); 
    std::ofstream outf(outfile);
    for(auto res=(*(results.get())).begin();res!=(*(results.get())).end();res++) {
      double 
	low = res->first,
	high = (res+1)->first;
      // Fix last bin
      if (res==((*results.get()).end()-1)) high=(pbins->operator[](length-1)).second;

      cout << low<<"\t" << high <<"\t" << (res->second)*factor << std::endl;
      outf << low<<"\t" << high <<"\t" << (res->second)*factor << std::endl;
    }
    outf.close();
    std::printf("Finished writing to %s\n", outfile);

	      //    for(auto cp=copies.begin();cp!=copies.end();cp++) delete *cp;
    
  }  
  
  
  
  
  
     
//OLD FUNCTIONS -- DON'T USE
//######################################################################################
     
     
          
  template<class T>
  double pheno::HistWriter<T>::writeHistCore(double lo, double hi, double acc, double prev)
  {
    //Integrate with Suave
    double res = (pobj->* pfunc)(lo, hi, acc, psmear, 2);
//     if(prev!=0)
//     {
//       double ratio1, ratio2, diff= res-prev;
//       if(diff>0)
//       {
//         ratio1 = diff/prev;
//         ratio2 = diff/res;
//       }else{
//         ratio1 = -diff/prev;
//         ratio2 = -diff/res;
//       }
// //       std::printf("Checking deviations %g %g\n", ratio1, ratio2);      //#############################################v DEBUG
//       //Check whether there is a hughe leap in the integral --> wrong convergence
//       if( ( ratio1 > reldiff) || ( ratio2 > reldiff) || ( -1*ratio1 > reldiff) || ( -1*ratio2 > reldiff) )
//       {
//         //Switch to Monte Carlo integration
//         std::printf("Match jumping criterion. Relative deviations are %g %g\n", ratio1, ratio2);      //#############################################v DEBUG
//         res = ((this->pobj)->* (this->pfunc))(lo, hi, 1e-2,  this->psmear, 2); 
//       }
//     }
    return res;
  }
  
  
  
  
  //Function that uses a binning functor to calculate bins
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

      double res=0;
      
      outf<<low<<"\t"<<res * factor <<"\n"; //Write the bin to file
      low = high; //Set new lower bound
      prev= res;  //Save the last value
      std::printf("Previous result %g\n", prev);      //#############################################v DEBUG
    }
    outf.close();
  }

  
  
  
 
 
}

#endif
