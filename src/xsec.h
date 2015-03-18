#ifndef XSEC_H
#define XSEC_H




#include "pheno.h"
#include <stdexcept>
#include <vector>
//PDF package
#include <mstwpdf.h>
//Numerical integration
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>


namespace pheno{
  
  
  const double GeV2fb =  0.3894e12;
  
  //*******************************************************//
  //            Class for partonic cross sections          //
  //*******************************************************//

  
  class PartonXSec{
    ///Class for calculcation of partonic cross sections of f_in f_in~ --> f_out f_out~
    private:
      int pdgin;
      //Pointer to the model in use
      pheno::ZpModel *  _model;      
      //Numerical factors for partial cross sections
      double numGam;
      double numZ;
      double numZp;
      double numGamZ;
      double numGamZp;
      double numZZp;
      //Partial cross sections (Zp part can be called->public)
    public: // ##############################     DEBUG  ####################3
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
      //Member function that fills EMPTY vector with sigSM, sigInt, sigSignal, sigTotal 
      //ALWAYS use this function if more than one of these cross sections is needed at a time
      //Strategy only needed for interface with SpectrumScanner -> NOT USED
      void crossSections (double Ecm, std::vector<double> * results, unsigned int int_strategy=1);
      //Class Constructor: Give model as reference &model
      PartonXSec(fundamental::fermionExt* f_in, fundamental::fermionExt* f_out, pheno::ZpModel* p_model);
      ///Give the two fermions and the model as reference to the constructor!
  };

  
  
  
  
  
  
  
  //*********************************************************//
  //            Classes for hadronic cross sections          //
  //*********************************************************//

  

  
  class HadronXSec{
    ///Class for calculcation of partonic cross sections of p p --> f_out f_out~
    private:
      size_t calls;
      double accuracy_goal;
      double Epp;
      //Internal parton cross sections: no top, as pdf negligible
  public: //####################################################################################### DEBUG ONLY ####################################
      PartonXSec* dxsec;
      PartonXSec* uxsec;
      PartonXSec* sxsec;
      PartonXSec* cxsec;
      PartonXSec* bxsec;
      c_mstwpdf* pdf;
      //Subroutine for pdf convoluted cross sections: Specify partial cross section in functor object
      template<class PartialCrossX> double pdfconvoluted( double Ecm, unsigned int int_strategy=1);
      //Calculate cross section bin bin between [el, eh]
      template<class PartialCrossX>  friend  double dSigdM(double E, void * p);
      template<class PartialCrossX> double binnedXsec(double el, double eh, double accuracy);
    public:
      void set_accuracy(double accuracy);
      void set_monte_calls(size_t int_calls);
      //Hadronic cross sections
      double sigSM    (double Ecm, unsigned int int_strategy=1);
      double sigInt   (double Ecm, unsigned int int_strategy=1);
      double sigSignal(double Ecm, unsigned int int_strategy=1);
      double sigTotal (double Ecm, unsigned int int_strategy=1);
      //Member function that fills EMPTY vector with sigSM, sigInt, sigSignal, sigTotal 
      //ALWAYS use this function if more than one of these cross sections is needed at a time
      void crossSections (double Ecm, std::vector<double> * results, unsigned int int_strategy=1);
      double totXsec(double el, double eh, double accuracy);
      double zpXsec(double el, double eh, double accuracy);
      //Constructor and destructor to take care of memory allocations
      HadronXSec(fundamental::fermionExt* f_out, pheno::ZpModel* p_model, char* pdf_grid_file, double Ecoll=8000.);
      ~HadronXSec();  
  };
  
  
  
  
  
  
  //IMPLEMENTATION OF MEMBER TEMPLATE  
  //---------------------------------------------------------------------------------------------------------
  
  
  

  //Parameter struct for integrable function
  struct parameters{ PartonXSec** ppx; int arr_size; c_mstwpdf* ppdf; double Ecm; double Ecoll; double* cross_sections; };
  
  
  
  //Implementation of gsl monte carlo integrable function: 
  //Numerical Integration
  inline double num_pdf(double x, void * p)
  {
    //Get initialization parameters
    struct parameters * pars = (struct parameters *)p;
    //Define Bjorken x
    double x1 = x;
    double x2 = pars->Ecm * pars->Ecm/(x1 * pars->Ecoll * pars->Ecoll);
    //Calculate sum of cross sections
    double sum = 0.;
    for(int i=0; i<pars->arr_size; ++i)
    {
      double pdg = (pars->ppx[i])->pdg_in();
      sum +=  pars->cross_sections[i]/(pars->Ecoll * pars->Ecoll) *  //Parton level cross section
              ( 
              pars->ppdf->parton( -1*pdg,  x1, pars->Ecm )/(x1*x1) * //PDFs with antiquark in first proton
              pars->ppdf->parton(    pdg,  x2, pars->Ecm )/x2
              + 
              pars->ppdf->parton(     pdg,  x1, pars->Ecm )/(x1*x1) * //Mirror: PDFs with antiquark in second proton   
              pars->ppdf->parton(  -1*pdg,  x2, pars->Ecm )/x2 
              );
    }
    return sum;
  }
  
  
  
  //Implementation of gsl monte carlo integrable function: 
  //Monte Carlo Integration
  inline double monte_pdf(double x[], size_t dim, void * p)
  {
    //Check that integration vector is one dimensional
    if(dim!=1)
      throw std::runtime_error("ERROR: Integration dim!=1");
    //Use implementation for regular numerical integration
    return num_pdf(x[0], p);
  }
  
  
  
  
  //Function that integrates cross section 
  template<class PartialCrossX> 
  double pheno::HadronXSec::pdfconvoluted( double Ecm, unsigned int int_strategy)
  {
    //Construct parameters
    PartialCrossX f;
    PartonXSec* pp[]= {dxsec, uxsec, sxsec, cxsec, bxsec};
    //Store partonic cross section of quarks at energy Ecm in array
    double cxs[] = {f(dxsec, Ecm),
                    f(uxsec, Ecm), 
                    f(sxsec, Ecm),
                    f(cxsec, Ecm),
                    f(bxsec, Ecm)};
                    
    struct parameters local_pars = {pp, 5, pdf, Ecm, Epp, cxs};
    
    //Parameter and result objects
    double result, error;
    double xl = Ecm*Ecm/(Epp*Epp);  //lower integration limit
    double xu = 1;  //Upper integration limit
    
    if(int_strategy==1)
    {
      // QAG adaptive integration
      //Define Integration function
      gsl_function F;
      F.function = &num_pdf;
      F.params = &local_pars;
      
      //Integration
      size_t bisection_lim =1000;
      gsl_integration_workspace * w = gsl_integration_workspace_alloc (bisection_lim);
      gsl_integration_qags (&F, xl, xu, 0, accuracy_goal, bisection_lim, w, &result, &error); 
      //Free memory of integration workspace
      gsl_integration_workspace_free(w);
    }
    else if(int_strategy==2)
    {
      // VEGAS Monte Carlo Integration
      //Define a gsl_monte_function to pass to the integrator
      gsl_monte_function F;
      F.f = &monte_pdf;
      F.dim = 1;
      F.params = &local_pars;      

      //Initialize random number gen for monte carlo
      gsl_rng_env_setup ();
      const gsl_rng_type *T;
      gsl_rng *r;
      T = gsl_rng_default;
      r = gsl_rng_alloc (T);

      //Integration
      double axl[1] = {xl};
      double axu[1] = {xu};
      gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);
      gsl_monte_vegas_integrate (&F, axl, axu, 1, calls, r, s, &result, &error);
      //Free integration memory
      gsl_monte_vegas_free(s);
    }
    else
    {
      throw std::runtime_error("ERROR: Undefined integration strategy code.");
    }
    
    return result;
  }
  
  
  
  
  template<class PartialCrossX>   
  double dSigdM(double E, void * p)
  {
    HadronXSec * phsec = (HadronXSec *) p;
    return 2*E * phsec->pdfconvoluted<PartialCrossX>(E, 1); 
  }
  
  
  
  template<class PartialCrossX> 
  double pheno::HadronXSec::binnedXsec(double el, double eh, double accuracy)
  {
    //Implement binwise integration scheme
    //Parameter and result objects
    double result, error;
    
    // QAG adaptive integration
    //Define Integration function
    gsl_function F;
    F.function = &dSigdM<PartialCrossX>;
    F.params = this;  
    
    //Integration
    size_t bisection_lim =1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (bisection_lim);
    gsl_integration_qags (&F, el, eh, 0, accuracy, bisection_lim, w, &result, &error); 
    //Free memory of integration workspace
    gsl_integration_workspace_free(w);
    
    return result;
  }

  
  
  
}

#endif