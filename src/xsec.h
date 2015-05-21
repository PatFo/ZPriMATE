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
//Cubature integration methods
#include <cubature.h>
#include <cuba.h>


#include <iostream> //DEBUG

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
      void crossSections (double Ecm, std::vector<double> * results);
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
      template<class PartialCrossX> double pdfconvoluted( double Ecm );
      //Full experimentally detectable cross section
      template<class PartialCrossX> double theoXsec(double el, double eh, double acc);
      template<class PartialCrossX> double detectedXsec(double el, double eh, double acc, int strategy, double (* psmear)(double, double));
    public:
      void set_accuracy(double accuracy);
      void set_monte_calls(size_t int_calls);
      //Hadronic differential cross sections dSig/dm
      double dsigSM    (double Ecm );
      double dsigInt   (double Ecm );
      double dsigSignal(double Ecm );
      double dsigTotal (double Ecm );
      //Member function that fills EMPTY vector with sigSM, sigInt, sigSignal, sigTotal 
      //ALWAYS use this function if more than one of these cross sections is needed at a time
      void crossSections (double Ecm, std::vector<double> * results);
      //Full pure hadronic cross sections
      double smXsec(double el, double eh, double accuracy, double (* psmear)(double,double)=NULL, int strategy=1);
      double zpXsec(double el, double eh, double accuracy, double (* psmear)(double,double)=NULL, int strategy=1);
      double totXsec(double el, double eh, double accuracy, double (* psmear)(double,double)=NULL, int strategy=1);
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
              pars->ppdf->parton( -1*pdg,  x1, pars->Ecoll )/(x1*x1) * //PDFs with antiquark in first proton
              pars->ppdf->parton(    pdg,  x2, pars->Ecoll )/x2
              + 
              pars->ppdf->parton(     pdg,  x1, pars->Ecoll )/(x1*x1) * //Mirror: PDFs with antiquark in second proton   
              pars->ppdf->parton(  -1*pdg,  x2, pars->Ecoll )/x2 
              );
    }
    return sum;
  }
  
  
  

  
  
  //Function that integrates cross section 
  template<class PartialCrossX> 
  double pheno::HadronXSec::pdfconvoluted( double Ecm )
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
    
    return result;
  }
  
  


  
  //Total binned cross section integration (theoretical) 
  //-------------------------------------------------------------------------------------------
  
  //Parameter struct for integrable function
  struct parameter_set{ int narr; PartonXSec** ppx; c_mstwpdf* ppdf; double (* psmear)(double, double); double Ecoll; };
  
  
  //Parameter struct for integrable function
  struct parameter_set_hahn{ parameter_set * pars; double*  low; double*  high;};
    

  
  
  
    //Cubature adaptor (Johnson)
  template<class PartialCrossX> 
  inline int integrand_th(unsigned ndim, const double *d, void *fdata, unsigned fdim, double *fval)
  {
    //Get initialization parameters
    PartialCrossX cross;
    struct parameter_set* pars = (struct parameter_set *)fdata;
    
    //Define Bjorken x
    double scoll = pars->Ecoll * pars->Ecoll;
    double dx1dy = (scoll - d[0]*d[0])/scoll;
    double x1 = dx1dy * d[1] + d[0]*d[0]/scoll; // Variable trafo from x1 to y=a*x1 + b
    double x2 = d[0]*d[0]/(x1 * scoll);
    
    //Calculate sum of cross sections
    fval[0]= 0.;
    for(int i=0; i<pars->narr; ++i)
    {
      double pdg = (pars->ppx[i])->pdg_in();
      fval[0] +=  dx1dy * 2 * d[0]/scoll *     //Prefactors 
              cross(pars->ppx[i], d[0]) *  //Parton level cross section
              ( 
              pars->ppdf->parton( -1*pdg,  x1, pars->Ecoll )/(x1*x1) * //PDFs with antiquark in first proton
              pars->ppdf->parton(    pdg,  x2, pars->Ecoll )/x2
              + 
              pars->ppdf->parton(    pdg,  x1, pars->Ecoll )/(x1*x1) * //Mirror: PDFs with antiquark in second proton   
              pars->ppdf->parton( -1*pdg,  x2, pars->Ecoll )/x2 );
    }    
    return 0;
  }
  
  
  
  
  
  template<class PartialCrossX> 
  double pheno::HadronXSec::theoXsec(double el, double eh, double acc)
  {    
    //Construct parameters
    PartonXSec* pp[]= {dxsec, uxsec, sxsec, cxsec, bxsec};                    
    struct parameter_set int_pars = {5, pp, pdf, NULL, Epp};  
    double epsabs(0.), epsrel(acc);
    const int dimint(2), dimres(1);

    //Integration variables
    double result, error;
    double xl[2] = {el, 0};
    double xu[2] = {eh, 1};    
 
    //Integration
    hcubature(dimres, &integrand_th<PartialCrossX>, &int_pars, dimint, xl, xu, 0, epsabs, epsrel, ERROR_INDIVIDUAL, &result, &error);
//     std::cout<<"Integral "<<result<<" Error: "<<error<<std::endl; //######################## DEBUG
    
    return result;
  }
  
  
  
  
  
  //Total binned cross section integration (detectable/smeared) 
  //-------------------------------------------------------------------------------------------
  
  
  
  
  
  //Implementation of gsl monte carlo integrable function: 
  //Numerical Integration
  template<class PartialCrossX> 
  inline double integrand(double d[], size_t dim, void * p)
  {
    //Get initialization parameters
    PartialCrossX cross;
    struct parameter_set* pars = (struct parameter_set *)p;
    
    //Define Bjorken x
    double scoll = pars->Ecoll * pars->Ecoll;
    double dx1dy = (scoll - d[1]*d[1])/scoll;
    double x1 = dx1dy * d[2] + d[1]*d[1]/scoll; // Variable trafo from x1 to y=a*x1 + b
    double x2 = d[1]*d[1]/(x1 * scoll);
    
    //Calculate sum of cross sections
    double sum = 0.;
    for(int i=0; i<pars->narr; ++i)
    {
      double pdg = (pars->ppx[i])->pdg_in();
      sum +=  dx1dy * 2 * d[1]/scoll *     //Prefactors 
              cross(pars->ppx[i], d[1]) *  //Parton level cross section
              ( 
              pars->ppdf->parton( -1*pdg,  x1, pars->Ecoll )/(x1*x1) * //PDFs with antiquark in first proton
              pars->ppdf->parton(    pdg,  x2, pars->Ecoll )/x2
              + 
              pars->ppdf->parton(    pdg,  x1, pars->Ecoll )/(x1*x1) * //Mirror: PDFs with antiquark in second proton   
              pars->ppdf->parton( -1*pdg,  x2, pars->Ecoll )/x2 
              ) *
              pars->psmear(d[1], d[0]);   //Smearing function
    }
    return sum;
  }
  
  
  //Cubature adaptor (Johnson)
  template<class PartialCrossX> 
  inline int integrand_cubature(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
  {
    double * x_temp= (double *) x;
//     std::cout<<"("<<x[0]<<", "<<x[1]<<", "<<x[2]<<")\n";   //#######################################################  DEBUG  ##################################
    size_t dim = ndim;
    struct parameter_set* pars = (struct parameter_set *)fdata;
//     std::cout<<pars->narr<<endl;   //#######################################################  DEBUG  ##################################
    fval[0]= integrand<PartialCrossX>(x_temp, dim, pars);    
//     std::cout<<fval[0]<<"\n"; //#######################################################  DEBUG  ##################################
    return 0;
  }
  
  
  
  
  //Cuba adaptor (Hahn)
  template<class PartialCrossX> 
  inline int integrand_cuba(const int* ndim, const double *x, const int* fdim, double *fval, void *fdata)
  {
    struct parameter_set_hahn* cpars = (struct parameter_set_hahn *)fdata;
    size_t dim = *ndim;
    double diff1 = cpars->high[0] - cpars->low[0];
    double diff2 = cpars->high[1] - cpars->low[1];
    double diff3 = cpars->high[2] - cpars->low[2];
    double point[3]={diff1*x[0] + cpars->low[0], diff2*x[1] + cpars->low[1], diff3*x[2] + cpars->low[2]};
    //Return integral
    fval[0] = integrand<PartialCrossX>( point, dim, cpars->pars ) * diff1 * diff2 * diff3 ;
//     std::cout<<fval[0]<<"\n"; //#######################################################  DEBUG  ##################################
    return 0;
  }
  
  
  
  
  
  template<class PartialCrossX> 
  double pheno::HadronXSec::detectedXsec(double el, double eh, double acc, int strat, double (* psmear)(double, double))
  {
    
    //Construct parameters
    PartonXSec* pp[]= {dxsec, uxsec, sxsec, cxsec, bxsec};                    
    struct parameter_set int_pars = {5, pp, pdf, psmear, Epp};  
    double epsabs(0.), epsrel(acc);
    const int dimint(3), dimres(1);

    //Integration variables
    double result, error, prev_res, diff;
    double xl[3] = {el, 1, 0};
    double xu[3] = {eh, Epp, 1};
    
    //CUBA parameters
    struct parameter_set_hahn int_pars_hahn = {&int_pars, xl, xu}; 
    int nregions, neval, fail;
    double integral[dimres], err[dimres], prob[dimres];
    
    if(strat==1)
    {
      // VEGAS Monte Carlo Integration
      //Define a gsl_monte_function to pass to the integrator
      gsl_monte_function F;
      F.f = &integrand<PartialCrossX>;
      F.dim = 3;
      F.params = &int_pars;      

      //Initialize random number gen for monte carlo
      gsl_rng_env_setup ();
      const gsl_rng_type *T;
      gsl_rng *r;
      T = gsl_rng_default;
      r = gsl_rng_alloc (T);
      
      //Integration
      gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(3);
      s->stage=0; // 0= new uniform grid
      s->mode=GSL_VEGAS_MODE_IMPORTANCE;  //Can pick between importance or stratified sampling
      gsl_monte_vegas_integrate (&F, xl, xu, dimint, 1000, r, s, &result, &error);
      s->stage=1; //1= keep grid from previous run
      diff=1.0; //set to 100% to start loop
      while(diff>epsrel)
      {
        prev_res=result;
        gsl_monte_vegas_integrate (&F, xl, xu, 3, calls, r, s, &result, &error);
        diff= fabs((prev_res-result)/result); //relative difference of last 2 iterations
        std::cout<<"Integral "<<result<<" Relative difference: "<<diff<<std::endl; //######################## DEBUG
      }
      //Free integration memory
      gsl_monte_vegas_free(s);
    }
    else if(strat==2)
    {
      hcubature(dimres, &integrand_cubature<PartialCrossX>, &int_pars, dimint, xl, xu, 0, epsabs, epsrel, ERROR_INDIVIDUAL, &result, &error);
      std::cout<<"Integral "<<result<<" Error: "<<error<<std::endl; //######################## DEBUG
    }
    else if(strat==3)
    {
//       std::cout<<"Before integration"<<std::endl; //######################## DEBUG
      Cuhre(dimint, dimres, &integrand_cuba<PartialCrossX>, &int_pars_hahn, 1, epsrel, epsabs, 0|4 , 0, 50000, 11, "", NULL, &nregions, &neval, &fail, integral, err, prob);
//       std::cout<<"After integration"<<std::endl; //######################## DEBUG
      result = integral[0];
      error = err[0];
      std::cout<<"Integral "<<result<<" Error: "<<error<<std::endl; //######################## DEBUG
    }
    else if(strat==4)
    {
      Suave(dimint, dimres, &integrand_cuba<PartialCrossX>, &int_pars_hahn, 1, epsrel, epsabs, 0 | 4, 0,   0, 50000, 1000, 2, 25, "", NULL,  &nregions, &neval, &fail, integral, err, prob);
      result = integral[0];
      error = err[0];
      std::cout<<"Integral "<<result<<" Error: "<<error<<std::endl; //######################## DEBUG
    }

    
    return result;
  }
  

  
}

#endif