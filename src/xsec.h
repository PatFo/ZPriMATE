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
//       template<class PartialCrossX> double binnedXsec(double el, double eh, double accuracy);
      template<class PartialCrossX> double binnedXsec(double el, double eh, double accuracy);
      //Full experimentally detectable cross section
      template<class PartialCrossX> double detectedXsec(double el, double eh, double acc, int strategy, double (* psmear)(double, double));
    public:
      void set_accuracy(double accuracy);
      void set_monte_calls(size_t int_calls);
      //Hadronic differential cross sections dSig/dm
      double dsigSM    (double Ecm, unsigned int int_strategy=1);
      double dsigInt   (double Ecm, unsigned int int_strategy=1);
      double dsigSignal(double Ecm, unsigned int int_strategy=1);
      double dsigTotal (double Ecm, unsigned int int_strategy=1);
      //Member function that fills EMPTY vector with sigSM, sigInt, sigSignal, sigTotal 
      //ALWAYS use this function if more than one of these cross sections is needed at a time
      void crossSections (double Ecm, std::vector<double> * results, unsigned int int_strategy=1);
      //Full pure hadronic cross sections
      double totXsec(double el, double eh, double accuracy, double (* psmear)(double,double)=NULL, int strategy=1);
      double zpXsec(double el, double eh, double accuracy, double (* psmear)(double,double)=NULL, int strategy=1);
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
    //If no smearing function was defined just integrate the pure hadron diff cross section
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
    size_t bisection_lim = 1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (bisection_lim);
    gsl_integration_qags (&F, el, eh, 0, accuracy, bisection_lim, w, &result, &error); 
    //Free memory of integration workspace
    gsl_integration_workspace_free(w);
    
    return result;
  }

  
  
  
  //Multi dim Monte Carlo implementation 
  //-------------------------------------------------------------------------------------------
  
  //Parameter struct for integrable function
  struct parameter_set{ int narr; PartonXSec** ppx; c_mstwpdf* ppdf; double (* psmear)(double, double); double Ecoll; };
  
  struct cuba_par{ double** interval; parameter_set pars; };
  
  
  
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
  inline int integrand_cuba(const int *ndim, const double x[], const int *ncomp, double integral[], void *userdata)
  {
    double * x_temp= (double *) x;
    size_t dim = (size_t) ndim;
    struct cuba_par* cpars = (struct cuba_par *) userdata;
    double** coeffs = cpars->interval;
    double jacobi =1; //Jacobi determinant of trafo
    //Linear trafo of coordinates
    for(unsigned int i=0; i<3; ++i)
    {
      double dydx = (coeffs[1][i]-coeffs[0][i]);
      jacobi *= dydx;
      x_temp[i] = dydx*x_temp[i] + coeffs[0][i]; 
    }
//     std::cout<<"("<<x_temp[0]<<", "<<x_temp[1]<<", "<<x_temp[2]<<")\n";   //#######################################################  DEBUG  ##################################
//     std::cout<<pars->narr<<endl;   //#######################################################  DEBUG  ##################################
    integral[0]= integrand<PartialCrossX>(x_temp, dim, &(cpars->pars)) * jacobi;
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
    const int dimint(3), dimres(1), mineval(0), maxeval(100000);

    //Integration variables
    double result, error, prev_res, diff;
    double xl[3] = {el, 1, 0};
    double xu[3] = {eh, Epp, 1};
    
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
  //     s->alpha=1.6;
      gsl_monte_vegas_integrate (&F, xl, xu, dimint, 1000, r, s, &result, &error);
  //     double chisqdiff =fabs (gsl_monte_vegas_chisq(s) - 1.0);
  //     std::cout<<"Integral "<<result<<" Reduced chisquare: "<<chisqdiff<<std::endl; //######################## DEBUG
      s->stage=1; //1= keep grid from previous run
      diff=1.0; //set to 100% to start loop
      while(diff>epsrel)
      {
        prev_res=result;
        gsl_monte_vegas_integrate (&F, xl, xu, 3, calls, r, s, &result, &error);
  //       chisqdiff= fabs (gsl_monte_vegas_chisq(s) - 1.0); //Check whether chisq/dof is consistent with 1 
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
    else if(strat==3)  //CUHRE integration
    {
      const int key(11);
      int nregions, neval, fail;
      cubareal integral[dimres], aerr[dimres], prob[dimres];
      //Constructing parameter space
      double* interval[2] ={xl, xu};
      struct cuba_par cpars = {interval, int_pars};
      Cuhre(dimint, dimres, integrand_cuba<PartialCrossX>, &cpars, 1, epsrel , epsabs, 0, mineval, maxeval, key, NULL, NULL, &nregions, &neval, &fail, integral, aerr, prob );
      result= integral[0];
      std::cout<<"Integral "<<result<<std::endl; //######################## DEBUG
    }
    else if(strat==4)  //SUAVE integration
    {
      int nregions, neval, fail;
      cubareal integral[dimres], aerr[dimres], prob[dimres];
      //Constructing parameter space
      double* interval[2] ={xl, xu};
      struct cuba_par cpars = {interval, int_pars};
      Suave(dimint, dimres, integrand_cuba<PartialCrossX>, &cpars, 1, epsrel , epsabs, 0, 0, mineval, maxeval, 10, 5, 0.4,  NULL, NULL, &nregions, &neval, &fail, integral, aerr, prob);
      result= integral[0];
      std::cout<<"Integral "<<result<<std::endl; //######################## DEBUG
    }
    
    return result;
  }
  

  
}

#endif