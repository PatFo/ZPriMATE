#ifndef XSEC_H
#define XSEC_H



#include <config.h>
#include <stdexcept>
#include <vector>
#include <iostream> //#########################DEBUG
//PDF package
#include <mstwpdf.h>

//If GSL is present use QAG integration
#if HAVE_LIBGSLCBLAS
  //Numerical integration
  #include <gsl/gsl_integration.h>
#endif

//Cuba integration methods
#include <cuba.h>

#include "pheno.h"



#define EPSABS 0.0

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
    // No constructor/destructor needed since ZpModel exists only once -> c++11 unique_ptr?

    
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
    HadronXSec(pheno::HadronXSec &copy);
    HadronXSec & operator= (pheno::HadronXSec &assignment);
    ~HadronXSec();  
  };
  
  
  
  
  
  
  
  
  
  //IMPLEMENTATION OF MEMBER TEMPLATE  
  //---------------------------------------------------------------------------------------------------------
  
  
  

  //Parameter struct for differential cross section function
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
  
  
  

  //Cuba integrand
  inline int num_cuba(const int* ndim, const double *x, const int* fdim, double *fval, void *fdata)
  {
    //Get initialization parameters
    struct parameters * pars = (struct parameters *)fdata;
    
    //Integration boundaries
    double xl = pars->Ecm*pars->Ecm/(pars->Ecoll*pars->Ecoll);  //lower integration limit
    double xu = 1;  //Upper integration limit    
    
    //Define Bjorken x
    double diff1 = xu - xl;
    double x1 = diff1*x[0] + xl;
    double x2 = xl/x1;
    
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
    fval[0] = sum * diff1;
    
    return 0;
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

//If GSL is present use QAG method    
#if HAVE_LIBGSLCBLAS    

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

//If GSL not present do fake Cuba 2d integration
#else    
    
    //CUBA parameters
    int nregions, neval, fail;
    const int dimint(2), dimres(1);
    double integral[dimres], err[dimres], prob[dimres];
    Cuhre(dimint, dimres, &num_cuba, &local_pars, 1, 0.02, EPSABS, 0|4 , 0, 10000, 9, "", NULL, &nregions, &neval, &fail, integral, err, prob);
//     Vegas(dimint, dimint, &num_cuba, &local_pars, 1, 1e-1, EPSABS, 1|4, 0, 0, 5000, 1000, 500, 1000, 0, NULL, 0,  &neval, &fail, integral, err, prob);
//     Suave(dimint, dimres, &num_cuba, &local_pars, 1, accuracy_goal, EPSABS, 0 | 4, 0,   0, 50000, 1000, 2, 3, "", NULL,  &nregions, &neval, &fail, integral, err, prob);
//     std::cout<<"Evaluations: "<<neval<<std::endl<<"Exit status: "<<fail<<std::endl; //######################## DEBUG
    result = integral[0];
    
#endif    
    
    return result;
  }
  
  


  
  //Total binned cross section integration (theoretical) 
  //-------------------------------------------------------------------------------------------
  
  //Parameter struct for integrable function
  struct parameter_set{ int narr; PartonXSec** ppx; c_mstwpdf* ppdf; double (* psmear)(double, double); double Ecoll; double*  low; double*  high;};
  
  

  
  
  //Integrand for theoretical hadron cross section (without smearing)
  template<class PartialCrossX> 
  inline int integrand_theo(const int* ndim, const double *x, const int* fdim, double *fval, void *fdata)
  {
    //Get initialization parameters
    PartialCrossX cross;
    struct parameter_set* pars = (struct parameter_set*) fdata;
    
    double diff1 = pars->high[0] - pars->low[0];
    double diff2 = pars->high[1] - pars->low[1];
    double d[2]={diff1*x[0] + pars->low[0], diff2*x[1] + pars->low[1]};
    
    
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
    fval[0] *= diff1*diff2;
    return 0;
  }
  
  
  
  
  template<class PartialCrossX> 
  double pheno::HadronXSec::theoXsec(double el, double eh, double acc)
  {    
    //Construct parameters
    PartonXSec* pp[]= {dxsec, uxsec, sxsec, cxsec, bxsec};                    
    //Integration boundaries
    double xl[2] = {el, 0};
    double xu[2] = {eh, 1};    
    
    //CUBA parameters
    struct parameter_set int_pars = {5, pp, pdf, NULL, Epp, xl, xu};  
    int nregions, neval, fail;
    const int dimint(2), dimres(1);
    double integral[dimres], err[dimres], prob[dimres];
 
    Cuhre(dimint, dimres, &integrand_theo<PartialCrossX>, &int_pars, 1, acc, EPSABS, 0|4 , 0, 50000, 11, "", NULL, &nregions, &neval, &fail, integral, err, prob);
    std::cout<<"Integral "<<integral[0]<<" Error: "<<err[0]<<std::endl; //######################## DEBUG

    return integral[0];
  }
  
  
  
  
  
  //Total binned cross section integration (detectable/smeared) 
  //-------------------------------------------------------------------------------------------
  
  
  
  
  //Cuba adaptor (Hahn)
  template<class PartialCrossX> 
  inline int integrand_cuba(const int* ndim, const double *x, const int* fdim, double *fval, void *fdata)
  {
    PartialCrossX cross;
    struct parameter_set* pars = (struct parameter_set*)fdata;
    size_t dim = *ndim;
    double diff1 = pars->high[0] - pars->low[0];
    double diff2 = pars->high[1] - pars->low[1];
    double diff3 = pars->high[2] - pars->low[2];
    double point[3]={diff1*x[0] + pars->low[0], diff2*x[1] + pars->low[1], diff3*x[2] + pars->low[2]};

    
    //Define Bjorken x
    double scoll = pars->Ecoll * pars->Ecoll;
    double dx1dy = (scoll - point[1]*point[1])/scoll;
    double x1 = dx1dy * point[2] + point[1]*point[1]/scoll; // Variable trafo from x1 to y=a*x1 + b
    double x2 = point[1]*point[1]/(x1 * scoll);
    
    //Calculate sum of cross sections
    double sum = 0.;
    for(int i=0; i<pars->narr; ++i)
    {
      double pdg = (pars->ppx[i])->pdg_in();
      sum +=  dx1dy * 2 * point[1]/scoll *     //Prefactors 
              cross(pars->ppx[i], point[1]) *  //Parton level cross section
              ( 
              pars->ppdf->parton( -1*pdg,  x1, pars->Ecoll )/(x1*x1) * //PDFs with antiquark in first proton
              pars->ppdf->parton(    pdg,  x2, pars->Ecoll )/x2
              + 
              pars->ppdf->parton(    pdg,  x1, pars->Ecoll )/(x1*x1) * //Mirror: PDFs with antiquark in second proton   
              pars->ppdf->parton( -1*pdg,  x2, pars->Ecoll )/x2 
              ) *
              pars->psmear(point[1], point[0]);   //Smearing function
    }
    fval[0] = sum * diff1 * diff2 * diff3 ;

    return 0;
  }
  
  
  
  
  
  template<class PartialCrossX> 
  double pheno::HadronXSec::detectedXsec(double el, double eh, double acc, int strat, double (* psmear)(double, double))
  {
    //Integration boundaries
    double xl[3] = {el, 1, 0};
    double xu[3] = {eh, Epp, 1};
    
    //CUBA parameters
    PartonXSec* pp[]= {dxsec, uxsec, sxsec, cxsec, bxsec};                    
    struct parameter_set int_pars = {5, pp, pdf, psmear, Epp, xl, xu};  
    //Construct parameters
    int nregions, neval, fail;
    const int dimint(3), dimres(1);
    double integral[dimres], err[dimres], prob[dimres];
    

    if(strat==1)
    {
      Cuhre(dimint, dimres, &integrand_cuba<PartialCrossX>, &int_pars, 1, acc, EPSABS, 0|4 , 0, 50000, 11, "", NULL, &nregions, &neval, &fail, integral, err, prob);
      std::cout<<"Integral "<<integral[0]<<" Error: "<<err[0]<<std::endl; //######################## DEBUG
    }
    else if(strat==2)
    {
      Suave(dimint, dimres, &integrand_cuba<PartialCrossX>, &int_pars, 1, acc, EPSABS, 0 | 4, 0,   0, 50000, 1000, 2, 25, "", NULL,  &nregions, &neval, &fail, integral, err, prob);
      std::cout<<"Integral "<<integral[0]<<" Error: "<<err[0]<<std::endl; //######################## DEBUG
    }
    
    return integral[0];
  }
  

  
}

#endif
