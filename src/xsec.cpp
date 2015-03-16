#define _USE_MATH_DEFINES

#include "xsec.h"
#include <complex>
#include <cstdio>
//PDF package
#include <mstwpdf.h>




//Function that calculates the kinematical prefactor for the cross section coming from the Breit-Wigner propagators
double kinematics(double Ecm, double m1, double m2, double w1, double w2, bool sum)
{
  std::complex<double> kin1(Ecm*Ecm-m1*m1,0);
  std::complex<double> kin2(Ecm*Ecm-m2*m2,0);
  std::complex<double> bw1(0,w1*m1);
  std::complex<double> bw2(0,w2*m2);
  std::complex<double> one(1,0);
  
  //Calculate kinematical factor
  std::complex<double> kinematics = one/( (kin1-bw1)*(kin2+bw2) );
  if(sum)
  {
    kinematics += one/( (kin1+bw1)*(kin2-bw2) );
  }
//   std::printf("Kinematics: (%g, %g)\n", std::real(kinematics), std::imag(kinematics));
  return std::real(kinematics);
}





//Function that sets the coupling structure
double couplings(double cinL1, double cinL2, double cinR1, double cinR2, double coutL1, double coutL2, double coutR1, double coutR2)
{
  return (cinL1*cinL2+cinR1*cinR2)*(coutL1*coutL2+coutR1*coutR2);
}







//*******************************************************//
//            Class for partonic cross sections          //
//*******************************************************//



//Constructor: Initialize the coupling-prefactors
pheno::PartonXSec::PartonXSec(fundamental::fermionExt* f_in, fundamental::fermionExt* f_out, pheno::ZpModel* p_model): _model(p_model)
{
  //Store for every partial cross section the prefactor consisting of couplings, color factors and phase space constants
  numGam   = f_in->Nc()*f_out->Nc()/(48*M_PI) * couplings(f_in->vecc().q_gam, 
                                                          f_in->vecc().q_gam, 
                                                          f_in->vecc().q_gam, 
                                                          f_in->vecc().q_gam, 
                                                          f_out->vecc().q_gam,
                                                          f_out->vecc().q_gam,
                                                          f_out->vecc().q_gam,
                                                          f_out->vecc().q_gam);
  
  numZ     = f_in->Nc()*f_out->Nc()/(48*M_PI) * couplings(f_in->vecc().q_zl, 
                                                          f_in->vecc().q_zl, 
                                                          f_in->vecc().q_zr, 
                                                          f_in->vecc().q_zr, 
                                                          f_out->vecc().q_zl,
                                                          f_out->vecc().q_zl,
                                                          f_out->vecc().q_zr,
                                                          f_out->vecc().q_zr);
  
  numZp    = f_in->Nc()*f_out->Nc()/(48*M_PI) * couplings(f_in->vecc().q_zpl,
                                                          f_in->vecc().q_zpl,
                                                          f_in->vecc().q_zpr,
                                                          f_in->vecc().q_zpr,
                                                          f_out->vecc().q_zpl,
                                                          f_out->vecc().q_zpl,
                                                          f_out->vecc().q_zpr,
                                                          f_out->vecc().q_zpr);
  
  numGamZ  = f_in->Nc()*f_out->Nc()/(48*M_PI) * couplings(f_in->vecc().q_gam, 
                                                          f_in->vecc().q_zl, 
                                                          f_in->vecc().q_gam, 
                                                          f_in->vecc().q_zr, 
                                                          f_out->vecc().q_gam,
                                                          f_out->vecc().q_zl,
                                                          f_out->vecc().q_gam,
                                                          f_out->vecc().q_zr);
  
  numGamZp = f_in->Nc()*f_out->Nc()/(48*M_PI) * couplings(f_in->vecc().q_gam, 
                                                          f_in->vecc().q_zpl, 
                                                          f_in->vecc().q_gam, 
                                                          f_in->vecc().q_zpr, 
                                                          f_out->vecc().q_gam,
                                                          f_out->vecc().q_zpl,
                                                          f_out->vecc().q_gam,
                                                          f_out->vecc().q_zpr);
  
  numZZp   = f_in->Nc()*f_out->Nc()/(48*M_PI) * couplings(f_in->vecc().q_zl, 
                                                          f_in->vecc().q_zpl, 
                                                          f_in->vecc().q_zr, 
                                                          f_in->vecc().q_zpr, 
                                                          f_out->vecc().q_zl,
                                                          f_out->vecc().q_zpl,
                                                          f_out->vecc().q_zr,
                                                          f_out->vecc().q_zpr);
  pdgin = int(f_in->get_pdg());
}






//Elementary Cross Section Pieces
//Calculate pure photonic cross section
double pheno::PartonXSec::sigGam(double Ecm)
{
  return GeV2fb*numGam*kinematics(Ecm, 0, 0, 0, 0, false)*Ecm*Ecm;
}


//Calculate pure Z cross section
double pheno::PartonXSec::sigZ(double Ecm)
{
  return GeV2fb*numZ*kinematics(Ecm, _model->mz_(), _model->mz_(), _model->wz_(), _model->wz_(), false)*Ecm*Ecm;
}


//Calculate pure Zp cross section
double pheno::PartonXSec::sigZp(double Ecm)
{
  return GeV2fb*numZp*kinematics(Ecm, _model->mzp_(), _model->mzp_(), _model->wzp_(), _model->wzp_(), false)*Ecm*Ecm;
}


//Calculate interference of Photon and Z 
double pheno::PartonXSec::sigGamZ(double Ecm)
{
  return GeV2fb*numGamZ*kinematics(Ecm, 0, _model->mz_(), 0, _model->wz_(), true)*Ecm*Ecm;
}


//Calculate interference of Photon and Z 
double pheno::PartonXSec::sigGamZp(double Ecm)
{
  return GeV2fb*numGamZp*kinematics(Ecm, 0, _model->mzp_(), 0, _model->wzp_(), true)*Ecm*Ecm;
}


//Calculate interference of Photon and Z 
double pheno::PartonXSec::sigZZp(double Ecm)
{
  return GeV2fb*numZZp*kinematics(Ecm, _model->mz_(), _model->mzp_(), _model->wz_(), _model->wzp_(), true)*Ecm*Ecm;
}




//Composite Cross Sections
//Calculate the cross section in the SM
double pheno::PartonXSec::sigSM(double Ecm)
{
  return sigGam(Ecm)+sigGamZ(Ecm)+sigZ(Ecm);
}


//Calculate the pure Zp induced interference
double pheno::PartonXSec::sigInt(double Ecm)
{
  return sigGamZp(Ecm)+sigZZp(Ecm);
}


//Calculate the total cross section in the Zp model
double pheno::PartonXSec::sigTot(double Ecm)
{
  return sigGam(Ecm)+sigGamZ(Ecm)+sigZ(Ecm)+sigGamZp(Ecm)+sigZZp(Ecm)+sigZp(Ecm);
}


//Fast method for multiple cross sections
void pheno::PartonXSec::crossSections(double Ecm, std::vector< double >* results, unsigned int int_strategy)
{
  double xsm = sigSM(Ecm);
  double xint = sigInt(Ecm);
  double xzp = sigZp(Ecm);
  
  results->push_back(xsm+xint+xzp);
  results->push_back(xsm+xint);
  results->push_back(xsm+xzp);
  results->push_back(xsm);
}




//Return pdg code of in-particle
int pheno::PartonXSec::pdg_in()
{
  return pdgin;
}










//*******************************************************//
//            Class for hadronic cross sections          //
//*******************************************************//




//Functors for partial cross section calculataions
//------------------------------------------------------------
//SM cross section
struct SigSM{
  double operator() (pheno::PartonXSec* pxsec, double Ecm)
  {
    return pxsec->sigSM(Ecm);
  }
};


//Interference cross section
struct SigInt{
  double operator() (pheno::PartonXSec* pxsec, double Ecm)
  {
    return pxsec->sigInt(Ecm);
  }
};


//Zp cross section
struct SigZp{
  double operator() (pheno::PartonXSec* pxsec, double Ecm)
  {
    return pxsec->sigZp(Ecm);
  }
};


//Total cross section
struct SigTot{
  double operator() (pheno::PartonXSec* pxsec, double Ecm)
  {
    return pxsec->sigSM(Ecm) + pxsec->sigInt(Ecm) + pxsec->sigZp(Ecm);
  }
};



//CLASS IMPLEMENTATION
//-------------------------------------------------------------
//Constructor
pheno::HadronXSec::HadronXSec(fundamental::fermionExt* f_out, pheno::ZpModel* p_model, char* pdf_grid_file, double Ecoll)
{
  accuracy_goal = 1e-2;  //Default numerica integ accuracy
  calls = 100; //Default value for calls per monte carlo integration point
  Epp = Ecoll;
  //Allocate parpxsec->sigSM(Ecm)tonic cross sections
  dxsec = new pheno::PartonXSec(&p_model->d, f_out, p_model);
  uxsec = new pheno::PartonXSec(&p_model->u, f_out, p_model);
  sxsec = new pheno::PartonXSec(&p_model->s, f_out, p_model);
  cxsec = new pheno::PartonXSec(&p_model->c, f_out, p_model);
  bxsec = new pheno::PartonXSec(&p_model->b, f_out, p_model);
  //Allocate pdf object
  pdf = new c_mstwpdf(pdf_grid_file);
}


pheno::HadronXSec::~HadronXSec()
{
  //Free all pointers
  delete dxsec;
  delete uxsec;
  delete sxsec;
  delete cxsec;
  delete bxsec;
  delete pdf;
}



//Set the number of calls of the cross section function per integration point
void pheno::HadronXSec::set_monte_calls(size_t int_calls)
{
  calls=int_calls;
}


//Set relative accuracy for numerical integration
void pheno::HadronXSec::set_accuracy(double accuracy)
{
  accuracy_goal=accuracy;
}




double pheno::HadronXSec::sigSM(double Ecm, unsigned int int_strategy)
{
  return pdfconvoluted<SigSM>(Ecm, int_strategy);
}



double pheno::HadronXSec::sigInt(double Ecm, unsigned int int_strategy)
{
  return sigSM(Ecm, int_strategy) + pdfconvoluted<SigInt>(Ecm, int_strategy);
}



double pheno::HadronXSec::sigSignal(double Ecm, unsigned int int_strategy)
{
  return sigSM(Ecm, int_strategy) + pdfconvoluted<SigZp>(Ecm, int_strategy);
}



double pheno::HadronXSec::sigTotal(double Ecm, unsigned int int_strategy)
{
  return sigInt(Ecm, int_strategy) +  pdfconvoluted<SigZp>(Ecm, int_strategy);
}


//Fast method for multiple cross sections
void pheno::HadronXSec::crossSections(double Ecm, std::vector< double >* results, unsigned int int_strategy)
{
  double xsm = pdfconvoluted<SigSM>(Ecm, int_strategy);
  double xint = pdfconvoluted<SigInt>(Ecm, int_strategy);
  double xzp = pdfconvoluted<SigZp>(Ecm, int_strategy);
  
  results->push_back(xsm+xint+xzp);
  results->push_back(xsm+xint);
  results->push_back(xsm+xzp);
  results->push_back(xsm);
}


//Calculate the total hadronic cross section in the bounds [el, eh]
double pheno::HadronXSec::totXsec(double el, double eh, double accuracy)
{
  return binnedXsec<SigTot>(el, eh, accuracy);
}


//Calculate the zp hadronic cross section in the bounds [el, eh]
double pheno::HadronXSec::zpXsec(double el, double eh, double accuracy)
{
  return binnedXsec<SigZp>(el, eh, accuracy);
}
