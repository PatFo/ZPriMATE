#define _USE_MATH_DEFINES

#include <stdexcept>
#include <iostream> //Only for debugging
#include <cmath>
#include "model.h"


using namespace fundamental;




//**************************************//
//         FERMION BASE CLASS           //
//**************************************//


//Implementation of constructor of fermion class
fermion::fermion(int n_pdg, double t3, double m, double emc, double xlc, double xrc)
{
  pdg=n_pdg;
  iso3=t3;
  
  mass=m;
  emcharge=emc;
  xlcharge=xlc;
  xrcharge=xrc;
}


//Read parameters of fermion object

int fermion::get_pdg()
{
  return pdg;
}

double fermion::get_iso3()
{
  return iso3;
}

double fermion::get_emcharge()
{
  return emcharge;
}

double fermion::get_xlcharge()
{
  return xlcharge;
}

double fermion::get_xrcharge()
{
  return xrcharge;
}

double fermion::get_mass()
{
  return mass;
}

//Set parameters of fermion object

void fermion::update_emcharge(double emc)
{
  emcharge=emc;
}

void fermion::update_xlcharge(double xlc)
{
  xlcharge=xlc;
}

void fermion::update_xrcharge(double xrc)
{
  xrcharge=xrc;
}

void fermion::change_mass(double new_mass)
{
  mass=new_mass;
}




//***************************************//
//            MODEL BASE CLASS           //
//***************************************//


void bsm_parameters::calc_fef()
{
  fef= mzp* cos(mixing_angle)/gx* sqrt( (pow(gz, 2.0)+pow((g1*tan(mixing_angle)), 2.0)-pow((2*mzp/vev), 2.0))/(2*pow(gz, 2.0)-8*pow((mzp/vev), 2)) );
}

void bsm_parameters::calc_xi()
{
  xi=atan( (2*g1*gz*tan(mixing_angle))/(8*pow(gx*fef/(vev*cos(mixing_angle)), 2) - pow(gz,2) + pow(g1*tan(mixing_angle),2)  ) )/2;
}

void bsm_parameters::update()
{
  calc_fef();
  calc_xi();
}


//Constructor
bsm_parameters::bsm_parameters(double cpl, double mass, double mix)
{
  //Calculate SM parameters
  e=0.313451;//sqrt(4*M_PI*aew);
  g1=0.358072;//e/sqrt(1-sw2);
  g2=0.648397;//e/sqrt(sw2);
  gz=g2/sqrt(1-sw2);
  vev=2*mw*sqrt(sw2)/e;
  //Set BSM parameters
  gx=cpl;
  mzp=mass;
  mixing_angle=mix;
  //Calculate BSM parameters
  //Value         -> (MZp Cos[cchi] Sqrt[4 MZp^2 - (gw/cw)^2 vev^2 - g1^2 vev^2 Tan[cchi]^2])/(gx Sqrt[8 MZp^2 - 2 (gw/cw)^2 vev^2]),


  update();
}

//Read parameters of model
//SM parameters

double bsm_parameters::e_()
{
  return e;
}

double bsm_parameters::g1_()
{
  return g1;
}

double bsm_parameters::g2_()
{
  return g2;
}

double bsm_parameters::gz_()
{
  return gz;
}

double bsm_parameters::vev_()
{
  return vev;
}

double bsm_parameters::aew_()
{
  return aew;
}

double bsm_parameters::as_()
{
  return as;
}

double bsm_parameters::mw_()
{
  return mw;
}

double bsm_parameters::mz_()
{
  return mz;
}

double bsm_parameters::wz_()
{
  return wz;
}

double bsm_parameters::sw2_()
{
  return sw2;
}



//BSM

double bsm_parameters::gx_()
{
  return gx;
}

double bsm_parameters::mzp_()
{
  return mzp;
}

double bsm_parameters::mixing_()
{
  return mixing_angle;
}

double bsm_parameters::fef_()
{
  return fef;
}

double bsm_parameters::xi_()
{
  return xi;
}



//Set parameters

void bsm_parameters::set_gx(double g)
{
  gx=g;
  update();
}

void bsm_parameters::set_mzp(double m)
{
  mzp=m;
  update();
}

void bsm_parameters::set_mixing(double mix)
{
  mixing_angle=mix;
  update();
}



//*****************************************************//
//          VECTOR COUPLING COEFFICIENT CLASS          //
//*****************************************************//


vcoeff::vcoeff(fermion f, bsm_parameters paras)
{
  //Calculate hypercharge of particle
  hypl= 2*(f.get_emcharge()-f.get_iso3());
  hypr= 2*f.get_emcharge();
  
  
  //Calculate couplings before mixing
  double czl  =  f.get_iso3() - f.get_emcharge()*paras.sw2_();
  double czr  = -f.get_emcharge() * paras.sw2_();
  double czpl = -paras.g1_()/2*hypl*tan(paras.mixing_()) + paras.gx_()*f.get_xlcharge()/cos(paras.mixing_()); //-g1/2*hypl*Tan[chi] + gx*QxqL/Cos[chi]
  double czpr = -paras.g1_()/2*hypr*tan(paras.mixing_()) + paras.gx_()*f.get_xrcharge()/cos(paras.mixing_()); 
  
  //Calculate full couplings after mixing
  q_gam = f.get_emcharge()*paras.e_();
  q_zl  = paras.gz_()*czl*cos(paras.xi_()) - czpl*sin(paras.xi_()); //gz*kZiL*Cos[xi] - kZPiL*Sin[xi]
  q_zr  = paras.gz_()*czr*cos(paras.xi_()) - czpr*sin(paras.xi_());
  q_zpl = paras.gz_()*czl*sin(paras.xi_()) + czpl*cos(paras.xi_()); //gz*kZiL*Sin[xi] + kZPiL*Cos[xi]
  q_zpr = paras.gz_()*czr*sin(paras.xi_()) + czpr*cos(paras.xi_());
}


double vcoeff::get_hypl()
{
  return hypl;
}


double vcoeff::get_hypr()
{
  return hypr;
}



//*************************************************************//
//          FERMION CLASS EXTENDED BY VECTOR COUPLING          //
//*************************************************************//


//Constructor of extended fermion class: initialize base class
fermionExt::fermionExt(bool massive, int n_pdg, double t3, double m, double emc, double xlc, double xrc, int n): fermion(n_pdg, t3, m, emc, xlc, xrc), pvecc() 
{
  nc=n;
  is_massive = massive;
}




//Destructor: make sure to free vcoef pointer if it has been assigned a value
fermionExt::~fermionExt()
{
  if(pvecc)
  {
//     std::cout<<"Deleting pointer of type vcoeff\n";
    delete pvecc;
  }
}



int fermionExt::Nc()
{
  return nc;
}


double fermionExt::m()
{
  if(is_massive)return get_mass();
  else return 0;
}


//Handling of vector couplings

void fermionExt::set_vecc(vcoeff* ptr)
{
  pvecc=ptr;
}



vcoeff fermionExt::vecc()
{
  if(pvecc)
  {
    return *pvecc;
  }else{
    throw std::runtime_error("ERROR: Trying to access uninitialized instance of type vcoeff in fermion class!");
  }  
}



