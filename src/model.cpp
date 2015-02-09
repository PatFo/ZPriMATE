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
fermion::fermion(int fam, float t3, float m, float emc, float xlc, float xrc)
{
  family=fam;
  iso3=t3;
  
  mass=m;
  emcharge=emc;
  xlcharge=xlc;
  xrcharge=xrc;
}


//Read parameters of fermion object

int fermion::get_family()
{
  return family;
}

float fermion::get_iso3()
{
  return iso3;
}

float fermion::get_emcharge()
{
  return emcharge;
}

float fermion::get_xlcharge()
{
  return xlcharge;
}

float fermion::get_xrcharge()
{
  return xrcharge;
}

float fermion::get_mass()
{
  return mass;
}

//Set parameters of fermion object

void fermion::update_emcharge(float emc)
{
  emcharge=emc;
}

void fermion::update_xlcharge(float xlc)
{
  xlcharge=xlc;
}

void fermion::update_xrcharge(float xrc)
{
  xrcharge=xrc;
}

void fermion::change_mass(float new_mass)
{
  mass=new_mass;
}




//***************************************//
//            MODEL BASE CLASS           //
//***************************************//


//Constructor
bsm_parameters::bsm_parameters(float cpl, float mass, float mix)
{
  //Calculate SM parameters
  e=sqrt(4*M_PI*aew);
  g1=e/sqrt(1-sw2);
  g2=e/sqrt(sw2);
  gz=g2/sqrt(1-sw2);
  vev=2*mw*sqrt(sw2)/e;
  //Set BSM parameters
  gx=cpl;
  mzp=mass;
  mixing_angle=mix;
  //Calculate BSM parameters
  //fef=(MZp Cos[chi] Sqrt[4 MZp^2 - (gw/cw)^2 vev^2 - g1^2 vev^2 Tan[chi]^2])/(gx Sqrt[8 MZp^2 - 2 (gw/cw)^2 vev^2])
  fef= mzp* cos(mixing_angle)* sqrt( 4* pow(mzp,2) - pow(gz,2) * pow(vev,2) - pow(g1,2)*pow(vev,2)*pow(tan(mixing_angle),2) )/( gx * sqrt(8*pow(mzp,2) - 2*pow(gz,2)*pow(vev,2)) );
  //If mixing angle is 0 also xi should be 0
  if(mix==0)
  {
    xi=0;
  }else{
    //xi=ArcTan[ (2 g1 gw/cw Tan[chi] vev^2) / (8 gx^2 Sec[chi]^2 fef^2 + vev^2 (-(gw/cw)^2 +g1^2 Tan[chi]^2)) ]/2
    xi=atan( (2*g1*gz*tan(mixing_angle)*pow(vev,2))/(8*pow(gx*fef/cos(mixing_angle),2) + pow(vev,2)*( - pow(gz,2)+pow(g1*tan(mixing_angle),2) ) ) )/2;
  }
}

//Read parameters of model
//SM parameters

float bsm_parameters::e_()
{
  return e;
}

float bsm_parameters::g1_()
{
  return g1;
}

float bsm_parameters::g2_()
{
  return g2;
}

float bsm_parameters::gz_()
{
  return gz;
}

float bsm_parameters::vev_()
{
  return vev;
}

float bsm_parameters::aew_()
{
  return aew;
}

float bsm_parameters::as_()
{
  return as;
}

float bsm_parameters::mw_()
{
  return mw;
}

float bsm_parameters::mz_()
{
  return mz;
}

float bsm_parameters::wz_()
{
  return wz;
}

float bsm_parameters::sw2_()
{
  return sw2;
}



//BSM

float bsm_parameters::gx_()
{
  return gx;
}

float bsm_parameters::mzp_()
{
  return mzp;
}

float bsm_parameters::mixing_()
{
  return mixing_angle;
}

float bsm_parameters::fef_()
{
  return fef;
}

float bsm_parameters::xi_()
{
  return xi;
}



//Set parameters

void bsm_parameters::set_gx(float g)
{
  gx=g;
}

void bsm_parameters::set_mzp(float m)
{
  mzp=m;
}

void bsm_parameters::set_mixing(float mix)
{
  mixing_angle=mix;
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
  float czl  =  f.get_iso3() - (hypl/2 + f.get_iso3())*paras.sw2_();
  float czr  = -hypr/2 * paras.sw2_();
  float czpl = -paras.g1_()/2*hypl*tan(paras.mixing_()) + paras.gx_()*f.get_xlcharge()/cos(paras.mixing_()); //-g1/2*hypl*Tan[chi] + gx*QxqL/Cos[chi]
  float czpr = -paras.g1_()/2*hypr*tan(paras.mixing_()) + paras.gx_()*f.get_xrcharge()/cos(paras.mixing_()); 
  
  //Calculate full couplings after mixing
  q_gam = f.get_emcharge()*paras.e_();
  q_zl  = paras.gz_()*czl*cos(paras.xi_()) - czpl*sin(paras.xi_()); //gz*kZiL*Cos[xi] - kZPiL*Sin[xi]
  q_zr  = paras.gz_()*czr*cos(paras.xi_()) - czpr*sin(paras.xi_());
  q_zpl = paras.gz_()*czl*sin(paras.xi_()) + czpl*cos(paras.xi_()); //gz*kZiL*Sin[xi] + kZPiL*Cos[xi]
  q_zpr = paras.gz_()*czr*sin(paras.xi_()) + czpr*cos(paras.xi_());
}




//*************************************************************//
//          FERMION CLASS EXTENDED BY VECTOR COUPLING          //
//*************************************************************//


//Constructor of extended fermion class: initialize base class
fermionExt::fermionExt(bool massive, int fam, float t3, float m, float emc, float xlc, float xrc): fermion(fam, t3, m, emc, xlc, xrc), pvecc() 
{
  if(!massive)
  {
    this->change_mass(0);
  }

}




//Destructor: make sure to free vcoef pointer if it has been assigned a value
fermionExt::~fermionExt()
{
  if(pvecc)
  {
    std::cout<<"Deleting pointer of type vcoeff\n";
    delete pvecc;
  }
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



