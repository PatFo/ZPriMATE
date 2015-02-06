#include "model.h"


using namespace fundamental;




//**************************************//
//         FERMION BASE CLASS           //
//**************************************//


//Implementation of constructor of fermion class
fermion::fermion(int fam, int t3, float m, float emc, float xlc, float xrc)
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

int fermion::get_iso3()
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
  gx=cpl;
  mzp=mass;
  mixing_angle=mix;
}

//Read parameters of model
//SM parameters

float bsm_parameters::get_aew()
{
  return aew;
}

float bsm_parameters::get_as()
{
  return as;
}

float bsm_parameters::get_mz()
{
  return mz;
}

float bsm_parameters::get_wz()
{
  return wz;
}

float bsm_parameters::get_sw2()
{
  return sw2;
}

//BSM

float bsm_parameters::get_gx()
{
  return gx;
}

float bsm_parameters::get_mzp()
{
  return mzp;
}

float bsm_parameters::get_mixing()
{
  return mixing_angle;
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







