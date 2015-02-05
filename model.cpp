#include "model.h"


using namespace fundamental;


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
//-----------------------------------------------
int fermion::get_family()
{
  return family;
}
//-----------------------------------------------
int fermion::get_iso3()
{
  return iso3;
}
//-----------------------------------------------
float fermion::get_emcharge()
{
  return emcharge;
}
//-----------------------------------------------
float fermion::get_xlcharge()
{
  return xlcharge;
}
//-----------------------------------------------
float fermion::get_xrcharge()
{
  return xrcharge;
}
//-----------------------------------------------
float fermion::get_mass()
{
  return mass;
}

//Set parameters of fermion object
//-----------------------------------------------
void fermion::update_emcharge(float emc)
{
  emcharge=emc;
}
//-----------------------------------------------
void fermion::update_xlcharge(float xlc)
{
  xlcharge=xlc;
}
//-----------------------------------------------
void fermion::update_xrcharge(float xrc)
{
  xrcharge=xrc;
}
//-----------------------------------------------
void fermion::change_mass(float new_mass)
{
  mass=new_mass;
}


