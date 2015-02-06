#include "pheno.h"

using namespace pheno;




//************************************//
//            QUARK CLASSES           //
//************************************//

///Implemenation of quark constructors

up::up(float xlc, float xrc): fermion(1, 1, 2.3e-3, 2./3, xlc, xrc) { }

down::down(float xlc, float xrc): fermion(1, -1, 4.8e-3, -1./3, xlc, xrc) { }

charm::charm(float xlc, float xrc): fermion(2, 1, 1.275, 2./3, xlc, xrc) { }

strange::strange(float xlc, float xrc): fermion(2, -1, 9.5e-2, -1./3, xlc, xrc) { }

top::top(float xlc, float xrc): fermion(3, 1, 173.5, 2./3, xlc, xrc) { }

bottom::bottom(float xlc, float xrc): fermion(3, -1, 4.18, -1./3, xlc, xrc) { }



//************************************//
//            LEPTON CLASSES          //
//************************************//


///Implemenation of lepton constructors

electron::electron(bool massless, float xlc, float xrc): fermion(1, -1, 5.11e-4, -1, xlc, xrc)
{
  if(massless)
  {
    this->change_mass(0);
  }
}

muon::muon(bool massless, float xlc, float xrc):fermion(2, -1, 0.1057, -1, xlc, xrc)
{
  if(massless)
  {
    this->change_mass(0);
  }
}

tauon::tauon(bool massless, float xlc, float xrc):fermion(3, -1, 1.777, -1, xlc, xrc)
{
  if(massless)
  {
    this->change_mass(0);
  }
}

nu_el::nu_el(float xlc, float xrc): fermion(1, 1, 0, 0, xlc, xrc) { }

nu_mu::nu_mu(float xlc, float xrc): fermion(2, 1, 0, 0, xlc, xrc) { }

nu_tau::nu_tau(float xlc, float xrc): fermion(3, 1, 0, 0, xlc, xrc) { }

