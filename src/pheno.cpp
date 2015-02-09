#include "pheno.h"
#include <iostream>

using namespace pheno;




//************************************//
//            QUARK CLASSES           //
//************************************//

///Implemenation of quark constructors

down::down(float xlc, float xrc, bool massive): fermionExt(massive, 1, -1./2, 4.8e-3, -1./3, xlc, xrc) { }

up::up(float xlc, float xrc, bool massive): fermionExt(massive, 1, 1./2, 2.3e-3, 2./3, xlc, xrc) { }

strange::strange(float xlc, float xrc, bool massive): fermionExt(massive, 2, -1./2, 9.5e-2, -1./3, xlc, xrc) { }

charm::charm(float xlc, float xrc, bool massive): fermionExt(massive, 2, 1./2, 1.275, 2./3, xlc, xrc) { }

bottom::bottom(float xlc, float xrc, bool massive): fermionExt(massive, 3, -1./2, 4.18, -1./3, xlc, xrc) { }

top::top(float xlc, float xrc, bool massive): fermionExt(massive, 3, 1./2, 173.5, 2./3, xlc, xrc) { }



//************************************//
//            LEPTON CLASSES          //
//************************************//


///Implemenation of lepton constructors

electron::electron(float xlc, float xrc, bool massive): fermionExt(massive, 1, -1./2, 5.11e-4, -1, xlc, xrc) { }

muon::muon(float xlc, float xrc, bool massive): fermionExt(massive, 2, -1./2, 0.1057, -1, xlc, xrc) { }

tauon::tauon(float xlc, float xrc, bool massive): fermionExt(massive, 3, -1./2, 1.777, -1, xlc, xrc) { }

nu_el::nu_el(float xlc, float xrc): fermionExt(false, 1, 1./2, 0, 0, xlc, xrc) { }

nu_mu::nu_mu(float xlc, float xrc): fermionExt(false, 2, 1./2, 0, 0, xlc, xrc) { }

nu_tau::nu_tau(float xlc, float xrc): fermionExt(false, 3, 1./2, 0, 0, xlc, xrc) { }




//************************************//
//            ZP-MODEL CLASS          //
//************************************//

zpmodel::zpmodel(const char* configfile): bsm_parameters(0.1, 1500, 0)
{

    conf_reader reader(configfile);
    dict init = reader.get_config();
    
    
  
}

