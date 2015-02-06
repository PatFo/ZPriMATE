#include "pheno.h"

using namespace pheno;


//Implemenation of quark constructors
up::up(float xlc, float xrc): fermion(1, 1, 2.3e-3, 2./3, xlc, xrc) { }

down::down(float xlc, float xrc): fermion(1, -1, 4.8e-3, -1./3, xlc, xrc) { }

charm::charm(float xlc, float xrc): fermion(2, 1, 1.275, 2./3, xlc, xrc) { }

strange::strange(float xlc, float xrc): fermion(2, -1, 9.5e-2, -1./3, xlc, xrc) { }

top::top(float xlc, float xrc): fermion(3, 1, 173.5, 2./3, xlc, xrc) { }

bottom::bottom(float xlc, float xrc): fermion(3, -1, 4.18, -1./3, xlc, xrc) { }

