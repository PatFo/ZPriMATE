#define _USE_MATH_DEFINES

#include "pheno.h"
#include <cstdio>
#include <cmath>
#include <iostream>
using namespace pheno;




//************************************//
//            QUARK CLASSES           //
//************************************//

///Implemenation of quark constructors

down::down(double xlc, double xrc, bool massive): fermionExt(massive, 1, -1./2, 4.8e-3, -1./3, xlc, xrc, 3) { }

up::up(double xlc, double xrc, bool massive): fermionExt(massive, 2, 1./2, 2.3e-3, 2./3, xlc, xrc, 3) { }

strange::strange(double xlc, double xrc, bool massive): fermionExt(massive, 3, -1./2, 9.5e-2, -1./3, xlc, xrc, 3) { }

charm::charm(double xlc, double xrc, bool massive): fermionExt(massive, 4, 1./2, 1.275, 2./3, xlc, xrc, 3) { }

bottom::bottom(double xlc, double xrc, bool massive): fermionExt(massive, 5, -1./2, 4.18, -1./3, xlc, xrc, 3) { }

top::top(double xlc, double xrc, bool massive): fermionExt(massive, 6, 1./2, 173.5, 2./3, xlc, xrc, 3) { }  //PDG=173.5



//************************************//
//            LEPTON CLASSES          //
//************************************//


///Implemenation of lepton constructors

electron::electron(double xlc, double xrc, bool massive): fermionExt(massive, 11, -1./2, 5.11e-4, -1, xlc, xrc, 1) { }

muon::muon(double xlc, double xrc, bool massive): fermionExt(massive, 13, -1./2, 0.1057, -1, xlc, xrc, 1) { }

tauon::tauon(double xlc, double xrc, bool massive): fermionExt(massive, 15, -1./2, 1.777, -1, xlc, xrc, 1) { }

nu_el::nu_el(double xlc, double xrc): fermionExt(false, 12, 1./2, 0, 0, xlc, xrc, 1) { }

nu_mu::nu_mu(double xlc, double xrc): fermionExt(false, 14, 1./2, 0, 0, xlc, xrc, 1) { }

nu_tau::nu_tau(double xlc, double xrc): fermionExt(false, 16, 1./2, 0, 0, xlc, xrc, 1) { }




//************************************//
//            ZP-MODEL CLASS          //
//************************************//



//Initializes the internal fermion list
void ZpModel::setup_flst()
{
  //Set up list of pointers to fermions for map iteration
  flst["up"]=&u; flst["charm"]=&c; flst["top"]=&t;
  flst["down"]=&d; flst["strange"]=&s; flst["bottom"]=&b;
  
  flst["electron"]=&el; flst["muon"]=&mu; flst["tauon"]=&tau;  
  flst["nu_el"]=&ne; flst["nu_mu"]=&nm; flst["nu_tau"]=&nt;
}


// Helper function for identification of fermion type to assign proper matrices
// Returns pair indicating <type, familyIndex>
std::tuple<std::string,unsigned int,int> ZpModel::fermionType(std::string name) {
  std::string type;
  unsigned int familyInd;
  int pdgID;
  
  if(name=="up" || name == "charm" || name == "top") {
    type = "UP";
  }
  else if (name=="down" || name == "strange" || name == "bottom") {
    type = "DOWN";
  }
  else if (name == "electron" || name == "muon" || name == "tauon") {
    type = "LEP";
  }
  else if (name == "nu_el" || name == "nu_mu" || name == "nu_tau") {
    type = "NEU";
  } else {
    throw std::runtime_error("ERROR: Couldn't assign fermion type to "+name+".");
  }

  if(name=="up" || name == "down" || name == "electron" || name == "nu_el") {
    familyInd = 0;
  }
  else if (name=="charm" || name == "strange" || name == "muon" || name == "nu_mu") {
    familyInd = 1;
  }
  else if (name == "top" || name == "bottom" || name == "tauon" || name == "nu_tau") {
    familyInd = 2;
  }
  else {
    throw std::runtime_error("ERROR: Couldn't assign fermion family index to "+name+".");
  }

  if (type=="DOWN"){
    pdgID = 1+2*familyInd;
  } else if( type=="UP"){
    pdgID = 2+2*familyInd;
  } else if( type == "LEP") {
    pdgID = 11 + 2*familyInd;
  } else if(type == "NEU") {
    pdgID = 12 + 2*familyInd;
  }
  
    
  std::tuple<std::string,unsigned int,int>result=std::make_tuple(type,familyInd,pdgID);
  return result;
}



//Constructor: Setup model with configuration file specifying couplings
ZpModel::ZpModel(const char* configfile): bsm_parameters(0.1, 1500, 0) /*partial_widths(),*/  //Default values if no parameters are specified in config file
{
  std::fprintf(stderr,"\n*** CONSTRUCTING MODEL ***\n");
  
  //Initialize fermion list
  setup_flst();

  //Get model configuration from config file
  
  std::fprintf(stderr,"\nReading configuration file %s\n", configfile);
  conf_reader reader(configfile);
  dict config = reader.get_config();


  
  //Set GEBERAK model parameters
  dict::iterator it= config.find("GENERAL");
  if(it == config.end())
  {
    std::fprintf(stderr,"No model parameters specified. Using default:\n");
  }
  else
  {
    std::fprintf(stderr,"Setting model parameters:\n");
    set_mzp( (it->second)[0] );
    set_gx( (it->second)[1] );
    set_mixing( (it->second)[2] );
    set_whid( (it->second)[3] );
    // bare mixing mass is optional
    // Not yet implemented
    if((it->second).size()==5) {
      set_dm( (it->second)[4] );
    }
  }
  
  std::fprintf(stderr,
	      "\n"
	      "\t%-10s %-10s\n"
	      "\t%-10s|%-10g\n"
	      "\t%-10s|%-10g\n"
	      "\t%-10s|%-10g\n"
	      "\t%-10s|%-10g\n\n",
	      "Parameter", "Value",
	      "mzp", mzp_(),
	      "gx", gx_(),
	      "mixing", mixing_(),
	      "whid", whid_()
	      );
    
  //Applying FERMION CONFIGURATION:
  //Iterate over the whole fermion list and check for initialization values passed in config file
  std::fprintf(stderr,"Calculating vector couplings:\n\n");
  
  std::fprintf(stderr,"\t%-10s %-8s %-5s %-5s %-10s %-10s %-10s %-10s %-10s\n",
	      "Fermion", "mass", "cxl", "cxr", "qgam/e", "qzl", "qzr", "qzpl", "qzpr"
	      );
  for (fermion_list::iterator ferms=flst.begin(); ferms!=flst.end(); ++ferms)
  {
    auto fermType = fermionType(ferms->first);
    std::string type = std::get<0>(fermType);
    unsigned int famIndex = std::get<1>(fermType);
    int pdgID = std::get<2>(fermType);

    // Left handed coupling
    it = config.find((type+"L").c_str());  //Fermion label(string)
    
    if(it != config.end()) //No start parameters found
      {
	//Set new fermion parameters
	(ferms->second)->update_xlcharge( (it->second)[famIndex] );
      }
    
    // Right handed coupling
    it = config.find((type+"R").c_str());  //Fermion label(string)
    
    if(it != config.end()) //No start parameters found
      {
	//Set new fermion parameters
	(ferms->second)->update_xrcharge( (it->second)[famIndex] );
     }
    
    //Initialize fermion vector couplings
    (ferms->second)->set_vecc( new fundamental::vcoeff( *(ferms->second), *this) );   
    
    //Print fermion parameters after initialization
    std::fprintf(stderr,"\t%-10s|%-8g|%-5g|%-5g|%-10g|%-10g|%-10g|%-10g|%-10g\n",
		ferms->first.c_str(), (ferms->second)->m(),
		(ferms->second)->get_xlcharge(), (ferms->second)->get_xrcharge(),
		((ferms->second)->vecc()).q_gam/e_(),
		((ferms->second)->vecc()).q_zl, ((ferms->second)->vecc()).q_zr,
		((ferms->second)->vecc()).q_zpl, ((ferms->second)->vecc()).q_zpr
		);     
  } 

  // Check if couplings are universal
  // This is only a temporary check until implemented
  
  if(!reader.couplingUniversal()){
    throw std::runtime_error("ERROR: Non-universal couplings are not yet supported.");
  }
  
  //Initialize widths to -1 ("not yet calculated")
  partial_fwidths=NULL;
  higgs_width=-1.;
  wzp=-1.;
  wzp = wzp_(); // Initialize Z' width to its value in initial model
  std::fprintf(stderr,"\n*** MODEL CONSTRUCTED ***\n\n");
}





//DEFAULT CONSTRUCTOR: Generates Sequantial Standard Model 
// with gx=0.1 mzp=1500, mixing=0 and whid=0
///WARNING: Only use for benchmarks
ZpModel::ZpModel(double mzp): bsm_parameters(0.1, mzp, 0)
{
  //Initialize fermion list
  setup_flst();
  
  
  std::fprintf(stderr,"Sequential Standard Model couplings:\n\n");
  std::fprintf(stderr,"\t%-10s %-8s %-10s %-10s %-10s %-10s %-10s\n","Fermion", "mass", "qgam/e", "qzl", "qzr", "qzpl", "qzpr");
  for (fermion_list::iterator ferms=flst.begin(); ferms!=flst.end(); ++ferms)
  {
    //Initialize fermion vector couplings with default values
    
    (ferms->second)->set_vecc( new fundamental::vcoeff( *(ferms->second), *this) );   
    
//     (ferms->second)->change_mass(0); //#########################################################################3333 Only Testing
    //Set Zp couplings to Z couplings (SSM)
    (ferms->second)->set_qzpl((ferms->second)->vecc().q_zl);
    (ferms->second)->set_qzpr((ferms->second)->vecc().q_zr);
    
    //Print fermion parameters after initialization
    std::fprintf(stderr,"\t%-10s|%-8g|%-10g|%-10g|%-10g|%-10g|%-10g\n"
                ,ferms->first.c_str(), (ferms->second)->m(), ((ferms->second)->vecc()).q_gam/e_(), ((ferms->second)->vecc()).q_zl
                ,((ferms->second)->vecc()).q_zr, ((ferms->second)->vecc()).q_zpl,((ferms->second)->vecc()).q_zpr);     
  }
  std::fprintf(stderr,"\n");
  //Initialize widths to -1 ("not yet calculated")
  partial_fwidths=NULL;
  higgs_width=-1.;
  set_whid(0);
  wzp=-1.;
  wzp = wzp_(); // Initialize Z' width to its value in initial model
}





//Calculate partial fermionic width of Zp
double ZpModel::calc_fwidth(fundamental::fermionExt& f)
{
  double ratio = pow(f.get_mass()/mzp_(), 2);
//   std::cout<<f.get_mass();
  if(ratio>0.5*0.5)
  {
    std::fprintf(stderr,"INFO: Following decay channel kinematically not possible: 2*mf>mzp !\nSet to zero:\n");
    return 0;
  }else{
    return mzp_()* double(f.Nc()) /(24*M_PI) * sqrt(1.-4*ratio) * ( (pow(f.vecc().q_zpl, 2) + pow(f.vecc().q_zpr, 2))*(1.-ratio) + 6*f.vecc().q_zpl*f.vecc().q_zpr*ratio );
  }
}


//Calculate partial width  Zp-> Z H
double ZpModel::calc_zhwidth()
{
  if(mh_()+mz_()>mzp_())
  {
    std::fprintf(stderr,"INFO: Zp -> Z H kinematically not possible: mz + mh > mzp !\nSet to zero:\n");
    return 0;
  }else{
    double geff = (gz_()*sin(xi_()) + g1_()*cos(xi_())*tan(mixing_()))*(gz_()*cos(xi_()) - g1_()*sin(xi_())*tan(mixing_()));
    double pref = pow(geff*vev_(),2)/(96*M_PI*pow(mzp_(),3));
    double fac1 = 1 + pow( (mh_()*mh_() - mz_()*mz_() - mzp_()*mzp_())/(mz_()*mzp_()) ,2)/8;
    double fac2 = sqrt( pow(mzp_()*mzp_()+mz_()*mz_() -mh_()*mh_(), 2) - pow(2*mz_()*mzp_(),2) );
    return pref*fac1*fac2;
  }
}


//Calculate partial width  Zp-> W W
double ZpModel::calc_wwidth()
{
  if(2*mw_()>mzp_())
  {
    std::fprintf(stderr,"INFO: Zp -> W+ W- kinematically not possible: 2*mw > mzp !\nSet to zero:\n");
    return 0;
  }else{
    double ratio = pow(mw_()/mzp_(), 2);
    double pref = (1-sw2_())*pow(g2_()*sin(xi_()), 2.0)/(192*M_PI*mzp_());
    double matr = -48*pow(mw_(),2)-68*pow(mzp_(),2)+16*pow(mzp_(),2)/ratio+pow(mzp_()/ratio,2);
    return pref*sqrt(1-4*ratio)*matr;
  }
}



//Calculate total Zp width
void ZpModel::update_width()
{
    std::fprintf(stderr,"\nCalculating Zp width:\n\n\t%-14s %-14s\n", "Decay channel", "Width [GeV]");
    wzp=0;
    for(fermion_list::iterator it=flst.begin(); it!=flst.end(); ++it)
    {
      double pwidth = calc_fwidth( *(it->second) );
      wzp+= pwidth;
      std::fprintf(stderr,"\t%-14s|%-14g\n", it->first.c_str(), pwidth);
    }
    //Higgs width
    higgs_width=calc_zhwidth();
    wzp+=higgs_width;
    std::fprintf(stderr,"\t%-14s|%-14g\n", "Higgs Z", higgs_width);
    
    //WW width
    ww_width=calc_wwidth();
    wzp+=ww_width;
    wzp+=whid_();
    std::fprintf(stderr,"\t%-14s|%-14g\n", "WW", ww_width);
    std::fprintf(stderr,"\t%-14s|%-14g\n", "W_hidden", whid_());

    std::fprintf(stderr,"\t%-14s:%-14g\n\n","Total width",  wzp);
}




//Return Zp width
double ZpModel::wzp_()
{
  if(wzp==-1)  update_width(); //-1 means not yet calculated

  return wzp;  
}

