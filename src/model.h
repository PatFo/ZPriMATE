#ifndef MODEL_H
#define MODEL_H


namespace fundamental{
  
  
  //**************************************//
  //         FERMION BASE CLASS           //
  //**************************************//
  
  class fermion{
    private:
      //Fermion attributes
      int pdg;
      double iso3;
      double mass;
      double emcharge;
      double xlcharge;
      double xrcharge;    
    public:
      //Get parameters
      int get_pdg();
      double get_iso3();
      double get_emcharge();
      double get_xlcharge();
      double get_xrcharge();
      double get_mass();
      //Set parameters
      void update_emcharge(double emc);
      void update_xlcharge(double xlc);
      void update_xrcharge(double xrc);
      void change_mass(double new_mass);
      //Constructor
      fermion(int n_pdg, double t3, double m, double emc, double xlc, double xrc);
  };
  
  
  //***************************************//
  //            MODEL BASE CLASS           //
  //***************************************//
  
  class bsm_parameters{
    private:
      //SM parameters from PDG (http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf)
      const static double aew=1./127.9;  // @q=MZ
      const static double as=0.1184;
      const static double mw=80.385;
      const static double mz=91.1876;
      const static double wz=2.4952;
      const static double sw2=0.233699;// PDG=0.23155;
      double e;
      double g1;
      double g2;
      double gz;
      double vev;
      //BSM parameters
      double gx;
      double mzp;
      double mixing_angle;
      double fef;
      double xi;
      void calc_fef();
      void calc_xi();
      void update();
    public:
      //Get SM parameters
      double e_();
      double g1_();
      double g2_();
      double gz_();
      double vev_();
      double aew_();
      double as_();
      double mw_();
      double mz_();
      double wz_();
      double sw2_();
      //Get BSM parameters
      double gx_();
      double mzp_();
      double mixing_();
      double fef_();
      double xi_();
      //Set parameters
      void set_gx(double g);
      void set_mzp(double m);
      void set_mixing(double mix);
      //Constructor
      bsm_parameters(double cpl, double mass, double mix=0);
  };
  
  
  
  //*****************************************************//
  //          VECTOR COUPLING COEFFICIENT CLASS          //
  //*****************************************************//
  
  class vcoeff{
    private:
      double hypl;
      double hypr;
    public:
      double q_gam;
      double q_zl;
      double q_zr;
      double q_zpl;
      double q_zpr;
      double get_hypl();
      double get_hypr();
      //Constructor
      vcoeff(fermion f, bsm_parameters paras);    
  };
  
  
  
  
  //*************************************************************//
  //          FERMION CLASS EXTENDED BY VECTOR COUPLING          //
  //*************************************************************//
  
  
  
  class fermionExt: public fermion{
    private:
      //Pointer for storage of vector coupling object
      vcoeff* pvecc;
      //Number of color
      bool is_massive;
      int nc;
    public:
      //Gives the user defined mass
      double m();
      //Returns nc
      int Nc();
      //Vector coupling handling
      void set_vecc(vcoeff* ptr);
      vcoeff vecc();
      //Destuctor should take care of vcoeff pointer
      ~fermionExt();
    protected:
      //Constructor can only be used by derived classes
      fermionExt(bool massive, int fam, double t3, double m, double emc, double xlc, double xrc, int n);
  };
  
}
#endif