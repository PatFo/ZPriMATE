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
      float iso3;
      float mass;
      float emcharge;
      float xlcharge;
      float xrcharge;    
    public:
      //Get parameters
      int get_pdg();
      float get_iso3();
      float get_emcharge();
      float get_xlcharge();
      float get_xrcharge();
      float get_mass();
      //Set parameters
      void update_emcharge(float emc);
      void update_xlcharge(float xlc);
      void update_xrcharge(float xrc);
      void change_mass(float new_mass);
      //Constructor
      fermion(int n_pdg, float t3, float m, float emc, float xlc, float xrc);
  };
  
  
  //***************************************//
  //            MODEL BASE CLASS           //
  //***************************************//
  
  class bsm_parameters{
    private:
      //SM parameters from PDG (http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf)
      const static float aew=1./128;  // @q=MZ
      const static float as=0.1184;
      const static float mw=80.385;
      const static float mz=91.1876;
      const static float wz=2.4952;
      const static float sw2=0.23155;
      float e;
      float g1;
      float g2;
      float gz;
      float vev;
      //BSM parameters
      float gx;
      float mzp;
      float mixing_angle;
      float fef;
      float xi;
    public:
      //Get SM parameters
      float e_();
      float g1_();
      float g2_();
      float gz_();
      float vev_();
      float aew_();
      float as_();
      float mw_();
      float mz_();
      float wz_();
      float sw2_();
      //Get BSM parameters
      float gx_();
      float mzp_();
      float mixing_();
      float fef_();
      float xi_();
      //Set parameters
      void set_gx(float g);
      void set_mzp(float m);
      void set_mixing(float mix);
      //Constructor
      bsm_parameters(float cpl, float mass, float mix=0);
  };
  
  
  
  //*****************************************************//
  //          VECTOR COUPLING COEFFICIENT CLASS          //
  //*****************************************************//
  
  class vcoeff{
    private:
      float hypl;
      float hypr;
    public:
      float q_gam;
      float q_zl;
      float q_zr;
      float q_zpl;
      float q_zpr;
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
      float m();
      //Returns nc
      int Nc();
      //Vector coupling handling
      void set_vecc(vcoeff* ptr);
      vcoeff vecc();
      //Destuctor should take care of vcoeff pointer
      ~fermionExt();
    protected:
      //Constructor can only be used by derived classes
      fermionExt(bool massive, int fam, float t3, float m, float emc, float xlc, float xrc, int n);
  };
  
}
#endif