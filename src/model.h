#ifndef MODEL_H
#define MODEL_H


namespace fundamental{
  
  
  //**************************************//
  //         FERMION BASE CLASS           //
  //**************************************//
  
  class fermion{
    private:
      int family;
      int iso3;
      float mass;
      float emcharge;
      float xlcharge;
      float xrcharge;    
    public:
      //Get parameters
      int get_family();
      int get_iso3();
      float get_emcharge();
      float get_xlcharge();
      float get_xrcharge();
      float get_mass();
      //Set parameters
      void update_emcharge(float emc);
      void update_xlcharge(float xlc);
      void update_xrcharge(float xrc);
      void change_mass(float new_mass);    
    protected:
      //Constructor
      fermion(int fam, int t3, float m, float emc, float xlc, float xrc);
  };
  
  
  //***************************************//
  //            MODEL BASE CLASS           //
  //***************************************//
  
  class bsm_parameters{
    private:
      //SM parameters from PDG (http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf)
      const static float aew=1./128;
      const static float as=0.1184;
      const static float mz=91.1876;
      const static float wz=2.4952;
      const static float sw2=0.23155;
      //BSM parameters
      float gx;
      float mzp;
      float mixing_angle;
    public:
      //Get parameters
      float get_aew();
      float get_as();
      float get_mz();
      float get_wz();
      float get_sw2();
      float get_gx();
      float get_mzp();
      float get_mixing();
      //Set parameters
      void set_gx(float g);
      void set_mzp(float m);
      void set_mixing(float mix);
    protected:
      //Constructor
      bsm_parameters(float cpl, float mass, float mix=0);
  };
  
  
}
#endif