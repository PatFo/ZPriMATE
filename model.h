#ifndef MODEL_H
#define MODEL_H


namespace fundamental{
  
  
  
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
  
  
}
#endif