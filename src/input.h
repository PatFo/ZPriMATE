#ifndef INPUT_H
#define INPUT_H

#include <map>
#include <string>


typedef  std::map<std::string, double>  parmap;
typedef std::map<std::string, parmap >  dict;

typedef std::map<std::string, std::string> inmap;




  //**************************************//
  //         CONFIG READER                //  
  //**************************************//

class conf_reader{
  private:
    dict config;
  
  public:
    dict get_config();
    //Constructor
    conf_reader(const char* filename);
  
};





  //**************************************//
  //         START FILE READER            //
  //**************************************//
  
  
  class settings{
  private:
    
  public:
    settings(const char* startfile);
  };



#endif
