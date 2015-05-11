#ifndef INPUT_H
#define INPUT_H

#include <map>
#include <string>
#include <vector>

typedef  std::map<std::string, double>  parmap;
typedef std::map<std::string, parmap >  dict;

typedef std::map<std::string, std::vector<std::string> > strmap;




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
    std::string _pdfset;
    std::string _binning;
    double _smin;
    double _smax;
    
  public:
    settings(const char* startfile, bool verbose);
  };



#endif
