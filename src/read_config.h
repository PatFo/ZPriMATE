#ifndef READ_CONFIG_H
#define READ_CONFIG_H

#include <map>
#include <string>


typedef  std::map<std::string, double>  parmap;
typedef std::map<std::string, parmap >  dict;




class conf_reader{
  private:
    dict config;
  
  public:
    dict get_config();
    //Constructor
    conf_reader(const char* filename);
  
};


#endif
