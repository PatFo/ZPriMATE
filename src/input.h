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
    bool _use_ssm;
    int _proc_id;    
    double _ebeam;
    double _mzssm;
    double _smax;
    double _smin;
    std::string _binning;
    std::string _efficiencies;
    std::string _limdir;
    std::string _model;
    std::string _pdfset;
  public:
    bool use_ssm();
    int proc_id();    
    double ebeam();
    double mzssm();
    double smax();
    double smin();
    std::string binning();
    std::string efficiencies();
    std::string limdir();
    std::string model();
    std::string pdfset();
    //Constructor
    settings(const char* startfile, bool verbose);
  };



#endif
