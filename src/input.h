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
    bool _verb;
    int _proc_id; 
    double _acc;
    double _ebeam;
    double _luminosity;
    double _mzssm;
    std::string _binning;
    std::string _efficiencies;
    std::string _limdir;
    std::string _model;
    std::string _odir;
    std::string _pdfset;
  public:
    bool use_ssm();
    bool verbose();
    int proc_id();    
    double int_acc();
    double ebeam();
    double luminosity();
    double mzssm();
    std::string binning();
    std::string efficiencies();
    std::string limdir();
    std::string model();
    std::string odir();
    std::string pdfset();
    //Constructor
    settings(const char* startfile);
  };



#endif
